import time
import argparse
import numpy as np
import faiss

# ---------- Encoding (exact, 4-bit one-hot per base) ----------
_BASE2IDX = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    _BASE2IDX[b] = i

def encode_batch_onehot_4bit(sequences: list) -> dict:
    """Encode multiple sequences into binary codes, grouped by length.
    
    Args:
        sequences: List of sequence strings
        
    Returns:
        dict: {length: (codes_array, indices)} where codes_array is (N, d_bytes) 
              and indices maps back to original sequence positions
    """
    # Group sequences by length
    length_groups = {}
    for seq_idx, seq in enumerate(sequences):
        s = seq.strip().upper()
        L = len(s)
        if L not in length_groups:
            length_groups[L] = {'sequences': [], 'indices': []}
        length_groups[L]['sequences'].append(s)
        length_groups[L]['indices'].append(seq_idx)
    
    # Encode each length group
    result = {}
    for L, group in length_groups.items():
        if L == 0:
            continue
            
        seqs = group['sequences']
        n_seqs = len(seqs)
        
        # Calculate dimensions
        nbits = 4 * L
        d_bits = ((nbits + 7) // 8) * 8  # Round up to multiple of 8
        d_bytes = d_bits // 8
        
        # Pre-allocate output array
        codes = np.zeros((n_seqs, d_bytes), dtype=np.uint8)
        
        # Convert all sequences to byte arrays
        seq_bytes = np.array([list(seq.encode('ascii')) for seq in seqs], dtype=np.uint8)
        
        # Vectorized base to index conversion
        indices = _BASE2IDX[seq_bytes]  # Shape: (n_seqs, L)
        
        # Check for invalid bases
        invalid_mask = indices < 0
        if np.any(invalid_mask):
            invalid_seqs = np.any(invalid_mask, axis=1)
            # Remove invalid sequences
            valid_mask = ~invalid_seqs
            indices = indices[valid_mask]
            group['indices'] = [group['indices'][i] for i in range(n_seqs) if valid_mask[i]]
            codes = codes[valid_mask]
            n_seqs = np.sum(valid_mask)
            
        if n_seqs == 0:
            continue
            
        # Vectorized bit position calculation
        # For each sequence, calculate bit indices: 4*position + base_index
        pos_array = np.arange(L, dtype=np.int64)  # [0, 1, 2, ..., L-1]
        bit_indices = 4 * pos_array[None, :] + indices.astype(np.int64)  # Shape: (n_seqs, L)
        
        # Calculate byte and bit positions
        byte_pos = bit_indices >> 3  # Divide by 8
        bit_pos = bit_indices & 7    # Modulo 8
        
        # Vectorized bit setting
        seq_indices = np.arange(n_seqs)[:, None]  # Shape: (n_seqs, 1)
        codes[seq_indices, byte_pos] |= (1 << bit_pos).astype(np.uint8)
        
        result[L] = {
            'codes': codes,
            'indices': group['indices'],
            'd_bits': d_bits
        }
    
    return result

def cluster_sequences_vectorized(sequences: list, headers: list, include_next_nearest=False):
    """Vectorized clustering of sequences using FAISS.
    
    Args:
        sequences: List of sequence strings
        headers: List of corresponding headers
        include_next_nearest: If True, cluster with up to 2 mismatches
        
    Returns:
        dict: cluster_id -> cluster_info
    """
    max_hamming = 4 if include_next_nearest else 2
    hamming_check = hamming_leq2_nt if include_next_nearest else hamming_leq1_nt
    
    # Encode all sequences grouped by length
    encoded_groups = encode_batch_onehot_4bit(sequences)
    
    cluster_data = {}
    cluster_id_counter = 0
    
    # Process each length group
    for L, group_data in encoded_groups.items():
        codes = group_data['codes']
        seq_indices = group_data['indices']
        d_bits = group_data['d_bits']
        n_seqs = len(codes)
        
        if n_seqs == 0:
            continue
            
        print(f"Processing {n_seqs} sequences of length {L}")
        
        # Create FAISS index - use MIH for better performance on radius searches
        try:
            # IndexBinaryMultiHash(d, nhash, bits_per_hash) where:
            # d = dimension in bits
            # nhash = number of hash tables (typically 4-8)
            # bits_per_hash = bits per hash table (typically d/nhash or close to it)
            nhash = min(8, max(4, d_bits // 32))  # Adaptive number of hash tables
            bits_per_hash = d_bits // nhash  # Bits per hash table
            index = faiss.IndexBinaryMultiHash(d_bits, nhash, bits_per_hash)
            print(f"  Using MIH with {nhash} hash tables, {bits_per_hash} bits each")
        except (RuntimeError, ValueError) as e:
            # Fall back to flat index if MIH fails
            print(f"  MIH failed ({e}), using flat index")
            index = faiss.IndexBinaryFlat(d_bits)
        
        # Track which sequences are representatives vs clustered
        cluster_assignments = np.full(n_seqs, -1, dtype=np.int32)  # -1 means unassigned
        representatives = []  # List of (rep_seq_idx, cluster_id)
        
        for seq_idx in range(n_seqs):
            if cluster_assignments[seq_idx] != -1:
                continue  # Already assigned
                
            orig_seq_idx = seq_indices[seq_idx]
            seq = sequences[orig_seq_idx]
            header = headers[orig_seq_idx]
            
            if index.ntotal == 0:
                # First sequence of this length - becomes a representative
                index.add(codes[seq_idx:seq_idx+1])
                cluster_id = cluster_id_counter
                cluster_id_counter += 1
                
                representatives.append((seq_idx, cluster_id))
                cluster_assignments[seq_idx] = cluster_id
                
                cluster_data[cluster_id] = {
                    'sequences': {seq: {'count': 1, 'header': header}},
                    'total_count': 1
                }
            else:
                # Search for nearest neighbor
                D, I = index.search(codes[seq_idx:seq_idx+1], k=1)
                d0 = int(D[0, 0])
                i0 = int(I[0, 0])
                
                if i0 != -1 and d0 <= max_hamming:
                    # Get the representative sequence index and cluster ID
                    rep_seq_idx, cluster_id = representatives[i0]
                    rep_orig_idx = seq_indices[rep_seq_idx]
                    rep_seq = sequences[rep_orig_idx]
                    
                    # Verify with nucleotide-level Hamming distance
                    if hamming_check(seq, rep_seq):
                        # Assign to existing cluster
                        cluster_assignments[seq_idx] = cluster_id
                        
                        if seq in cluster_data[cluster_id]['sequences']:
                            cluster_data[cluster_id]['sequences'][seq]['count'] += 1
                        else:
                            cluster_data[cluster_id]['sequences'][seq] = {'count': 1, 'header': header}
                        cluster_data[cluster_id]['total_count'] += 1
                        continue
                
                # Create new cluster (new representative)
                index.add(codes[seq_idx:seq_idx+1])
                cluster_id = cluster_id_counter
                cluster_id_counter += 1
                
                representatives.append((seq_idx, cluster_id))
                cluster_assignments[seq_idx] = cluster_id
                
                cluster_data[cluster_id] = {
                    'sequences': {seq: {'count': 1, 'header': header}},
                    'total_count': 1
                }
    
    return cluster_data

# ---------- Exact Hamming check on nucleotides (fast, short-circuit) ----------
def hamming_leq2_nt(seq_a: str, seq_b: str) -> bool:
    if len(seq_a) != len(seq_b):
        return False
    diff = 0
    for ca, cb in zip(seq_a, seq_b):
        if ca != cb:
            diff += 1
            if diff > 2:
                return False
    return True

def hamming_leq1_nt(seq_a: str, seq_b: str) -> bool:
    if len(seq_a) != len(seq_b):
        return False
    diff = 0
    for ca, cb in zip(seq_a, seq_b):
        if ca != cb:
            diff += 1
            if diff > 1:
                return False
    return True

def is_valid_sequence(seq):
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    Args:
        seq (str): The sequence to check.

    Returns:
        bool: True if sequence contains only A, G, C, T.
    """
    return all(c in 'ACGT' for c in seq)

def cluster_sequences(fastq_file, include_next_nearest=False):
    """
    Vectorized FAISS-based clustering of sequences from FASTQ file.

    Args:
        fastq_file (str): Path to the FASTQ file.
        include_next_nearest (bool): If True, cluster with up to 2 mismatches.

    Returns:
        list: List of clusters, each containing sequences and their info.
    """
    # Load all sequences and headers into memory first
    sequences = []
    headers = []
    skipped_sequences = 0
    processed_sequences = 0

    print("Loading sequences...")
    with open(fastq_file, 'r', encoding='utf-8') as f:
        while True:
            header_line = f.readline().strip()
            if not header_line:
                break

            # Remove the '@' from the header
            header = header_line[1:] if header_line.startswith('@') else header_line
            seq = f.readline().strip()
            f.readline()  # Skip '+' line
            f.readline()  # Skip quality line

            # Skip sequences with invalid characters
            if not is_valid_sequence(seq):
                skipped_sequences += 1
                continue

            sequences.append(seq)
            headers.append(header)
            processed_sequences += 1

            # Progress reporting during loading
            if processed_sequences % 100000 == 0:
                print(f"Loaded {processed_sequences:,} sequences...")

    print(f"Loaded {processed_sequences:,} valid sequences")
    print(f"Skipped {skipped_sequences:,} sequences with invalid characters")
    
    if not sequences:
        return []

    print("Starting vectorized clustering...")
    cluster_data = cluster_sequences_vectorized(sequences, headers, include_next_nearest)

    # Convert to list format expected by write function
    clusters = list(cluster_data.values())
    return clusters

def write_clustered_fasta(clusters, output_file):
    """
    Write clustered sequences to FASTA file.

    Args:
        clusters (list): List of clusters from cluster_sequences function.
        output_file (str): Path to output FASTA file.
    """
    # Sort clusters by total count (descending)
    sorted_clusters = sorted(clusters, key=lambda x: x['total_count'], reverse=True)

    with open(output_file, 'w', encoding='utf-8') as f:
        for cluster in sorted_clusters:
            # Find the most abundant sequence in this cluster as representative
            rep_seq = max(cluster['sequences'].items(), key=lambda x: x[1]['count'])
            seq, seq_info = rep_seq

            header = seq_info['header']
            total_count = cluster['total_count']

            f.write(f">{header} N:{total_count}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description='FAISS-based clustering of FASTQ sequences.'
    )
    parser.add_argument('-i', '--input', default='merged_lenfilt.fq',
                       help='Path to input FASTQ file (default: merged_lenfilt.fq)')
    parser.add_argument('-o', '--output', default='merged_lenfilt_cluster.fasta',
                       help='Path to output FASTA file (default: merged_lenfilt_cluster.fasta)')
    parser.add_argument('--include_next_nearest', action='store_true',
                       help='Include next-nearest neighbors (2 mismatches) in clustering')

    args = parser.parse_args()

    print(f"Reading FASTQ file: {args.input}")
    if args.include_next_nearest:
        print("Using FAISS with next-nearest neighbors (up to 2 mismatches)")
    else:
        print("Using FAISS with nearest neighbors only (up to 1 mismatch)")

    clusters = cluster_sequences(args.input, args.include_next_nearest)

    print(f"Found {len(clusters):,} clusters")
    print(f"Writing clustered sequences to: {args.output}")
    write_clustered_fasta(clusters, args.output)

    print("Done!")

if __name__ == "__main__":
    start = time.time()
    main()
    print(f"Total time taken: {time.time() - start:.3g} seconds")
