import time
import argparse
import gc
import numpy as np
import faiss

# ---------- Encoding (exact, 4-bit one-hot per base) ----------
_BASE2IDX = np.full(256, -1, dtype=np.int8)
for i, b in enumerate(b"ACGT"):
    _BASE2IDX[b] = i

def encode_onehot_4bit(seq: str) -> np.ndarray:
    """Return packed binary code (uint8 array) with 4 bits per nucleotide.
       Exactly 1 bit set per position; length = 4*L bits, padded to multiple of 8."""
    s = np.frombuffer(seq.strip().upper().encode("ascii"), dtype=np.uint8)
    idx = _BASE2IDX[s]
    if np.any(idx < 0):
        raise ValueError("Only A/C/G/T are allowed")
    L = len(idx)
    nbits = 4 * L
    # Pad to multiple of 8 bits for FAISS binary index compatibility
    nbits_padded = ((nbits + 7) // 8) * 8
    nbytes = nbits_padded // 8
    code = np.zeros(nbytes, dtype=np.uint8)

    # For position p with base index k in {0,1,2,3}, set bit at bidx = 4*p + k
    bidx = 4 * np.arange(L, dtype=np.int64) + idx.astype(np.int64)

    byte_pos = bidx >> 3
    bit_pos  = bidx & 7
    code[byte_pos] |= (1 << bit_pos).astype(np.uint8)
    return code

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

# ---------- FAISS-based clustering for pipeline integration ----------
class FaissClustering:
    """
    FAISS-based clustering that maintains cluster information compatible with the pipeline.
    """
    def __init__(self, include_next_nearest=False):
        self.include_next_nearest = include_next_nearest
        self.max_hamming = 4 if include_next_nearest else 2  # 2 nt diff = 4 bit diff in one-hot
        self._perL = {}   # L -> dict(index=faiss.IndexBinary*, sequences=[dict])
        
    def _get_index_for_L(self, L: int):
        entry = self._perL.get(L)
        if entry is not None:
            return entry
        # Calculate bits and ensure it's a multiple of 8 for FAISS binary index
        nbits = 4 * L
        d_bits = ((nbits + 7) // 8) * 8  # Round up to nearest multiple of 8
        # Use flat binary index for simplicity and reliability
        index = faiss.IndexBinaryFlat(d_bits)
        entry = {"index": index, "sequences": []}
        self._perL[L] = entry
        return entry

    def process_sequence(self, seq: str, header: str):
        """
        Process one sequence and assign it to a cluster.
        Returns the cluster index or creates a new one.
        """
        s = seq.strip().upper()
        L = len(s)
        if L == 0:
            return None
            
        entry = self._get_index_for_L(L)
        index = entry["index"]
        sequences = entry["sequences"]

        code = encode_onehot_4bit(s)[None, :]  # shape (1, d_bits/8)

        if index.ntotal > 0:
            # Find nearest neighbor
            D, I = index.search(code, k=1)
            d0 = int(D[0, 0])
            i0 = int(I[0, 0])
            
            if i0 != -1 and d0 <= self.max_hamming:
                # Verify on nucleotides to guard against rare bit collisions
                rep_seq = sequences[i0]['representative']
                hamming_check = hamming_leq2_nt if self.include_next_nearest else hamming_leq1_nt
                
                if hamming_check(s, rep_seq):
                    # Add to existing cluster
                    cluster_idx = (L, i0)  # Use (length, index) as cluster identifier
                    return cluster_idx

        # Create new cluster (new representative)
        new_idx = len(sequences)
        sequences.append({
            'representative': s,
            'sequences': {},
            'total_count': 0
        })
        index.add(code)
        
        cluster_idx = (L, new_idx)
        return cluster_idx

def cluster_sequences(fastq_file, include_next_nearest=False):
    """
    FAISS-based clustering of sequences from FASTQ file.

    Args:
        fastq_file (str): Path to the FASTQ file.
        include_next_nearest (bool): If True, cluster with up to 2 mismatches.

    Returns:
        list: List of clusters, each containing sequences and their info.
    """
    clusterer = FaissClustering(include_next_nearest)
    cluster_data = {}  # cluster_idx -> cluster info
    
    # Statistics
    skipped_sequences = 0
    processed_sequences = 0

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

            processed_sequences += 1

            # Progress reporting
            if processed_sequences % 100000 == 0:
                print(f"Processed {processed_sequences:,} sequences, {len(cluster_data)} clusters")
                # Force garbage collection periodically
                gc.collect()

            try:
                cluster_idx = clusterer.process_sequence(seq, header)
                if cluster_idx is None:
                    continue
                    
                # Initialize cluster if new
                if cluster_idx not in cluster_data:
                    cluster_data[cluster_idx] = {
                        'sequences': {},
                        'total_count': 0
                    }
                
                # Add sequence to cluster
                if seq in cluster_data[cluster_idx]['sequences']:
                    cluster_data[cluster_idx]['sequences'][seq]['count'] += 1
                else:
                    cluster_data[cluster_idx]['sequences'][seq] = {'count': 1, 'header': header}
                cluster_data[cluster_idx]['total_count'] += 1
                
            except ValueError:
                # Skip sequences that can't be encoded (shouldn't happen with valid sequences)
                skipped_sequences += 1

    print(f"Processed {processed_sequences:,} valid sequences")
    print(f"Skipped {skipped_sequences:,} sequences with invalid characters")

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
