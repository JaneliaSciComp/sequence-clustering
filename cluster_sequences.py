import time
import argparse

def generate_neighbors(seq, include_next_nearest=False):
    """
    Generate all possible sequences that are 1 nucleotide different from the given sequence.
    Optionally include next-nearest neighbors (2 nucleotide differences).

    Args:
        seq (str): The input sequence.
        include_next_nearest (bool): If True, also generate 2-mismatch neighbors.

    Yields:
        str: A sequence that is 1 or 2 nucleotides different from the input sequence.
    """
    nucleotides = 'ATCG'

    # Generate 1-mismatch neighbors
    for i in range(len(seq)):
        for nucleotide in nucleotides:
            if nucleotide != seq[i]:
                yield seq[:i] + nucleotide + seq[i+1:]

    # Generate 2-mismatch neighbors if requested
    if include_next_nearest:
        for i in range(len(seq)):
            for j in range(i+1, len(seq)):
                for nuc1 in nucleotides:
                    if nuc1 != seq[i]:
                        for nuc2 in nucleotides:
                            if nuc2 != seq[j]:
                                neighbor = list(seq)
                                neighbor[i] = nuc1
                                neighbor[j] = nuc2
                                yield ''.join(neighbor)

def is_valid_sequence(seq):
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    Args:
        seq (str): The sequence to check.

    Returns:
        bool: True if sequence contains only A, G, C, T.
    """
    return all(c in 'ATCG' for c in seq)

def cluster_sequences(fastq_file, include_next_nearest=False):
    """
    Cluster sequences from FASTQ file using 1-mismatch clustering (optionally 2-mismatch).

    Args:
        fastq_file (str): Path to the FASTQ file.
        include_next_nearest (bool): If True, cluster with up to 2 mismatches.

    Returns:
        list: List of clusters, each containing sequences and their info.
    """
    # Dictionary mapping sequence -> cluster_id
    seq_to_cluster = {}
    # List of clusters, each cluster is a dict with sequences and their counts/headers
    clusters = []
    # Counter for generating unique cluster IDs
    next_cluster_id = 0
    # Statistics
    skipped_sequences = 0
    processed_sequences = 0

    with open(fastq_file, 'r') as f:
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

            # Check if this sequence is already assigned to a cluster
            if seq in seq_to_cluster:
                cluster_id = seq_to_cluster[seq]
                # Increment count for this sequence in the cluster
                if seq in clusters[cluster_id]['sequences']:
                    clusters[cluster_id]['sequences'][seq]['count'] += 1
                else:
                    clusters[cluster_id]['sequences'][seq] = {'count': 1, 'header': header}
                clusters[cluster_id]['total_count'] += 1
            else:
                # Check if any neighbor of this sequence is in a cluster
                found_cluster = None
                for neighbor in generate_neighbors(seq, include_next_nearest):
                    if neighbor in seq_to_cluster:
                        found_cluster = seq_to_cluster[neighbor]
                        break

                if found_cluster is not None:
                    # Add to existing cluster
                    cluster_id = found_cluster
                    clusters[cluster_id]['sequences'][seq] = {'count': 1, 'header': header}
                    clusters[cluster_id]['total_count'] += 1

                    # Map this sequence and all its neighbors to this cluster
                    seq_to_cluster[seq] = cluster_id
                    for neighbor in generate_neighbors(seq, include_next_nearest):
                        if neighbor not in seq_to_cluster:
                            seq_to_cluster[neighbor] = cluster_id
                else:
                    # Create new cluster
                    cluster_id = next_cluster_id
                    next_cluster_id += 1

                    clusters.append({
                        'sequences': {seq: {'count': 1, 'header': header}},
                        'total_count': 1
                    })

                    # Map this sequence and all its neighbors to this new cluster
                    seq_to_cluster[seq] = cluster_id
                    for neighbor in generate_neighbors(seq, include_next_nearest):
                        if neighbor not in seq_to_cluster:
                            seq_to_cluster[neighbor] = cluster_id

    print(f"Processed {processed_sequences} valid sequences")
    print(f"Skipped {skipped_sequences} sequences with invalid characters")

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

    with open(output_file, 'w') as f:
        for cluster in sorted_clusters:
            # Find the most abundant sequence in this cluster as representative
            rep_seq = max(cluster['sequences'].items(), key=lambda x: x[1]['count'])
            seq, seq_info = rep_seq

            header = seq_info['header']
            total_count = cluster['total_count']

            f.write(f">{header} N:{total_count}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Cluster FASTQ sequences with 1-mismatch tolerance and output representative sequences.'
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
        print("Including next-nearest neighbors (up to 2 mismatches)")
    else:
        print("Using nearest neighbors only (up to 1 mismatch)")

    clusters = cluster_sequences(args.input, args.include_next_nearest)

    print(f"Found {len(clusters)} clusters")
    print(f"Writing clustered sequences to: {args.output}")
    write_clustered_fasta(clusters, args.output)

    print("Done!")

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Total time taken: {end - start:.3g} seconds")
