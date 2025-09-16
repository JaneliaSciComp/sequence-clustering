import time
import argparse
from dataclasses import dataclass

import numpy as np


INVALID_NUCLEOTIDE = 255  # Using uint8, so 255 is a convenient invalid value

@dataclass
class UniqueSequence:
    sequence: str
    count: int
    header: str

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


def is_valid_sequence(seq: str) -> bool:
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    :param seq: The sequence to check.
    :returns: True if sequence contains only A, G, C, T.
    """
    return all(c in 'ATCG' for c in seq)


def read_sequences(fastq_file: str) -> tuple[dict[int, list[UniqueSequence]], int]:
    """
    Generator to read sequences from a FASTQ file.

    :param fastq_file: Path to the FASTQ file.
    :returns: A tuple containing the unique sequences separated by length, the
        total count of processed sequences, and the count of skipped invalid
        sequences.
    """
    sequences = {}
    total = 0
    skipped = 0

    with open(fastq_file, 'r', encoding='ascii') as f:
        while True:
            header_line = f.readline().strip()
            if not header_line:
                break
            seq = f.readline().strip()

            f.readline()  # Skip '+' line
            f.readline()  # Skip quality line

            if not is_valid_sequence(seq):
                skipped += 1
                continue
            total += 1

            seqlen = len(seq)
            if seqlen not in sequences:
                sequences[seqlen] = {}
            sequences_for_length = sequences[seqlen]

            if seq in sequences_for_length:
                sequences_for_length[seq].count += 1
            else:
                sequences_for_length[seq] = UniqueSequence(
                    sequence=seq,
                    count=1,
                    header=header_line[1:] if header_line.startswith('@') else header_line
                )

    sequences = {length: list(seq_dict.values()) for length, seq_dict in sequences.items()}

    return sequences, total, skipped


def sequences_to_numpy(sequences: list[UniqueSequence]) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert list of UniqueSequence to a numpy array representation.

    :param sequences: List of UniqueSequence objects.
    :returns: A tuple containing a 2D numpy array of sequences and a 1D array of counts.
    """
    sequences = [
        np.frombuffer(s.sequence.encode("ascii", "strict"), dtype="uint8")
        for s in sequences
    ]
    counts = np.array([s.count for s in sequences], dtype='int64')

    return np.array(sequences, dtype='uint8'), counts


def cluster_sequences(fastq_file: str, include_next_nearest: bool = False) -> list[dict]:
    """
    Memory-optimized clustering of sequences from FASTQ file.

    :param fastq_file: Path to the FASTQ file.
    :param include_next_nearest: If True, cluster with up to 2 mismatches.
    :returns: List of clusters, each containing sequences and their info.
    """
    # Only store actual sequences seen, not all possible neighbors
    seq_to_cluster = {}
    # List of actual sequences in each cluster (for neighbor checking)
    cluster_sequences = []
    # Cluster metadata
    clusters = []
    # Counter for generating unique cluster IDs
    next_cluster_id = 0
    # Statistics
    skipped_sequences = 0
    processed_sequences = 0

    def find_matching_cluster(seq):
        """Find if sequence matches any existing cluster within distance threshold."""
        for neighbor in generate_neighbors(seq, include_next_nearest):
            if neighbor in seq_to_cluster:
                return seq_to_cluster[neighbor]
        return None

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

            # Progress reporting
            if processed_sequences % 1000000 == 0:
                print(f"Processed {processed_sequences:,} sequences, {len(clusters)} clusters")

            # Check if this exact sequence is already assigned to a cluster
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
                found_cluster = find_matching_cluster(seq)

                if found_cluster is not None:
                    # Add to existing cluster
                    cluster_id = found_cluster
                    clusters[cluster_id]['sequences'][seq] = {'count': 1, 'header': header}
                    clusters[cluster_id]['total_count'] += 1

                    # Only map the actual sequence we've seen
                    seq_to_cluster[seq] = cluster_id
                    cluster_sequences[cluster_id].add(seq)
                else:
                    # Create new cluster
                    cluster_id = next_cluster_id
                    next_cluster_id += 1

                    clusters.append({
                        'sequences': {seq: {'count': 1, 'header': header}},
                        'total_count': 1
                    })

                    # Only map the actual sequence we've seen
                    seq_to_cluster[seq] = cluster_id
                    cluster_sequences.append({seq})

    print(f"Processed {processed_sequences:,} valid sequences")
    print(f"Skipped {skipped_sequences:,} sequences with invalid characters")

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Memory-optimized clustering of FASTQ sequences.'
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

    start = time.time()
    sequences, total, skipped = read_sequences(args.input)
    read_time = time.time()
    total_unique = sum(len(v) for v in sequences.values())
    print(
        f"Read {total_unique:,} unique sequences ({total:,} total, {skipped:,} skipped) "
        f"in {read_time - start:.3g} seconds"
    )

    for seqlen in sorted(sequences.keys()):
        seqs = sequences[seqlen]
        sequences[seqlen] = sorted(seqs, key=lambda x: x.count, reverse=True)
        total_counts = sum(s.count for s in sequences[seqlen])
        print(f"Length {seqlen}: Found {len(seqs):,} unique sequences ({total_counts:,} total)")
    dedup_time = time.time()
    print(f"Sorted unique sequences by frequency in {dedup_time - read_time:.3g} seconds")

    # clusters = cluster_sequences(args.input, args.include_next_nearest)
    # cluster_time = time.time()
    # print(f"Clustering took: {cluster_time - read_time:.3g} seconds")

    # print(f"Found {len(clusters):,} clusters")
    # print(f"Writing clustered sequences to: {args.output}")
    # write_clustered_fasta(clusters, args.output)
    # write_time = time.time()
    # print(f"Writing took: {write_time - cluster_time:.3g} seconds")

    print("Done!")

    print(f"Total time taken: {time.time() - start:.3g} seconds")
