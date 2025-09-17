import time
import argparse
from dataclasses import dataclass
from collections import defaultdict

from sequence_clustering import check_levenshtein_distance_same_length


@dataclass
class UniqueSequence:
    """
    A unique sequence with the number of duplicate reads and header info of
    one of the reads.
    """
    sequence: str
    count: int
    header: str


def is_valid_sequence(seq: str) -> bool:
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    :param seq: The sequence to check.
    :returns: True if sequence contains only A, G, C, T.
    """
    return all(c in 'ATCG' for c in seq)


def read_sequences(fastq_file: str) -> tuple[dict[int, list[UniqueSequence]], int]:
    """
    Read unique sequences from a FASTQ file.

    :param fastq_file: Path to the FASTQ file.
    :returns: A tuple containing the unique sequences separated by length, the
        total count of processed sequences, and the count of skipped invalid
        sequences.
    """
    sequences = defaultdict(dict)  # length -> {sequence: UniqueSequence}
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

            sequences_for_length = sequences[len(seq)]
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


def generate_partitions(n: int, k: int) -> list[tuple[int, int]]:
    """
    Partition n items into k contiguous sets of roughly equal size.

    :param n: Total number of items.
    :param k: Number of subsets.
    :returns: A list of all partitions.
    """
    base_size = n // k
    remainder = n % k

    partitions = []
    start = 0
    for i in range(k):
        size = base_size + (1 if i < remainder else 0)
        end = start + size
        partitions.append((start, end))
        start = end

    return partitions


def connect_sequences_same_length(sequences: list[UniqueSequence], n_edits: int) -> list[tuple[int, int]]:
    """
    Cluster sequences of the same length based on edit distance of their representatives.
    
    :param sequences: List of UniqueSequence objects of the same length.
    :param n_edits: Number of allowed edits (typically 1 or 2).
    :returns: List of clusters.
    """
    # Create n_edits + 1 partitions, so that at least one partition is guaranteed to match
    sequence_length = len(sequences[0].sequence)
    partitions = generate_partitions(sequence_length, n_edits + 1)

    edges = []

    for start, end in partitions:
        # Bucket sequences by hashing the current partition
        hash_to_bucket = defaultdict(list)
        for idx, seq in enumerate(sequences):
            sequence_partition = seq.sequence[start:end]
            hash_to_bucket[sequence_partition].append(idx)

        # Compare sequences only within the same bucket
        buckets = list(hash_to_bucket.values())
        for bucket in buckets:
            if len(bucket) < 2:
                continue

            while bucket:
                i = bucket.pop()
                seq_i = sequences[i].sequence

                for j in bucket:
                    seq_j = sequences[j].sequence

                    if check_levenshtein_distance_same_length(seq_i, seq_j, n_edits):
                        edges.append((i, j))

    return edges


def write_clustered_fasta(clusters, output_file):
    """
    Write clustered sequences to FASTA file.

    Args:
        clusters (list): List of clusters from cluster_sequences function.
        output_file (str): Path to output FASTA file.
    """
    # Sort clusters by total count (descending)
    sorted_clusters = sorted(clusters, key=lambda x: x['total_count'], reverse=True)

    with open(output_file, 'w', encoding='ascii') as f:
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

    start_time = time.time()
    sequences, total, skipped = read_sequences(args.input)
    read_time = time.time()
    total_unique = sum(len(v) for v in sequences.values())
    print(
        f"Read {total_unique:,} unique sequences ({total:,} total, {skipped:,} skipped) "
        f"in {read_time - start_time:.3g} seconds"
    )

    edges = {}
    for seqlen in sorted(sequences.keys()):
        seqs = sequences[seqlen]
        total_counts = sum(s.count for s in sequences[seqlen])
        print(f"  Length {seqlen}: Found {len(seqs):,} unique sequences ({total_counts:,} total)")
        N_EDITS = 2 if args.include_next_nearest else 1
        edges[seqlen] = connect_sequences_same_length(seqs, N_EDITS)
    cluster_time = time.time()
    total_edges = sum(len(e) for e in edges.values())
    print(f"Found {total_edges:,} edges in {cluster_time - read_time:.3g} seconds")

    # print(f"Writing clustered sequences to: {args.output}")
    # write_clustered_fasta(clusters, args.output)
    # write_time = time.time()
    # print(f"Writing took: {write_time - cluster_time:.3g} seconds")

    print("Done!")

    print(f"Total time taken: {time.time() - start_time:.3g} seconds")
