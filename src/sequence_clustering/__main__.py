import time
import argparse
from dataclasses import dataclass
from collections import defaultdict

from polyleven import levenshtein

from .dsu import DisjointSetUnion
from .io import FastQReader, FastAWriter
from .utils import generate_partitions, fill_buckets, compare_buckets


@dataclass
class UniqueSequence:
    """
    A unique sequence with the number of duplicate reads and header info of
    one of the reads.
    """
    sequence: str
    count: int
    header: str


def read_sequences(fastq_file_name: str) -> tuple[dict[int, list[UniqueSequence]], int]:
    """
    Read unique sequences from a FASTQ file.

    :param fastq_file_name: Path to the FASTQ file.
    :returns: A tuple containing the unique sequences separated by length, the
        total count of processed sequences, and the count of skipped invalid
        sequences.
    """
    sequences = defaultdict(dict)  # length -> {sequence: UniqueSequence}

    with FastQReader(fastq_file_name) as file_contents:
        for header, sequence in file_contents:
            sequences_for_length = sequences[len(sequence)]
            if sequence in sequences_for_length:
                sequences_for_length[sequence].count += 1
            else:
                sequences_for_length[sequence] = UniqueSequence(
                    sequence=sequence, count=1, header=header
                )

    sequences = {length: list(seq_dict.values()) for length, seq_dict in sequences.items()}

    return sequences, file_contents.total, file_contents.skipped


def write_clustered_fasta(clusters: list[UniqueSequence], output_file_name: str):
    """
    Write clustered sequences to FASTA file.

    Args:
        clusters (list): List of clusters from cluster_sequences function.
        output_file (str): Path to output FASTA file.
    """
    # Sort clusters by total count (descending)
    sorted_clusters = sorted(clusters, key=lambda c: c.count, reverse=True)

    with FastAWriter(output_file_name) as out_file:
        for cluster in sorted_clusters:
            out_file.write(cluster.header, cluster.sequence, cluster.count)


def connect_sequences_same_length(
    sequences: list[UniqueSequence], n_edits: int
) -> list[tuple[int, int]]:
    """
    Connect sequences of the same length based on edit distance of their representatives.
    
    :param sequences: List of UniqueSequence objects of the same length.
    :param n_edits: Number of allowed edits (typically 1 or 2).
    :returns: List of edges.
    """
    # Create n_edits + 1 partitions, so that at least one partition is guaranteed to match
    sequence_length = len(sequences[0].sequence)
    partitions = generate_partitions(sequence_length, n_edits + 1)

    edges = []

    for start, end in partitions:
        # Bucket sequences by hashing the current partition
        seed_to_bucket = fill_buckets(sequences, start, end)

        # Compare sequences only within the same bucket
        for bucket in seed_to_bucket.values():
            if len(bucket) < 2:
                continue

            compare_buckets(
                bucket, bucket,
                sequences, sequences,
                n_edits, edges
            )

    return edges


def connect_sequences_different_length(
    sequences_a: list[UniqueSequence], sequences_b: list[UniqueSequence], n_edits: int
) -> list[tuple[int, int]]:
    """
    Connect sequences of different length based on edit distance of their representatives.

    :param sequences_a: List of UniqueSequence objects.
    :param sequences_b: List of UniqueSequence objects of different length.
    :param n_edits: Number of allowed edits (typically 1 or 2).
    :returns: List of edges.
    """
    # Make sure sequences_a is the shorter list
    len_a = len(sequences_a[0].sequence)
    len_b = len(sequences_b[0].sequence)
    if len_a > len_b:
        edges = connect_sequences_different_length(sequences_b, sequences_a, n_edits)
        return [(j, i) for i, j in edges]

    # Create n_edits + 1 partitions, so that at least one partition is guaranteed to match
    partitions = generate_partitions(len_a, n_edits + 1)

    edges = []
    max_shift = min(n_edits, len_b - len_a)

    for start, end in partitions:
        # Bucket sequences by hashing the current partition
        seed_to_bucket_a = fill_buckets(sequences_a, start, end)

        # Shift partition for sequences_b
        seed_to_bucket_b = {}
        for shift in range(max_shift + 1):
            stb_for_shift = fill_buckets(sequences_b, start + shift, end + shift)
            seed_to_bucket_b.update(stb_for_shift)

        # Compare sequences only within the same bucket
        for seed, bucket_a in seed_to_bucket_a.items():
            bucket_b = seed_to_bucket_b.get(seed, [])
            if len(bucket_a) == 0 or len(bucket_b) == 0:
                continue

            compare_buckets(
                bucket_a, bucket_b,
                sequences_a, sequences_b,
                n_edits, edges
            )

    return edges


def main():
    """Main function to read, cluster, and write sequences."""
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
    sequences = dict(sorted(sequences.items()))
    read_time = time.time()
    n_sequences = [len(v) for v in sequences.values()]
    total_unique = sum(n_sequences)
    print(
        f"Read {total_unique:,} unique sequences ({total:,} total, {skipped:,} skipped) "
        f"in {read_time - start_time:.3g} seconds"
    )

    # Compute offsets of sequences into a global index space
    offsets = [0]
    for n in n_sequences:
        offsets.append(offsets[-1] + n)

    # Compute edges and cluster sequences based on those
    dsu = DisjointSetUnion(total_unique)
    lengths = list(sequences.keys())
    n_edits = 2 if args.include_next_nearest else 1

    for i in range(len(sequences)):
        seqlen = lengths.pop(0)
        seqs = sequences[seqlen]
        total_counts = sum(s.count for s in sequences[seqlen])
        print(f"  Length {seqlen}: Found {len(seqs):,} unique sequences ({total_counts:,} total)")
        edges = connect_sequences_same_length(seqs, n_edits)
        print(f"    Found {len(edges):,} links among length {seqlen} sequences")
        dsu.update(edges, offset_a=offsets[i], offset_b=offsets[i])

        for j, other_len in enumerate(lengths, start=i + 1):
            if abs(other_len - seqlen) > n_edits:
                continue
            other_seqs = sequences[other_len]
            edges = connect_sequences_different_length(seqs, other_seqs, n_edits)
            print(f"    Found {len(edges):,} links with length {other_len} sequences")
            dsu.update(edges, offset_a=offsets[i], offset_b=offsets[j])

    # Flatten the list of lists into a single list and report cluster representatives
    clusters = []
    sequences = [seq for seqs in sequences.values() for seq in seqs]
    for components in dsu.get_components():
        representative_idx = max(components, key=lambda i: sequences[i].count)
        representative = sequences[representative_idx]
        representative.count = sum(sequences[i].count for i in components)
        clusters.append(representative)

    cluster_time = time.time()
    print(f"Found {len(clusters):,} clusters in {cluster_time - read_time:.3g} seconds")

    print(f"Writing clustered sequences to: {args.output}")
    write_clustered_fasta(clusters, args.output)
    write_time = time.time()
    print(f"Writing took: {write_time - cluster_time:.3g} seconds")

    print("Done!")

    print(f"Total time taken: {time.time() - start_time:.3g} seconds")

if __name__ == '__main__':
    main()
