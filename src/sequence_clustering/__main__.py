import time
import argparse
from dataclasses import dataclass
from collections import defaultdict

from polyleven import levenshtein


@dataclass
class UniqueSequence:
    """
    A unique sequence with the number of duplicate reads and header info of
    one of the reads.
    """
    sequence: str
    count: int
    header: str

class DisjointSetUnion:
    """
    Disjoint Set Union (Union-Find) data structure with path compression and union by rank.
    """
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        """Find the root of the set containing x with path compression."""
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int):
        """Union the sets containing a and b. Returns True if merged, False if
        already in the same set."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True

    def update(self, edges: list[tuple[int, int]], offset_a: int = 0, offset_b: int = 0):
        """Union all pairs in edges."""
        for u, v in edges:
            self.union(u + offset_a, v + offset_b)

    def get_components(self) -> list[list[int]]:
        """Get all components as lists of nodes."""
        components = defaultdict(list)
        for i in range(len(self.parent)):
            r = self.find(i)
            components[r].append(i)
        return list(components.values())


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

                    if levenshtein(seq_i, seq_j, n_edits) <= n_edits:
                        edges.append((i, j))

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
        for shift in range(max_shift + 1):
            # Bucket sequences by hashing the current partition
            hash_to_bucket = defaultdict(lambda: ([], []))
            for idx, seq in enumerate(sequences_a):
                sequence_partition = seq.sequence[start:end]
                hash_to_bucket[sequence_partition][0].append(idx)

            # Shift partition for sequences_b
            start += shift
            end += shift
            for idx, seq in enumerate(sequences_b):
                sequence_partition = seq.sequence[start:end]
                hash_to_bucket[sequence_partition][1].append(idx)

            # Compare sequences only within the same bucket
            buckets = list(hash_to_bucket.values())
            for bucket_a, bucket_b in buckets:
                if len(bucket_a) == 0 or len(bucket_b) == 0:
                    continue

                for i in bucket_a:
                    seq_i = sequences_a[i].sequence

                    for j in bucket_b:
                        seq_j = sequences_b[j].sequence

                        if levenshtein(seq_i, seq_j, n_edits) <= n_edits:
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
    lengths = sorted(sequences.keys())
    n_lengths = len(lengths)
    N_EDITS = 2 if args.include_next_nearest else 1
    # Flatten the list of lists into a single list
    all_sequences = [seq for l in lengths for seq in sequences[l]]

    for i in range(n_lengths):
        seqlen = lengths.pop(0)
        seqs = sequences[seqlen]
        total_counts = sum(s.count for s in sequences[seqlen])
        print(f"  Length {seqlen}: Found {len(seqs):,} unique sequences ({total_counts:,} total)")
        edges = connect_sequences_same_length(seqs, N_EDITS)
        dsu.update(edges, offset_a=offsets[i], offset_b=offsets[i])

        for j, other_len in enumerate(lengths):
            if abs(other_len - seqlen) > N_EDITS:
                continue
            other_seqs = sequences[other_len]
            edges = []
            edges = connect_sequences_different_length(seqs, other_seqs, N_EDITS)
            dsu.update(edges, offset_a=offsets[i], offset_b=offsets[j])

    clusters = []
    for components in dsu.get_components():
        representative_idx = max(components, key=lambda i: all_sequences[i].count)
        representative = all_sequences[representative_idx].sequence
        clusters.append(representative)

    cluster_time = time.time()
    print(f"Found {len(clusters):,} clusters in {cluster_time - read_time:.3g} seconds")

    # print(f"Writing clustered sequences to: {args.output}")
    # write_clustered_fasta(clusters, args.output)
    # write_time = time.time()
    # print(f"Writing took: {write_time - cluster_time:.3g} seconds")

    print("Done!")

    print(f"Total time taken: {time.time() - start_time:.3g} seconds")
