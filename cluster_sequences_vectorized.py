import time
import argparse
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class UniqueSequence:
    """
    A unique sequence with the number of duplicate reads and header info of
    one of the reads.
    """
    sequence: str
    count: int
    header: str


@dataclass
class Cluster:
    """
    A cluster of similar sequences.
    """
    representative: UniqueSequence
    members: list[UniqueSequence]

    @property
    def total_count(self) -> int:
        """Return the total count of sequences in the cluster."""
        return sum(member.count for member in self.members) + self.representative.count

    def absorb(self, other: 'Cluster'):
        """Absorb another cluster into this one."""
        self.members.append(other.representative)
        self.members.extend(other.members)

        other.members.clear()
        other.representative = None  # type: ignore

    def __len__(self) -> int:
        return 1 + len(self.members)

    def __iter__(self):
        yield self.representative
        yield from self.members


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


def are_sequences_similar(seq1: str, seq2: str, n_edits: int) -> bool:
    """
    Check if two sequences are similar within a given number of edits (Hamming distance).

    :param seq1: First sequence.
    :param seq2: Second sequence.
    :param n_edits: Maximum number of allowed edits.
    :returns: True if sequences are similar within n_edits, False otherwise.
    """
    differences = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    return differences <= n_edits


def cluster_sequences_same_length(sequences: list[UniqueSequence], n_edits: int) -> list[Cluster]:
    # Create initial single-member clusters (start with most abundant sequences)
    sequences = sorted(seqs, key=lambda x: x.count, reverse=True)
    clusters = [Cluster(representative=seq, members=[]) for seq in sequences]

    # Create n_edits + 1 partitions, so that at least one partition is guaranteed to match
    sequence_length = len(clusters[0].representative.sequence)
    partitions = generate_partitions(sequence_length, n_edits + 1)

    for start, end in partitions:
        hash_to_candidates = defaultdict(list)
        for cluster in clusters:
            sequence_partition = cluster.representative.sequence[start:end]
            hash_to_candidates[sequence_partition].append(cluster)

        candidate_lists = list(hash_to_candidates.values())
        for candidates in candidate_lists:
            if len(candidates) < 2:
                continue

            for i, candidate_i in enumerate(candidates):
                if candidate_i.representative is None:
                    continue
                seq_i = candidate_i.representative.sequence

                for j in range(i + 1, len(candidates)):
                    candidate_j = candidates[j]
                    if candidate_j.representative is None:
                        continue
                    seq_j = candidate_j.representative.sequence

                    difference = sum(ci != cj for ci, cj in zip(seq_i, seq_j))
                    if difference <= n_edits:
                        candidate_i.absorb(candidate_j)

        # Remove empty clusters
        clusters = [cluster for cluster in clusters if cluster.representative is not None]

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

    start_time = time.time()
    sequences, total, skipped = read_sequences(args.input)
    read_time = time.time()
    total_unique = sum(len(v) for v in sequences.values())
    print(
        f"Read {total_unique:,} unique sequences ({total:,} total, {skipped:,} skipped) "
        f"in {read_time - start_time:.3g} seconds"
    )

    for seqlen in sorted(sequences.keys()):
        seqs = sequences[seqlen]
        total_counts = sum(s.count for s in sequences[seqlen])
        print(f"Length {seqlen}: Found {len(seqs):,} unique sequences ({total_counts:,} total)")
        clusters = cluster_sequences_same_length(seqs, 2 if args.include_next_nearest else 1)
    cluster_time = time.time()
    print(f"Clustered in {cluster_time - read_time:.3g} seconds")

    # print(f"Found {len(clusters):,} clusters")
    # print(f"Writing clustered sequences to: {args.output}")
    # write_clustered_fasta(clusters, args.output)
    # write_time = time.time()
    # print(f"Writing took: {write_time - cluster_time:.3g} seconds")

    print("Done!")

    print(f"Total time taken: {time.time() - start_time:.3g} seconds")
