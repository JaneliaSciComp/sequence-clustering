from collections import defaultdict
from dataclasses import dataclass

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


def is_valid_sequence(seq: str) -> bool:
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    :param seq: The sequence to check.
    :returns: True if sequence contains only A, G, C, T.
    """
    return all(c in 'ATCG' for c in seq)


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

def fill_buckets(sequences: list[UniqueSequence], start: int, end: int) -> dict[str, list[int]]:
    """
    Fill buckets by hashing sequences based on a substring defined by start and end indices.

    :param sequences: List of sequences.
    :param start: Start index of the substring.
    :param end: End index of the substring.
    :returns: A dictionary of buckets, where each bucket is a list of sequence indices.
    """
    seed_to_bucket = defaultdict(list)
    for idx, seq in enumerate(sequences):
        seed = seq.sequence[start:end]
        seed_to_bucket[seed].append(idx)
    return seed_to_bucket


def compare_buckets(
    bucket_a: list[int],
    bucket_b: list[int],
    sequences_a: list[UniqueSequence],
    sequences_b: list[UniqueSequence],
    n_edits: int,
    edges: list[tuple[int, int]],
):
    """
    Compare sequences in two buckets and record edges for sequences within n_edits.

    :param bucket_a: List of indices for the first bucket.
    :param bucket_b: List of indices for the second bucket.
    :param sequences_a: List of sequences corresponding to indices in bucket_a.
    :param sequences_b: List of sequences corresponding to indices in bucket_b.
    :param n_edits: Maximum number of edits allowed to consider sequences as connected.
    :param edges: List to record edges between connected sequences.
    """
    while bucket_a:
        i = bucket_a.pop()
        seq_i = sequences_a[i].sequence

        for j in bucket_b:
            seq_j = sequences_b[j].sequence

            if levenshtein(seq_i, seq_j, n_edits) <= n_edits:
                edges.append((i, j))

