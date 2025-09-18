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


def is_valid_sequence(str seq):
    """
    Check if sequence contains only valid nucleotides (A, G, C, T).

    :param seq: The sequence to check.
    :returns: True if sequence contains only A, G, C, T.
    """
    cdef int i
    cdef int seq_len = len(seq)
    cdef str c
    
    for i in range(seq_len):
        c = seq[i]
        if c != 'A' and c != 'T' and c != 'C' and c != 'G':
            return False
    return True


def generate_partitions(int n, int k):
    """
    Partition n items into k contiguous sets of roughly equal size.
    Cython-optimized with typed variables.

    :param n: Total number of items.
    :param k: Number of subsets.
    :returns: A list of all partitions.
    """
    cdef int base_size = n // k
    cdef int remainder = n % k
    cdef int start = 0
    cdef int size, end
    cdef int i
    
    partitions = []
    for i in range(k):
        size = base_size + (1 if i < remainder else 0)
        end = start + size
        partitions.append((start, end))
        start = end

    return partitions


def fill_buckets(sequences: list[UniqueSequence], int start, int end):
    """
    Fill buckets by hashing sequences based on a substring defined by start and end indices.

    :param sequences: List of sequences.
    :param start: Start index of the substring.
    :param end: End index of the substring.
    :returns: A dictionary of buckets, where each bucket is a list of sequence indices.
    """
    seed_to_bucket = defaultdict(list)
    cdef int idx
    cdef object seq  # Use object instead of specific type
    cdef str seed
    
    for idx, seq in enumerate(sequences):
        seed = seq.sequence[start:end]
        seed_to_bucket[seed].append(idx)
    return seed_to_bucket


def compare_buckets(
    bucket_a: list[int],
    bucket_b: list[int],
    sequences_a: list[UniqueSequence],
    sequences_b: list[UniqueSequence],
    int n_edits,
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
    cdef int i, j
    cdef str seq_i, seq_j
    cdef int distance
    
    while bucket_a:
        i = bucket_a.pop()
        seq_i = sequences_a[i].sequence

        for j in bucket_b:
            seq_j = sequences_b[j].sequence
            
            # Use the levenshtein function with early termination
            distance = levenshtein(seq_i, seq_j, n_edits)
            if distance <= n_edits:
                edges.append((i, j))
