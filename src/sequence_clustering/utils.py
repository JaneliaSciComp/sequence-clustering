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
