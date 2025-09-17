def check_levenshtein_distance(a: str, b: str, max_distance: int) -> bool:
    """
    Check if Levenshtein distance <= max_distance using banded DP.

    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed Levenshtein distance.
    """
    # Quick length bound
    n, m = len(a), len(b)

    if abs(n - m) > max_distance:
        return False
    if n == 0 or m == 0:
        return True

    # Dynamic programming over a band of width 2 * max_distance + 1
    inf = 10**9
    band = 2 * max_distance + 1

    # prev[d] = cost at diagonal d in [-max_distance..max_distance]
    prev = [inf] * band
    prev[max_distance] = 0

    for i in range(1, n + 1):
        curr = [inf] * band
        j_lo = max(1, i - max_distance)
        j_hi = min(m, i + max_distance)

        # preload insertion (on 'a'): moving right in b within this row
        # We'll fill curr from j_lo..j_hi
        for j in range(j_lo, j_hi + 1):
            d = (j - i) + max_distance  # diagonal index in [0..2 * max_distance]
            sub = prev[d] + (a[i - 1] != b[j - 1]) if 0 <= d < band else inf
            ins = curr[d - 1] + 1 if d - 1 >= 0 else inf           # insert into a (advance j)
            dele = prev[d + 1] + 1 if d + 1 < band else inf        # delete from a (advance i)
            curr[d] = sub if sub <= ins and sub <= dele else (ins if ins <= dele else dele)

        if min(curr) > max_distance:   # whole band exceeded max_distance -> early reject
            return False

        prev = curr

    return prev[(m - n) + max_distance] <= max_distance


def check_levenshtein_distance_same_length(a: str, b: str, max_distance: int) -> bool:
    """
    Check if same-length strings differ by at most max_distance
    substitutions.

    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed substitutions.
    :raises ValueError: If strings are of different lengths.
    """
    if len(a) != len(b):
        raise ValueError("Strings must be of the same length")

    differences = sum(1 for i in range(len(a)) if a[i] != b[i])
    return differences <= max_distance
