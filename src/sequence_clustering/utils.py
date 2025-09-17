def check_levenshtein_distance(a: str, b: str, max_distance: int = 0) -> bool:
    """
    Check if the Levenshtein distance between two strings is less than or equal
    to a given maximum distance.

    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed Levenshtein distance.
    :returns: True if the distance is less than or equal to max_distance, False
        otherwise.
    """
    # This uses a dynamic programming approach over a band of width 2k+1
    # centered at diagonal (j-i)

    # Quick length check
    n, m = len(a), len(b)
    if abs(n - m) > max_distance:
        return False

    inf = 10**9
    band = 2 * max_distance + 1

    # prev[d] = cost at diagonal d in [-k..k], index with d+k
    prev = [inf] * band
    prev[max_distance] = 0  # empty vs empty

    for i in range(1, n + 1):
        curr = [inf] * band
        j_lo = max(1, i - max_distance)
        j_hi = min(m, i + max_distance)
        # preload insertion (on 'a'): moving right in b within this row
        # We'll fill curr from j_lo..j_hi
        for j in range(j_lo, j_hi + 1):
            d = (j - i) + max_distance  # diagonal index in [0..2k]
            sub = prev[d] + (a[i - 1] != b[j - 1]) if 0 <= d < band else inf
            ins = curr[d - 1] + 1 if d - 1 >= 0 else inf           # insert into a (advance j)
            dele = prev[d + 1] + 1 if d + 1 < band else inf        # delete from a (advance i)
            curr[d] = sub if sub <= ins and sub <= dele else (ins if ins <= dele else dele)
        if min(curr) > max_distance:   # whole band exceeded k → early reject
            return False
        prev = curr

    return prev[(m - n) + max_distance] <= max_distance
