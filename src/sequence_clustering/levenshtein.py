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
    inf = max_distance + 1  # Use smaller sentinel value
    band = 2 * max_distance + 1

    # Cost at diagonal d in [-max_distance..max_distance]
    prev = [inf] * band
    curr = [inf] * band
    prev[max_distance] = 0

    for i in range(1, n + 1):
        # Reset curr array
        for k in range(band):
            curr[k] = inf

        j_lo = max(1, i - max_distance)
        j_hi = min(m, i + max_distance)
        min_curr = inf

        for j in range(j_lo, j_hi + 1):
            d = (j - i) + max_distance  # diagonal index in [0..2 * max_distance]

            # Calculate costs directly without redundant conditionals
            sub_cost = prev[d] + (a[i - 1] != b[j - 1])
            ins_cost = curr[d - 1] + 1 if d > 0 else inf
            del_cost = prev[d + 1] + 1 if d + 1 < band else inf

            # Find minimum more efficiently
            if sub_cost <= ins_cost:
                curr[d] = sub_cost if sub_cost <= del_cost else del_cost
            else:
                curr[d] = ins_cost if ins_cost <= del_cost else del_cost

            # Track minimum for early exit
            if curr[d] < min_curr:
                min_curr = curr[d]

        # Early rejection if all values exceed max_distance
        if min_curr > max_distance:
            return False

        # Swap arrays instead of reassigning
        prev, curr = curr, prev

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

    for c1, c2 in zip(a, b):
        if c1 != c2:
            max_distance -= 1
            if max_distance < 0:
                return False
    return True
