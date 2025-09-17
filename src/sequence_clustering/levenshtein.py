from functools import partial
from typing import Callable


def create_levenshtein_check(
    max_distance: int, *, same_length_only: bool = False
) -> Callable[[str, str], bool]:
    """
    Initialize the Levenshtein checker.

    :param max_distance: Maximum allowed Levenshtein distance.
    :param same_length_only: If True, optimize for same-length strings only.
    """
    if same_length_only or max_distance == 0:
        return partial(_levenshtein_same_length, max_distance=max_distance)
    elif max_distance == 1:
        return _levenshtein_check_1
    elif max_distance == 2:
        return _levenshtein_check_2
    else:
        return partial(_levenshtein_check_general, max_distance=max_distance)


def check_levenshtein_distance(
    a: str, b: str, max_distance: int, *, same_length_only: bool = False
) -> bool:
    """
    Check if Levenshtein distance between a and b is <= max_distance.
    Note that this is a convenience function that creates a new checker each
    time. For better performance, create a checker with
    :func:`create_levenshtein_check` and reuse it.

    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed Levenshtein distance.
    :param same_length_only: If True, optimize for same-length strings only.
    """
    check = create_levenshtein_check(max_distance, same_length_only=same_length_only)
    return check(a, b)


def _check_one_insertion(shorter: str, longer: str, n: int, m: int) -> bool:
    """Check if longer string can be made from shorter by one insertion."""
    i = j = 0
    found_diff = False
    while i < n and j < m:
        if shorter[i] != longer[j]:
            if found_diff:
                return False
            found_diff = True
            j += 1  # Skip character in longer string
        else:
            i += 1
            j += 1
    return True


def _levenshtein_check_1(a: str, b: str) -> bool:
    """
    Optimized implementation to check if Levenshtein distance <= 1.

    :param a: First string.
    :param b: Second string.
    """
    n, m = len(a), len(b)

    # Quick length check
    if abs(n - m) > 1:
        return False

    # Handle edge cases
    if n == 0:
        return m <= 1
    if m == 0:
        return n <= 1

    if n == m:
        # Same length - at most one substitution allowed
        differences = sum(1 for i in range(n) if a[i] != b[i])
        return differences <= 1
    elif n == m - 1:
        # b is one character longer - at most one insertion allowed
        return _check_one_insertion(a, b, n, m)
    else:  # n == m + 1
        # a is one character longer - at most one deletion allowed
        return _check_one_insertion(b, a, m, n)


def _levenshtein_check_2(a: str, b: str) -> bool:
    """
    Optimized implementation to check if Levenshtein distance <= 2.

    :param a: First string.
    :param b: Second string.
    """
    n, m = len(a), len(b)

    # Quick length check
    if abs(n - m) > 2:
        return False

    # Handle edge cases
    if n == 0:
        return m <= 2
    if m == 0:
        return n <= 2

    # Use standard DP approach with early termination
    prev = list(range(m + 1))

    for i in range(1, n + 1):
        curr = [i]

        for j in range(1, m + 1):
            if a[i - 1] == b[j - 1]:
                cost = prev[j - 1]
            else:
                cost = 1 + min(
                    prev[j],  # deletion
                    curr[j - 1],  # insertion
                    prev[j - 1],  # substitution
                )
            curr.append(cost)

        # Early termination if all values > 2
        if min(curr) > 2:
            return False

        prev = curr

    return prev[m] <= 2


def _levenshtein_check_general(a: str, b: str, max_distance: int) -> bool:
    """
    Check if Levenshtein distance <= max_distance using banded DP.
    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed Levenshtein distance.
    """
    n, m = len(a), len(b)

    # Quick length check
    if abs(n - m) > max_distance:
        return False

    # Handle edge cases
    if n == 0:
        return m <= max_distance
    if m == 0:
        return n <= max_distance

    # Use banded algorithm
    band = 2 * max_distance + 1

    # Initialize previous row
    prev = [float("inf")] * band
    prev[max_distance] = 0  # empty vs empty

    # Fill initial costs for empty string vs prefixes of b
    for j in range(1, min(m + 1, max_distance + 1)):
        if max_distance - j >= 0:
            prev[max_distance - j] = j

    for i in range(1, n + 1):
        curr = [float("inf")] * band

        # Calculate valid range for j
        j_lo = max(1, i - max_distance)
        j_hi = min(m, i + max_distance)

        for j in range(j_lo, j_hi + 1):
            d = (j - i) + max_distance  # diagonal index
            if d < 0 or d >= band:
                continue

            # Calculate costs
            cost = 0 if a[i - 1] == b[j - 1] else 1

            # Substitution
            sub_cost = float("inf")
            if 0 <= d < band:
                sub_cost = prev[d] + cost

            # Insertion (advance j)
            ins_cost = float("inf")
            if d - 1 >= 0:
                ins_cost = curr[d - 1] + 1

            # Deletion (advance i)
            del_cost = float("inf")
            if d + 1 < band:
                del_cost = prev[d + 1] + 1

            curr[d] = min(sub_cost, ins_cost, del_cost)

        # Early termination if all values > max_distance
        if all(val > max_distance for val in curr if val != float("inf")):
            return False

        prev = curr

    # Check final result
    final_d = (m - n) + max_distance
    if 0 <= final_d < band:
        return prev[final_d] <= max_distance
    return False


# Same-length specialized implementation
def _levenshtein_same_length(a: str, b: str, max_distance: int) -> bool:
    """
    Check if same-length strings differ by at most max_distance
    substitutions.

    :param a: First string.
    :param b: Second string.
    :param max_distance: Maximum allowed substitutions.
    """
    if len(a) != len(b):
        return False

    differences = sum(1 for i in range(len(a)) if a[i] != b[i])
    return differences <= max_distance
