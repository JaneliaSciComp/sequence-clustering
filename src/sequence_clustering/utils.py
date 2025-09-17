from abc import ABC, abstractmethod


class LevenshteinCheck(ABC):
    """
    Abstract base class for Levenshtein distance checking with various optimizations.
    """

    def __init__(self, max_distance: int):
        """
        Initialize the Levenshtein checker.

        :param max_distance: Maximum allowed Levenshtein distance.
        """
        self.max_distance = max_distance

    @abstractmethod
    def __call__(self, a: str, b: str) -> bool:
        """
        Check if the Levenshtein distance between two strings is <= max_distance.

        :param a: First string.
        :param b: Second string.
        :returns: True if distance <= max_distance, False otherwise.
        """
        raise NotImplementedError("Subclasses must implement __call__")

    @staticmethod
    def create(max_distance: int, same_length_only: bool = False) -> "LevenshteinCheck":
        """
        Factory method to create the optimal Levenshtein checker implementation.

        :param max_distance: Maximum allowed Levenshtein distance.
        :param same_length_only: If True, optimize for same-length strings only.
        :returns: Optimal LevenshteinCheck implementation.
        """
        if same_length_only:
            if max_distance == 0:
                return LevenshteinCheckSameLength0(max_distance)
            elif max_distance == 1:
                return LevenshteinCheckSameLength1(max_distance)
            elif max_distance == 2:
                return LevenshteinCheckSameLength2(max_distance)
            else:
                return LevenshteinCheckSameLengthGeneral(max_distance)
        else:
            if max_distance == 0:
                return LevenshteinCheck0(max_distance)
            elif max_distance == 1:
                return LevenshteinCheck1(max_distance)
            elif max_distance == 2:
                return LevenshteinCheck2(max_distance)
            else:
                return LevenshteinCheckGeneral(max_distance)


class LevenshteinCheck0(LevenshteinCheck):
    """Optimized implementation for max_distance = 0 (exact match only)."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if strings are identical."""
        return a == b


class LevenshteinCheck1(LevenshteinCheck):
    """Optimized implementation for max_distance = 1."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if Levenshtein distance <= 1."""
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
            return self._check_one_insertion(a, b, n, m)
        else:  # n == m + 1
            # a is one character longer - at most one deletion allowed
            return self._check_one_insertion(b, a, m, n)

    def _check_one_insertion(self, shorter: str, longer: str, n: int, m: int) -> bool:
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


class LevenshteinCheck2(LevenshteinCheck):
    """Optimized implementation for max_distance = 2."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if Levenshtein distance <= 2."""
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


class LevenshteinCheckGeneral(LevenshteinCheck):
    """General implementation for any max_distance using banded algorithm."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if Levenshtein distance <= max_distance using banded DP."""
        n, m = len(a), len(b)

        # Quick length check
        if abs(n - m) > self.max_distance:
            return False

        # Handle edge cases
        if n == 0:
            return m <= self.max_distance
        if m == 0:
            return n <= self.max_distance

        # Use banded algorithm
        band = 2 * self.max_distance + 1

        # Initialize previous row
        prev = [float("inf")] * band
        prev[self.max_distance] = 0  # empty vs empty

        # Fill initial costs for empty string vs prefixes of b
        for j in range(1, min(m + 1, self.max_distance + 1)):
            if self.max_distance - j >= 0:
                prev[self.max_distance - j] = j

        for i in range(1, n + 1):
            curr = [float("inf")] * band

            # Calculate valid range for j
            j_lo = max(1, i - self.max_distance)
            j_hi = min(m, i + self.max_distance)

            for j in range(j_lo, j_hi + 1):
                d = (j - i) + self.max_distance  # diagonal index
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
            if all(val > self.max_distance for val in curr if val != float("inf")):
                return False

            prev = curr

        # Check final result
        final_d = (m - n) + self.max_distance
        if 0 <= final_d < band:
            return prev[final_d] <= self.max_distance
        return False


# Same-length specialized implementations
class LevenshteinCheckSameLength0(LevenshteinCheck):
    """Optimized for same-length strings with max_distance = 0."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if same-length strings are identical."""
        return a == b


class LevenshteinCheckSameLength1(LevenshteinCheck):
    """Optimized for same-length strings with max_distance = 1."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if same-length strings differ by at most 1 substitution."""
        if len(a) != len(b):
            return False

        differences = sum(1 for i in range(len(a)) if a[i] != b[i])
        return differences <= 1


class LevenshteinCheckSameLength2(LevenshteinCheck):
    """Optimized for same-length strings with max_distance = 2."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if same-length strings differ by at most 2 substitutions."""
        if len(a) != len(b):
            return False

        differences = sum(1 for i in range(len(a)) if a[i] != b[i])
        return differences <= 2


class LevenshteinCheckSameLengthGeneral(LevenshteinCheck):
    """Optimized for same-length strings with any max_distance."""

    def __call__(self, a: str, b: str) -> bool:
        """Check if same-length strings differ by at most max_distance substitutions."""
        if len(a) != len(b):
            return False

        differences = sum(1 for i in range(len(a)) if a[i] != b[i])
        return differences <= self.max_distance


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
    checker = LevenshteinCheck.create(max_distance, same_length_only=False)
    return checker(a, b)
