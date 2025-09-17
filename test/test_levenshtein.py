import pytest

from sequence_clustering import (
    check_levenshtein_distance,
    check_levenshtein_distance_same_length,
)


def test_identical_strings():
    """Test that identical strings have distance 0."""
    assert check_levenshtein_distance("abc", "abc", max_distance=0) is True
    assert check_levenshtein_distance("", "", max_distance=0) is True
    assert check_levenshtein_distance("hello", "hello", max_distance=0) is True
    assert check_levenshtein_distance("ATCG", "ATCG", max_distance=0) is True


def test_empty_strings():
    """Test edge cases with empty strings."""
    assert check_levenshtein_distance("", "", max_distance=0) is True
    assert check_levenshtein_distance("", "a", max_distance=0) is False
    assert check_levenshtein_distance("", "a", max_distance=1) is True
    assert check_levenshtein_distance("a", "", max_distance=1) is True
    assert check_levenshtein_distance("", "abc", max_distance=2) is False
    assert check_levenshtein_distance("", "abc", max_distance=3) is True


def test_single_character_differences():
    """Test strings differing by one character."""
    # Substitution
    assert check_levenshtein_distance("a", "b", max_distance=0) is False
    assert check_levenshtein_distance("a", "b", max_distance=1) is True
    assert check_levenshtein_distance("abc", "axc", max_distance=1) is True
    assert check_levenshtein_distance("abc", "axc", max_distance=0) is False

    # Insertion
    assert check_levenshtein_distance("abc", "abxc", max_distance=1) is True
    assert check_levenshtein_distance("abc", "abxc", max_distance=0) is False

    # Deletion
    assert check_levenshtein_distance("abxc", "abc", max_distance=1) is True
    assert check_levenshtein_distance("abxc", "abc", max_distance=0) is False


def test_multiple_differences():
    """Test strings with multiple differences."""
    # Two substitutions
    assert check_levenshtein_distance("abc", "xyz", max_distance=2) is False
    assert check_levenshtein_distance("abc", "xyz", max_distance=3) is True

    # Mixed operations
    assert check_levenshtein_distance("cat", "dog", max_distance=3) is True
    assert check_levenshtein_distance("cat", "dog", max_distance=2) is False

    # Longer strings
    assert check_levenshtein_distance("kitten", "sitting", max_distance=3) is True
    assert check_levenshtein_distance("kitten", "sitting", max_distance=2) is False


def test_length_difference_optimization():
    """Test the length difference optimization."""
    # If length difference > max_distance, should return False quickly
    assert check_levenshtein_distance("a", "abcde", max_distance=3) is False
    assert check_levenshtein_distance("a", "abcde", max_distance=4) is True
    assert check_levenshtein_distance("abc", "abcdefgh", max_distance=4) is False
    assert check_levenshtein_distance("abc", "abcdefgh", max_distance=5) is True


def test_dna_sequences():
    """Test with DNA-like sequences (relevant to the project context)."""
    # Single nucleotide differences
    assert check_levenshtein_distance("ATCG", "TTCG", max_distance=1) is True
    assert check_levenshtein_distance("ATCG", "TTCG", max_distance=0) is False

    # Multiple differences
    assert check_levenshtein_distance("ATCGATCG", "TTCGATCC", max_distance=2) is True
    assert check_levenshtein_distance("ATCGATCG", "TTCGATCC", max_distance=1) is False

    # Insertions/deletions in DNA
    assert check_levenshtein_distance("ATCG", "ATCCG", max_distance=1) is True
    assert check_levenshtein_distance("ATCCG", "ATCG", max_distance=1) is True


def test_longer_sequences():
    """Test with longer sequences to verify performance and correctness."""
    seq1 = "ATCGATCGATCGATCGATCG"
    seq2 = "ATCGATCGATCGATCGATCG"
    assert check_levenshtein_distance(seq1, seq2, max_distance=0) is True

    seq2 = "TTCGATCGATCGATCGATCG"  # 1 difference
    assert check_levenshtein_distance(seq1, seq2, max_distance=1) is True
    assert check_levenshtein_distance(seq1, seq2, max_distance=0) is False

    seq2 = "TTCGATCGATCGATCGATCT"  # 2 differences
    assert check_levenshtein_distance(seq1, seq2, max_distance=2) is True
    assert check_levenshtein_distance(seq1, seq2, max_distance=1) is False


def test_edge_cases_max_distance():
    """Test edge cases with max_distance parameter."""
    # max_distance = 0 should only match identical strings
    assert check_levenshtein_distance("abc", "abc", max_distance=0) is True
    assert check_levenshtein_distance("abc", "abd", max_distance=0) is False

    # Large max_distance should handle most cases
    assert check_levenshtein_distance("abc", "xyz", max_distance=10) is True
    assert check_levenshtein_distance("", "abcdefghij", max_distance=10) is True


def test_symmetric_property():
    """Test that the function is symmetric (distance(a,b) == distance(b,a))."""
    test_cases = [
        ("abc", "def", 3),
        ("hello", "world", 4),
        ("ATCG", "TTCC", 2),
        ("", "abc", 3),
        ("short", "verylongstring", 10),
    ]

    for a, b, max_dist in test_cases:
        result_ab = check_levenshtein_distance(a, b, max_dist)
        result_ba = check_levenshtein_distance(b, a, max_dist)
        assert result_ab == result_ba, f"Symmetry failed for '{a}' and '{b}'"


def test_transitive_property_cases():
    """Test some cases that help verify the algorithm correctness."""
    # If a->b has distance 1 and b->c has distance 1,
    # then a->c has distance at most 2
    a, b, c = "abc", "axc", "ayc"
    assert check_levenshtein_distance(a, b, 1) is True
    assert check_levenshtein_distance(b, c, 1) is True
    assert check_levenshtein_distance(a, c, 2) is True


def test_boundary_conditions():
    """Test boundary conditions for the algorithm."""
    # Test strings where the band just covers the solution
    assert check_levenshtein_distance("abcd", "efgh", max_distance=4) is True
    assert check_levenshtein_distance("abcd", "efgh", max_distance=3) is False

    # Test with strings of very different lengths
    short = "ab"
    long_str = "abcdefghijklmn"
    diff = len(long_str) - len(short)
    assert check_levenshtein_distance(short, long_str, max_distance=diff - 1) is False
    assert check_levenshtein_distance(short, long_str, max_distance=diff) is True


def test_known_examples():
    """Test with well-known Levenshtein distance examples."""
    # Classic example: kitten -> sitting
    # k->s, e->i, insert t = 3 operations
    assert check_levenshtein_distance("kitten", "sitting", max_distance=3) is True
    assert check_levenshtein_distance("kitten", "sitting", max_distance=2) is False

    # Saturday -> Sunday (3 operations)
    assert check_levenshtein_distance("saturday", "sunday", max_distance=3) is True
    assert check_levenshtein_distance("saturday", "sunday", max_distance=2) is False


def test_class_interface_correctness():
    """Test that class interface produces same results as function interface."""
    test_cases = [
        ("", "", 0),
        ("abc", "abc", 0),
        ("abc", "axc", 1),
        ("abc", "xyz", 3),
        ("kitten", "sitting", 3),
        ("ATCG", "TTCG", 1),
        ("ATCGATCG", "TTCGATCC", 2),
    ]

    for a, b, max_dist in test_cases:
        # Test function interface
        func_result = check_levenshtein_distance(a, b, max_dist)

        # Test class interface
        class_result = check_levenshtein_distance(a, b, max_dist)

        assert func_result == class_result, (
            f"Mismatch for '{a}', '{b}', max_dist={max_dist}"
        )


def test_same_length_optimization():
    """Test same-length optimization produces correct results."""
    same_length_cases = [
        ("abc", "abc", 0, True),
        ("abc", "axc", 1, True),
        ("abc", "xyz", 2, False),  # 3 substitutions needed
        ("abc", "xyz", 3, True),
        ("ATCG", "TTCG", 1, True),
        ("ATCG", "TTCC", 2, True),
    ]

    for a, b, max_dist, expected in same_length_cases:
        result = check_levenshtein_distance_same_length(a, b, max_dist)
        assert result == expected, (
            f"Same-length check failed for '{a}', '{b}', max_dist={max_dist}"
        )


def test_same_length_rejects_different_lengths():
    """Test that same-length checkers reject different-length strings."""
    # Different lengths should return False
    with pytest.raises(ValueError):
        check_levenshtein_distance_same_length("abc", "ab", 2)
    with pytest.raises(ValueError):
        check_levenshtein_distance_same_length("ab", "abc", 2)
    with pytest.raises(ValueError):
        check_levenshtein_distance_same_length("", "a", 2)
    with pytest.raises(ValueError):
        check_levenshtein_distance_same_length("a", "", 2)


@pytest.mark.parametrize("max_distance", [0, 1, 2, 5])
def test_all_implementations_consistency(max_distance):
    """Test that all implementations produce consistent results."""
    test_pairs = [
        ("", ""),
        ("a", "a"),
        ("a", "b"),
        ("abc", "axc"),
        ("abc", "xyz"),
        ("hello", "world"),
    ]

    for a, b in test_pairs:
        # Function interface
        func_result = check_levenshtein_distance(a, b, max_distance)

        # Regular class interface
        class_result = check_levenshtein_distance(a, b, max_distance)

        assert func_result == class_result, (
            f"Inconsistency for '{a}', '{b}', max_dist={max_distance}"
        )

        # Same-length interface (only test if lengths match)
        if len(a) == len(b):
            same_len_result = check_levenshtein_distance_same_length(a, b, max_distance)
            assert func_result == same_len_result, (
                f"Same-length inconsistency for '{a}', '{b}', max_dist={max_distance}"
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
