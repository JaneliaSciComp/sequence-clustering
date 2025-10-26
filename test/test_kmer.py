import pytest

from sequence_clustering.kmer import (
    compute_kmer_counts,
    reverse_complement_permutation,
)


def _encode_kmer(kmer: str) -> int:
    idx, _ = compute_kmer_counts(kmer, k=len(kmer))
    return idx[0]


def test_compute_kmer_counts_basic() -> None:
    idxs, cnts = compute_kmer_counts("AACGT", k=3)
    assert idxs.tolist() == [_encode_kmer(k) for k in ("AAC", "ACG", "CGT")]
    assert cnts.tolist() == [1, 1, 1]


def test_compute_kmer_counts_handles_duplicates() -> None:
    idxs, cnts = compute_kmer_counts("AAAAA", k=2)
    assert idxs.tolist() == [_encode_kmer("AA")]
    assert cnts.tolist() == [4]


def test_reverse_complement_permutation_round_trip() -> None:
    perm = reverse_complement_permutation(k=3)
    acg_idx = _encode_kmer("ACG")
    rc_idx = perm[acg_idx]
    assert rc_idx == _encode_kmer("CGT")
    # RC of RC returns original index
    assert perm[rc_idx] == acg_idx


def test_compute_kmer_counts_invalid_sequence() -> None:
    with pytest.raises(ValueError):
        compute_kmer_counts("ABCD", k=3)
