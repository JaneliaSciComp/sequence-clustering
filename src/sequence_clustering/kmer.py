"""K-mer utilities shared across stages."""

from collections import Counter
from typing import Iterable

import numpy as np

KMER_ALPHABET = {"A": 0, "C": 1, "G": 2, "T": 3}
BITS_PER_NUCLEOTIDE = 2


def _normalize_sequence(seq: str) -> str:
    normalized = seq.strip().upper()
    if not normalized:
        raise ValueError("sequence is empty")
    if any(base not in KMER_ALPHABET for base in normalized):
        raise ValueError(f"sequence contains invalid nucleotide: {seq}")
    return normalized


def _sliding_kmer_indices(seq: str, k: int) -> Iterable[int]:
    mask = (1 << (BITS_PER_NUCLEOTIDE * k)) - 1
    current = 0
    for char in seq[:k]:
        current = (current << BITS_PER_NUCLEOTIDE) | KMER_ALPHABET[char]
    yield current
    for char in seq[k:]:
        current = ((current << BITS_PER_NUCLEOTIDE) & mask) | KMER_ALPHABET[char]
        yield current


def compute_kmer_counts(sequence: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted sparse k-mer counts for a sequence."""
    if k <= 0:
        raise ValueError("k must be positive")
    seq = _normalize_sequence(sequence)
    if len(seq) < k:
        return np.array([], dtype=np.uint32), np.array([], dtype=np.uint16)

    counts: Counter[int] = Counter(_sliding_kmer_indices(seq, k))
    if not counts:
        return np.array([], dtype=np.uint32), np.array([], dtype=np.uint16)
    items = sorted(counts.items())
    idxs = np.fromiter((idx for idx, _ in items), dtype=np.uint32, count=len(items))
    cnts = np.fromiter((cnt for _, cnt in items), dtype=np.uint16, count=len(items))
    return idxs, cnts


def reverse_complement_permutation(k: int) -> np.ndarray:
    """Return the permutation that maps k-mer index to its reverse complement."""
    if k <= 0:
        raise ValueError("k must be positive")
    dimension = 4 ** k
    perm = np.empty(dimension, dtype=np.uint32)
    for idx in range(dimension):
        perm[idx] = _reverse_complement_index(idx, k)
    return perm


def _reverse_complement_index(index: int, k: int) -> int:
    result = 0
    for _ in range(k):
        digit = index & 0b11
        comp = 0b11 - digit
        result = (result << BITS_PER_NUCLEOTIDE) | comp
        index >>= BITS_PER_NUCLEOTIDE
    return result
