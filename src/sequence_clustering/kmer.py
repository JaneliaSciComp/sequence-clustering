"""K-mer utilities shared across stages."""

from __future__ import annotations

import numpy as np
from numba import njit

# DNA alphabet to numeric encoding (two bits per nucleotide).
KMER_ALPHABET = {"A": 0, "C": 1, "G": 2, "T": 3}
BITS_PER_NUCLEOTIDE = 2

# Lookup table for fast char -> int conversion inside Numba-compiled functions.
_ALPHABET_LOOKUP = np.full(128, -1, dtype=np.int8)
for _base, _value in KMER_ALPHABET.items():
    _ALPHABET_LOOKUP[ord(_base)] = _value


def _normalize_sequence(sequence: str) -> str:
    """Return uppercase DNA sequence, validating input."""
    normalized = sequence.strip().upper()
    if not normalized:
        raise ValueError("sequence is empty")
    for base in normalized:
        if base not in KMER_ALPHABET:
            raise ValueError(f"sequence contains invalid nucleotide: {sequence}")
    return normalized


@njit(cache=True)
def _encode_sequence(seq: str) -> np.ndarray:
    encoded = np.empty(len(seq), dtype=np.uint8)
    for i in range(len(seq)):
        val = _ALPHABET_LOOKUP[ord(seq[i])]
        if val < 0:
            raise ValueError("sequence contains invalid nucleotide")
        encoded[i] = np.uint8(val)
    return encoded


@njit(cache=True)
def _sliding_kmer_indices(encoded_seq: np.ndarray, k: int) -> np.ndarray:
    num_kmers = len(encoded_seq) - k + 1
    indices = np.empty(num_kmers, dtype=np.uint32)
    mask_bits = BITS_PER_NUCLEOTIDE * k
    if mask_bits >= 32:
        mask = np.uint32(0xFFFFFFFF)
    else:
        mask = np.uint32((1 << mask_bits) - 1)

    current = np.uint32(0)
    for i in range(k):
        current = (current << BITS_PER_NUCLEOTIDE) | np.uint32(encoded_seq[i])
    indices[0] = current

    for pos in range(1, num_kmers):
        incoming = np.uint32(encoded_seq[pos + k - 1])
        current = ((current << BITS_PER_NUCLEOTIDE) & mask) | incoming
        indices[pos] = current
    return indices


@njit(cache=True)
def _compute_kmer_counts_numba(seq: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    if len(seq) < k:
        return np.empty(0, dtype=np.uint32), np.empty(0, dtype=np.uint16)

    encoded = _encode_sequence(seq)
    indices = _sliding_kmer_indices(encoded, k)

    sorted_indices = np.sort(indices)
    unique_count = 1
    for i in range(1, sorted_indices.size):
        if sorted_indices[i] != sorted_indices[i - 1]:
            unique_count += 1

    unique_indices = np.empty(unique_count, dtype=np.uint32)
    counts = np.empty(unique_count, dtype=np.uint16)

    current_idx = sorted_indices[0]
    current_count = 1
    out_pos = 0

    for i in range(1, sorted_indices.size):
        idx = sorted_indices[i]
        if idx == current_idx:
            current_count += 1
        else:
            unique_indices[out_pos] = current_idx
            counts[out_pos] = current_count
            out_pos += 1
            current_idx = idx
            current_count = 1

    unique_indices[out_pos] = current_idx
    counts[out_pos] = current_count
    return unique_indices, counts


def compute_kmer_counts(sequence: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted sparse k-mer counts for a sequence."""
    if k <= 0:
        raise ValueError("k must be positive")
    seq = _normalize_sequence(sequence)
    return _compute_kmer_counts_numba(seq, k)


@njit(cache=True)
def _reverse_complement_index(index: int, k: int) -> int:
    result = 0
    for _ in range(k):
        digit = index & 0b11
        comp = 0b11 - digit
        result = (result << BITS_PER_NUCLEOTIDE) | comp
        index >>= BITS_PER_NUCLEOTIDE
    return result


@njit(cache=True)
def reverse_complement_permutation(k: int) -> np.ndarray:
    """Return the permutation that maps k-mer index to its reverse complement."""
    if k <= 0:
        raise ValueError("k must be positive")
    dimension = 4**k
    perm = np.empty(dimension, dtype=np.uint32)
    for idx in range(dimension):
        perm[idx] = _reverse_complement_index(idx, k)
    return perm
