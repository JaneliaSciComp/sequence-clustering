"""Utility script to load example sequences and compute cosine similarities."""

from pathlib import Path
import time

import numpy as np
from scipy.sparse import csr_matrix

from sequence_clustering.io import read_sequences_table
from sequence_clustering.kmer import compute_kmer_counts
from sequence_clustering.dsu import DisjointSetUnion

K = 5
DATA_PATH = Path("data/example_1e5.csv")
THRESHOLD = 0.4


def load_sequences_and_counts(path: Path) -> tuple[list[str], np.ndarray]:
    """Return deduplicated sequences and their read counts."""
    records = read_sequences_table(path)
    sequences = [record.sequence for record in records]
    counts = np.fromiter((record.count for record in records), dtype=np.int64, count=len(records))
    return sequences, counts


def sequences_to_kmer_matrix(sequences: list[str], k: int) -> csr_matrix:
    """Convert sequences into a sparse matrix of k-mer counts."""
    n_rows = len(sequences)
    dim = 4**k
    data: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    indptr = np.zeros(n_rows + 1, dtype=np.int64)

    nnz = 0
    for row, seq in enumerate(sequences, start=1):
        idxs, cnts = compute_kmer_counts(seq, k)
        data.append(cnts.astype(np.float64, copy=False))
        indices.append(idxs.astype(np.int64, copy=False))
        nnz += len(idxs)
        indptr[row] = nnz

    if nnz == 0:
        return csr_matrix((n_rows, dim), dtype=np.float64)

    values = np.concatenate(data)
    col_indices = np.concatenate(indices)
    return csr_matrix((values, col_indices, indptr), shape=(n_rows, dim), dtype=np.float64)


def cosine_similarity(matrix: csr_matrix) -> csr_matrix:
    """Compute cosine similarity between all rows of a CSR matrix."""
    matrix = matrix.todense()
    norms = np.linalg.norm(matrix, axis=1)
    norms[norms == 0.0] = 1.0
    normalized = matrix / norms[:, None]
    return normalized @ normalized.T


def main() -> None:
    start_time = time.time()
    sequences, counts = load_sequences_and_counts(DATA_PATH)
    load_time = time.time()
    print(f"Loaded {len(sequences)} unique sequences in {load_time - start_time:.3f} seconds")

    kmer_matrix = sequences_to_kmer_matrix(sequences, K)
    kmer_time = time.time()
    nnz = kmer_matrix.nnz / np.prod(kmer_matrix.shape) * 100
    print(f"Constructed k-mer matrix ({kmer_matrix.shape}, nnz={nnz:.3f}%) in {kmer_time - load_time:.3f} seconds")

    cosine_matrix = cosine_similarity(kmer_matrix)#.toarray()
    n_links = (np.sum(cosine_matrix > THRESHOLD) - len(sequences)) // 2
    cosine_time = time.time()
    print(f"Computed cosine similarity matrix in {cosine_time - kmer_time:.3f} seconds")
    print(f"Found {n_links} pairs with similarity > {THRESHOLD}")

    # Use DSU to find connected components
    edges = np.argwhere(np.triu(cosine_matrix > THRESHOLD, k=1))
    dsu = DisjointSetUnion(len(sequences))
    for u, v in edges:
        dsu.union(u, v)
    components = dsu.get_components()
    dsu_time = time.time()
    print(f"Found {len(components)} connected components in {dsu_time - cosine_time:.3f} seconds")



if __name__ == "__main__":
    main()
