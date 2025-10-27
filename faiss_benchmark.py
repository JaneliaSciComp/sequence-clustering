"""Utility script to load example sequences and compute cosine similarities."""

from pathlib import Path
import time

import faiss  # type: ignore
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
    """Convert sequences into a sparse CSR matrix of k-mer counts."""
    n_rows = len(sequences)
    dim = 4**k
    data: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    indptr = np.zeros(n_rows + 1, dtype=np.int64)

    nnz = 0
    for row, seq in enumerate(sequences, start=1):
        idxs, cnts = compute_kmer_counts(seq, k)
        data.append(cnts.astype(np.float32, copy=False))
        indices.append(idxs.astype(np.int64, copy=False))
        nnz += len(idxs)
        indptr[row] = nnz

    if nnz == 0:
        return csr_matrix((n_rows, dim), dtype=np.float32)

    values = np.concatenate(data)
    col_indices = np.concatenate(indices)
    return csr_matrix((values, col_indices, indptr), shape=(n_rows, dim), dtype=np.float32)


def cosine_similarity_numpy(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute cosine similarity via NumPy."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    normalized = (matrix / norms).astype(np.float32, copy=False)

    similarities = normalized @ normalized.T
    indices = np.argsort(-similarities, axis=1)
    return similarities, indices

def cosine_similarity_faiss(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute cosine similarity via FAISS IndexFlatIP."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    normalized = (matrix / norms).astype(np.float32, copy=False)

    index = faiss.IndexFlatIP(normalized.shape[1])
    index.add(normalized)
    similarities, indices = index.search(normalized, normalized.shape[0])
    return similarities, indices


def main() -> None:
    start_time = time.time()
    sequences, counts = load_sequences_and_counts(DATA_PATH)
    load_time = time.time()
    print(f"Loaded {len(sequences)} unique sequences in {load_time - start_time:.3f} seconds")

    kmer_matrix = sequences_to_kmer_matrix(sequences, K)
    kmer_time = time.time()
    density = kmer_matrix.nnz / np.prod(kmer_matrix.shape) * 100
    print(f"Constructed k-mer matrix ({kmer_matrix.shape}, nnz={density:.3f}%) in {kmer_time - load_time:.3f} seconds")

    dense_matrix = kmer_matrix.toarray()
    similarities, neighbors = cosine_similarity_faiss(dense_matrix)
    cosine_time = time.time()
    print(f"Computed cosine similarity (FAISS flat index) in {cosine_time - kmer_time:.3f} seconds")

    dsu = DisjointSetUnion(len(sequences))
    edge_count = 0
    for row, (row_scores, row_neighbors) in enumerate(zip(similarities, neighbors), start=0):
        for score, col in zip(row_scores, row_neighbors):
            if col <= row:
                continue
            if score > THRESHOLD:
                dsu.union(row, col)
                edge_count += 1

    print(f"Found {edge_count} pairs with similarity > {THRESHOLD}")
    components = dsu.get_components()
    dsu_time = time.time()
    print(f"Found {len(components)} connected components in {dsu_time - cosine_time:.3f} seconds")



if __name__ == "__main__":
    main()
