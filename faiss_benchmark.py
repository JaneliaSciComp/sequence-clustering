"""Utility script to load example sequences and compute cosine similarities."""

from pathlib import Path
import time

import faiss  # type: ignore
import pandas as pd
import numpy as np

from sequence_clustering.kmer import compute_kmer_counts
from sequence_clustering.dsu import DisjointSetUnion

K = 5
DATA_PATH = Path("data/example_1e5.csv")
THRESHOLD = 0.4


def load_sequences_and_counts(path: Path) -> tuple[list[str], np.ndarray]:
    """Return deduplicated sequences and their read counts."""
    records = pd.read_csv(path, comment="#", sep="\t")
    return records["sequence"].tolist(), records["count"].to_numpy()


def sequences_to_kmer_matrix(sequences: list[str], k: int) -> np.ndarray:
    """Convert sequences into a dense matrix of k-mer counts."""
    n_rows = len(sequences)
    dim = 4**k

    # Build the k-mer count matrix
    matrix = np.zeros((n_rows, dim), dtype=np.float32)
    for row, seq in enumerate(sequences):
        idxs, cnts = compute_kmer_counts(seq, k)
        matrix[row, idxs] = cnts.astype(np.float32, copy=False)

    # Normalize rows to unit length
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0.0] = 1.0
    matrix /= norms

    return matrix


def cosine_similarity_brute_force(matrix: np.ndarray, threshold: float) -> np.ndarray:
    """Return all index pairs whose cosine similarity exceeds the threshold."""
    start_time = time.time()
    similarities = matrix @ matrix.T
    sim_time = time.time()
    print(f"Computed cosine similarity matrix in {sim_time - start_time:.3f} seconds")

    row_ids, col_ids = np.nonzero(similarities > threshold)
    mask = row_ids < col_ids
    rows = row_ids[mask]
    cols = col_ids[mask]
    edges = np.column_stack((rows.astype(np.int64), cols.astype(np.int64)))
    edge_time = time.time()
    print(f"Extracted edges in {edge_time - sim_time:.3f} seconds")
    return edges


def cosine_similarity_faiss(matrix: np.ndarray, threshold: float) -> np.ndarray:
    """Return all index pairs whose cosine similarity exceeds the threshold."""
    start_time = time.time()
    n_neighbors = 64

    # Flat index
    index = faiss.IndexFlatIP(matrix.shape[1])

    index.train(matrix)
    index.add(matrix)
    add_time = time.time()
    print(f"Added vectors to FAISS index in {add_time - start_time:.3f} seconds")

    # Compared to knn search, ranged search is exactly what we need here
    similarities, indices = index.search(matrix, k=n_neighbors)
    search_time = time.time()
    print(f"Searched FAISS index in {search_time - add_time:.3f} seconds")

    row_ids = np.arange(matrix.shape[0])[:, None]
    mask = (indices > row_ids) & (similarities > threshold)
    rows, neighbor_positions = np.nonzero(mask)
    cols = indices[rows, neighbor_positions]
    if rows.size == 0:
        return np.empty((0, 2), dtype=np.int64)
    edges = np.column_stack((rows.astype(np.int64), cols.astype(np.int64)))
    edge_time = time.time()
    print(f"Extracted edges in {edge_time - search_time:.3f} seconds")
    return edges



def main() -> None:
    start_time = time.time()
    sequences, _ = load_sequences_and_counts(DATA_PATH)
    load_time = time.time()
    print(f"Loaded {len(sequences)} unique sequences in {load_time - start_time:.3f} seconds")

    kmer_matrix = sequences_to_kmer_matrix(sequences, K)
    kmer_time = time.time()
    density = np.count_nonzero(kmer_matrix) / np.prod(kmer_matrix.shape) * 100
    print(
        f"Constructed k-mer matrix ({kmer_matrix.shape}, nnz={density:.3f}%) "
        f"in {kmer_time - load_time:.3f} seconds"
    )

    edges = cosine_similarity_faiss(kmer_matrix, THRESHOLD)
    cosine_time = time.time()
    print(f"Computed cosine similarity (FAISS flat index) in {cosine_time - kmer_time:.3f} seconds")

    dsu = DisjointSetUnion(len(sequences))
    for u, v in edges:
        dsu.union(int(u), int(v))

    print(f"Found {len(edges)} pairs with similarity > {THRESHOLD}")
    components = dsu.get_components()
    dsu_time = time.time()
    print(f"Found {len(components)} connected components in {dsu_time - cosine_time:.3f} seconds")



if __name__ == "__main__":
    main()
