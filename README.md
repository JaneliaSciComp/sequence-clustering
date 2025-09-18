# Sequence Clustering Scripts

A small set of scripts for clustering sequences using the following algorithm:
- Find all unique sequences in the dataset, recording their absolute frequency.
- Find all pairs of unique sequences that are within max_edits (1 or 2) from each other; for this
  - Use a bucketing algorithm that's guaranteed to not have false negatives (multi-index hashing with max_edits + 1 seeds, shifting the seeds for sequences of different lengths)
  - Compute the Levenshtein distance all pairs of sequences within a bucket (the buckets are pretty small on average)
  - Record all pairs where the distance is <= max_edits
- Find the topological components of the graph represented by these pairs -> these are the clusters (two sequences in such a cluster can have a distance more than max_edits, but they are always linked by a chain of at most max_edits edits; sequences of two distinct clusters are separated by more than max_edits edits)
- Use the sequence with the highest frequency as a representative of the cluster and write the representatives to a .fasta file

## Setup

The project is implemented in Python and Cython and uses `uv` for dependency management and execution.

1. Download [`uv`](https://github.com/astral-sh/uv) and install it.

2. Install the package: `uv sync` (optional, since running the scripts via `uv run <script>` will automatically install dependencies).

3. Run the code: `uv run cluster_sequences -i <input.fq> -o <output.fasta> [--include_next_nearest]`


## Experiments

### 2025-09-15

Ran some preliminary scaling tests that suggest that the FAISS-based approach doesn't really scale well (the non-constant look up time grows too fast).

| Num sequences | `cluster_sequences.py` [s] | `cluster_sequences2.py` [s] | `cluster_sequences_faiss.py` [s] |
|---------------|----------------------------|-----------------------------|----------------------------------|
| 10^4          | 4.07                       | 1.52                        | 0.057                            |
| 10^6          | 155                        | 53                          | 47                               |

### 2025-09-17

Implemented a new approach that uses bucketing based on substrings to reduce the number of pairwise comparisons and implemented some parts in Cython for speedup.

| Num sequences | `__main__.py` [s] |
|---------------|-------------------|
| 10^4          | 0.020              |
| 10^6          | 4.13               |
