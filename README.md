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

3. Run the code: `uv run -m sequence_clustering -i <input.fq> -o <output.fasta> [--include_next_nearest]`


## Experiments

### 2025-09-15

Ran some preliminary scaling tests that suggest that an approach using the FAISS library doesn't really bring the expected speedup.
Given that FAISS indexing doesn't scale linearly (although with a very small exponent), this approach likely won't scale to very large datasets.
See an [older commit](https://github.com/JaneliaSciComp/sequence-clustering/commit/e5e4a3fcd3e834a2f8fa338f4b94c7cd2bd64351) to reproduce these results.

| Num sequences | `cluster_sequences.py` [s] | `cluster_sequences2.py` [s] | `cluster_sequences_faiss.py` [s] |
|---------------|----------------------------|-----------------------------|----------------------------------|
| 10^4          | 4.07                       | 1.52                        | 0.057                            |
| 10^6          | 155                        | 53                          | 47                               |

### 2025-09-17

Implemented a new approach that uses bucketing based on substrings to reduce the number of pairwise comparisons and implemented some parts in Cython for speedup.

| Num sequences | `uv run -m cluster_sequences` [s] |
|---------------|-----------------------------------|
| 10^4          | 0.020                             |
| 10^6          | 4.13                              |

### 2025-09-18

Compared the new approach to the previous Python implementation in terms of runtime and clustering results.
The number of clusters is reported in the following table, with the number of clusters that are present in the new implementation but not in the original ones in parentheses.
The new implementation seems to produce a good clustering (maybe even slightly better), but is significantly faster.

Note:
- The original implementations only consider substitutions, not indels.
- `cluster_sequences.py` produces all 2-neighbors of each sequence and compares them to all 2-neighbors of all other sequences, thereby maybe creating links of distance 4 that are not present in the data.
This is probably why it produces slightly different results than the new implementation.
- `cluster_sequences2.py` only considers distances up to 2, but only compares sequences by order of appearance in the input file, thereby potentially missing some chained links.
This is probably why it produces significantly more clusters than the new implementation.

| Num sequences | `cluster_sequences.py` | `cluster_sequences2.py` | `uv run -m sequence_clustering` |
|---------------|------------------------|-------------------------|---------------------------------|
| 10^4          | 2,253 (59)             | 2,350 (0)               | 2,124                           |
| 10^6          | 52,534 (7,459)         | 73,166 (0)              | 46,492                          |
