# Sequence Clustering Scripts

To run, use
```
python cluster_sequences.py -i merged_lenfilt.fq -o merged_lenfilt_cluster.fasta [--include_next_nearest]
```

## Questions

- Is this a well-defined problem? If three sequences A, B, C are such that A and B differ in 1 position, and A and C differ in 1 position, but B and C differ in 2 positions, the number of clusters depends on the order in which the sequences are processed. If they are processed in the order A, B, C, they will be grouped into one cluster. If they are processed in the order B, C, A, they will be grouped into two clusters (B and C into separate clusters, A into the cluster of whichever sequence is generated first by the perturbation function).
- The library seems to have sequences of different lengths. Do we only cluster based on substitutions, or do we also consider deletes / insertions?
- What will the ultimate input look like? One file with all the sequences, or will the sequences be split into multiple files? What is the expected number of sequences per file?


## Experiments

### 2025-09-15

Ran some preliminary scaling tests that suggest that the FAISS-based approach doesn't really scale well (the non-constant look up time grows too fast).

| Num sequences | `cluster_sequences.py` [s] | `cluster_sequences2.py` [s] | `cluster_sequences_faiss.py` [s] |
|---------------|----------------------------|-----------------------------|----------------------------------|
| 10^4          | 4.07                       | 1.52                        | 0.057                            |
| 10^6          | 155                        | 53                          | 47                               |