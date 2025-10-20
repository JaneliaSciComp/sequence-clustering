import csv
import time
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import zarr
from dask.distributed import Client, LocalCluster, as_completed

from .dsu import DisjointSetUnion
from .types import UniqueSequence
from .io import (
    write_sequences_table,
    read_sequences_table,
    FastQReader,
    ZarrStoreByLength,
)
from .utils import (
    compare_buckets,
    fill_buckets,
    generate_partitions,
)


@dataclass
class TileSpec:
    """Specification for a tile of sequences with specific length."""
    sequence_length: int
    offset: int
    start: int
    end: int


def run_unique(args) -> None:
    """Extract unique sequences from a FASTQ file."""
    start = time.time()
    fastq_path = Path(args.fastq)
    output_path = Path(args.output)
    sequences, total_reads, skipped = collect_unique_sequences(fastq_path)
    write_sequences_table(sequences, output_path)
    print(
        f"Found {len(sequences):,} unique sequences "
        f"({total_reads:,} total reads, {skipped:,} skipped)."
    )
    print(f"Wrote unique sequence table to {output_path}")
    print(f"Time elapsed: {time.time() - start:.2g} seconds")


def collect_unique_sequences(fastq_path: Path) -> tuple[list[UniqueSequence], int, int]:
    """Return unique sequences, total reads, and skipped reads from a FASTQ file."""
    counts: dict[str, int] = defaultdict(int)
    total_reads = 0
    skipped_reads = 0

    with FastQReader(fastq_path) as handle:
        for sequence in handle:
            counts[sequence] += 1

        total_reads = handle.total
        skipped_reads = handle.skipped

    ordered = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    uniques = [
        UniqueSequence(sequence=seq, count=count)
        for (seq, count) in ordered
    ]
    return uniques, total_reads, skipped_reads


def run_cluster(args) -> None:
    """Build clusters by computing edges with Dask and unioning them locally."""
    start = time.time()
    unique_path = Path(args.input)
    length_store = args.length_store or unique_path.parent / "by_length.zarr"
    output_path = Path(args.output)
    n_edits = args.distance

    # Split sequences by length to make accessing them easier
    print(f"Loading unique sequences from '{unique_path}'...")
    split_by_length(
        unique_path,
        length_store,
        chunk_size=args.chunk_size,
        sequence_column=args.sequence_column,
        count_column=args.count_column,
    )
    print(f"Wrote per-length Zarr store to '{length_store}'")

    # Generate all sequence pairs to compare
    length_to_total_counts = load_total_counts(length_store)
    total_count = sum(length_to_total_counts.values())
    pairs = generate_length_pairs(length_to_total_counts, n_edits, args.tile_size)
    if not pairs:
        print("No length pairs within the requested distance.")
        return

    dsu = DisjointSetUnion(total_count)
    n_edges = 0

    try:
        # Start a local Dask cluster
        cluster = LocalCluster(
            n_workers=args.workers or None,
            threads_per_worker=args.threads_per_worker or None,
            dashboard_address=None,
        )
        client = Client(cluster)
        nthreads = client.nthreads()
        print(
            f"Started local Dask cluster with {len(nthreads)} workers "
            f"and thread distribution {sorted(nthreads.values())}."
        )

        # Submit all length pairs as separate tasks (in tiles)
        futures = [
            client.submit(
                compute_edges_for_pair,
                length_store,
                tile_spec_a,
                tile_spec_b,
                n_edits,
            )
            for tile_spec_a, tile_spec_b in pairs
        ]

        # Collect results as they complete and aggregate edges
        for future, edges in as_completed(futures, with_results=True):
            for global_i, global_j in edges:
                dsu.union(global_i, global_j)
            n_edges += len(edges)
            future.release()

    finally:
        for future in futures:
            future.release()
        futures.clear()
        if client is not None:
            client.close()
        if cluster is not None:
            cluster.close()

    # Load all read counts for cluster assembly
    zarr_store = ZarrStoreByLength(length_store)
    counts = np.zeros(total_count, dtype=np.int64)
    sequences_flat: list[str] = []
    start_idx = 0
    for length in sorted(length_to_total_counts):
        local_counts = zarr_store.load_counts(length)
        local_sequences = zarr_store.load_sequences(length)
        if len(local_counts) != len(local_sequences):
            raise ValueError(
                f"Mismatched sequences/counts for length {length} in {length_store}"
            )
        end_idx = start_idx + len(local_counts)
        counts[start_idx:end_idx] = local_counts
        sequences_flat.extend(local_sequences)
        start_idx = end_idx

    if start_idx != total_count:
        raise ValueError(
            f"Expected {total_count} total sequences but reconstructed {start_idx}"
        )

    # Assemble clusters from the union-find structure
    components = dsu.get_components()
    clusters: list[tuple[int, int, int]] = []
    for component in components:
        component_total = sum(counts[idx] for idx in component)
        representative_idx = max(component, key=lambda idx: counts[idx])
        clusters.append((representative_idx, len(component), component_total))

    # Dereference representative sequences (and sort by total count)
    del dsu
    clusters = [
        (sequences_flat[rep_idx], cluster_size, total_count)
        for rep_idx, cluster_size, total_count in clusters
    ]
    clusters.sort(key=lambda item: item[2], reverse=True)

    # Write out cluster representatives
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "cluster_size", "total_count"])
        for representative, cluster_size, total_count in clusters:
            writer.writerow([representative, str(cluster_size), str(total_count)])

    elapsed = time.time() - start
    print(
        f"Processed {len(sequences_flat):,} sequences with {n_edges:,} edges "
        f"into {len(clusters):,} clusters."
    )
    print(f"Wrote cluster representatives to {output_path}")
    print(f"Time elapsed: {elapsed:.2g} seconds")


def split_by_length(
    input_file: Path,
    output_store: Path,
    chunk_size: int,
    sequence_column: str,
    count_column: str,
) -> None:
    """Write per-length tables into a Zarr store with a group per length."""
    # Read all sequences from the input csv file
    sequences = read_sequences_table(
        input_file,
        sequence_column,
        count_column,
    )

    # Write sequences to Zarr store by length
    ZarrStoreByLength.write(sequences, output_store, chunk_size)


def generate_length_pairs(
    lengths_to_total_counts: dict[int, int], max_distance: int, tile_size: int,
) -> list[tuple[TileSpec, TileSpec]]:
    """
    Return all length pairs (a <= b) within the given distance, tiled if
    there are too many sequences.
    """
    pairs: list[tuple[TileSpec, TileSpec]] = []
    lengths = list(sorted(lengths_to_total_counts.keys()))

    # Compute offsets into the global sequence list by length
    length_to_offset: dict[int, int] = {}
    offset = 0
    for length in lengths:
        length_to_offset[length] = offset
        offset += lengths_to_total_counts[length]

    # Generate all length pairs within the given constraints
    for i, a in enumerate(lengths):
        for b in lengths[i:]:
            if abs(a - b) > max_distance:
                continue

            offset_a = length_to_offset[a]
            offset_b = length_to_offset[b]
            total_counts_a = lengths_to_total_counts[a]
            total_counts_b = lengths_to_total_counts[b]

            for i in range(0, total_counts_a, tile_size):
                tile_spec_i = TileSpec(
                    sequence_length=a,
                    offset=offset_a + i,
                    start=i,
                    end=min(i + tile_size, total_counts_a),
                )
                for j in range(0, total_counts_b, tile_size):
                    if a == b and i > j:
                        continue  # Avoid duplicate tiles for same-length pairs

                    tile_spec_j = TileSpec(
                        sequence_length=b,
                        offset=offset_b + j,
                        start=j,
                        end=min(j + tile_size, total_counts_b),
                    )
                    pairs.append((tile_spec_i, tile_spec_j))

    return pairs


def load_total_counts(length_store: Path) -> dict[int, int]:
    """Load per-length sequence numbers."""
    store = zarr.open_group(str(length_store), mode="r")

    length_to_total_counts: dict[int, int] = {}
    for _, group in store.groups():
        length = int(group.attrs["sequence_length"])
        total_counts = int(group.attrs["unique_sequences"])
        length_to_total_counts[length] = total_counts

    return length_to_total_counts


def compute_edges_for_pair(
    length_store: Path,
    tile_spec_a: TileSpec,
    tile_spec_b: TileSpec,
    n_edits: int,
) -> list[tuple[int, int]]:
    """Return local index pairs within edit distance between two length buckets."""
    # Read sequences of the given lengths
    zarr_store = ZarrStoreByLength(length_store)
    sequences_a = zarr_store.load_sequences(
        tile_spec_a.sequence_length,
        slice(tile_spec_a.start, tile_spec_a.end)
    )
    sequences_b = zarr_store.load_sequences(
        tile_spec_b.sequence_length,
        slice(tile_spec_b.start, tile_spec_b.end)
    )
    partitions = generate_partitions(tile_spec_a.sequence_length, n_edits + 1)
    edges: list[tuple[int, int]] = []
    max_shift = min(n_edits, tile_spec_b.sequence_length - tile_spec_a.sequence_length) + 1

    # Generate buckets and compare within each bucket
    total_buckets = 0
    start_time = time.time()
    for start, end in partitions:
        seed_to_bucket_a = fill_buckets(sequences_a, start, end)
        seed_to_bucket_b: dict[str, list[int]] = defaultdict(list)
        for shift in range(max_shift):
            shifted = fill_buckets(sequences_b, start + shift, end + shift)
            for key, value in shifted.items():
                seed_to_bucket_b[key].extend(value)

        for seed, bucket_a in seed_to_bucket_a.items():
            bucket_b = seed_to_bucket_b.get(seed)
            if not bucket_b:
                continue
            total_buckets += len(bucket_a) * len(bucket_b)
            compare_buckets(
                bucket_a, bucket_b,
                sequences_a, sequences_b,
                n_edits, edges
            )

    total_pairwise = len(sequences_a) * len(sequences_b)
    elapsed = time.time() - start_time

    # Remove duplicate edges and apply offset
    edges = [(a + tile_spec_a.offset, b + tile_spec_b.offset) for (a, b) in edges]
    edges = DisjointSetUnion.deduplicate_edges(edges)
    print(
        f"Lengths {tile_spec_a.sequence_length} and {tile_spec_b.sequence_length}: "
        f"Found {len(edges):,} edges (performed {total_buckets:,} of {total_pairwise:,} "
        f"comparisons) in {elapsed:.2f} seconds."
    )
    return edges
