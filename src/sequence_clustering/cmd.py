import csv
import time
import sys
import logging
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import zarr
from dask.distributed import Client, as_completed, LocalCluster
from dask_jobqueue import LSFCluster

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


# Set up logging
logger = logging.getLogger(__name__)
if not logger.handlers:
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(levelname)s: %(message)s",
            "%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(handler)
    logger.propagate = False


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
    logger.info(
        "Found %d unique sequences (%d total reads, %d skipped).",
        len(sequences),
        total_reads,
        skipped,
    )
    logger.info("Wrote unique sequence table to %s", output_path)
    logger.info("Time elapsed: %.2g seconds", time.time() - start)


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
    length_store = (
        Path(args.length_store)
        if args.length_store
        else unique_path.parent / "by_length.zarr"
    )
    output_path = Path(args.output)
    n_edits = args.distance

    # Split sequences by length to make accessing them easier
    logger.info("Loading unique sequences from %s...", unique_path)
    split_by_length(
        unique_path,
        length_store,
        chunk_size=args.tile_size,
        sequence_column=args.sequence_column,
        count_column=args.count_column,
    )
    logger.info("Wrote per-length Zarr store to %s", length_store)

    # Generate all sequence pairs to compare
    length_to_total_counts = load_total_counts(length_store)
    total_count = sum(length_to_total_counts.values())
    pairs = generate_length_pairs(length_to_total_counts, n_edits, args.tile_size)
    if not pairs:
        logger.info("No length pairs within the requested distance.")
        return

    # Summarize some tiles to only have ~10 * workers tasks
    pair_groups = []
    chunk_size = max(1, len(pairs) // (10 * args.workers))
    for i in range(0, len(pairs), chunk_size):
        pair_groups.append(pairs[i:i + chunk_size])

    dsu = DisjointSetUnion(total_count)
    n_edges = 0

    # Start a local Dask cluster to process all length pairs
    with create_dask_cluster(args) as dask_cluster, Client(dask_cluster) as client:
        # Submit all length pairs as separate tasks (in tiles)
        futures = [
            client.submit(
                compute_edges_for_all_pairs,
                length_store,
                pairs,
                n_edits,
            )
            for pairs in pair_groups
        ]

        # Collect results as they complete and aggregate edges
        for i, (future, edges) in enumerate(as_completed(futures, with_results=True)):
            logger.info("Finished task %d / %d", i + 1, len(futures))
            for global_i, global_j in edges:
                dsu.union(global_i, global_j)
            n_edges += len(edges)
            future.release()

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
    logger.info("Processed %d sequences with %d edges into %d clusters.",
        len(sequences_flat), n_edges, len(clusters)
    )
    logger.info("Wrote cluster representatives to %s", output_path)
    logger.info("Time elapsed: %.2g seconds", elapsed)


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


def create_dask_cluster(args):
    """Create a Dask cluster based on the specified parallelism strategy."""
    if args.parallel == "local":
        cluster = LocalCluster(
            n_workers=args.workers,
            threads_per_worker=args.threads_per_worker,
        )
        logger.info("Started local Dask cluster with %d workers.", args.workers)
        return cluster
    elif args.parallel == "lsf":
        cluster = LSFCluster(
            queue="local",
            project="das",
            n_workers=args.workers,
            cores=args.threads_per_worker,
            log_directory="dask-logs",
            memory=f"{15 * args.threads_per_worker}GB",  # ignored by scheduler, but not by nanny
            walltime="24:00",  # set a reasonable walltime
            job_script_prologue=["export PYTHONUNBUFFERED=1"],  # unbuffer Python stdio
        )
        logger.info("Started LSF Dask cluster with %d workers.", args.workers)
        return cluster
    else:
        raise ValueError(f"Unknown parallelism strategy: {args.parallel}")


def load_total_counts(length_store: Path) -> dict[int, int]:
    """Load per-length sequence numbers."""
    store = zarr.open_group(str(length_store), mode="r")

    length_to_total_counts: dict[int, int] = {}
    for _, group in store.groups():
        length = int(group.attrs["sequence_length"])
        total_counts = int(group.attrs["unique_sequences"])
        length_to_total_counts[length] = total_counts

    return length_to_total_counts


def compute_edges_for_all_pairs(
    length_store: Path,
    pairs: list[tuple[TileSpec, TileSpec]],
    n_edits: int,
) -> list[tuple[int, int]]:
    """Return all edges within edit distance for the given length tile pairs."""
    all_edges: list[tuple[int, int]] = []
    for tile_spec_a, tile_spec_b in pairs:
        edges = compute_edges_for_pair(
            length_store,
            tile_spec_a,
            tile_spec_b,
            n_edits,
        )
        edges = [(a + tile_spec_a.offset, b + tile_spec_b.offset) for (a, b) in edges]
        all_edges.extend(edges)

    # Remove duplicate edges and apply offset
    old_edges_count = len(all_edges)
    all_edges = DisjointSetUnion.deduplicate_edges(all_edges)
    logger.info("Deduplicated %d edges to %d edges.", old_edges_count, len(all_edges))

    return all_edges

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
    logger.info(
        "Lengths %d and %d: Found %d edges (performed %s of %s comparisons) in %.2f seconds.",
        tile_spec_a.sequence_length,
        tile_spec_b.sequence_length,
        len(edges),
        format(total_buckets, ","),
        format(total_pairwise, ","),
        elapsed,
    )

    return edges
