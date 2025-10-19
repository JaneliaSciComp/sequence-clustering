import csv
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Sequence

import numpy as np
import zarr
from dask.distributed import Client, LocalCluster, as_completed

from .dsu import DisjointSetUnion
from .types import UniqueSequence
from .io import (
    read_sequences_table,
    write_sequences_table,
    FastQReader,
)
from .utils import (
    compare_buckets,
    fill_buckets,
    generate_partitions,
)

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


def read_sequences_table_with_columns(
    path: Path,
    sequence_column: str,
    count_column: str,
) -> list[UniqueSequence]:
    """Load unique sequences from a delimited file with configurable columns."""
    sequences: list[UniqueSequence] = []

    with path.open("r", encoding="utf8") as handle:
        # Detect dialect (in particular, the delimiter)
        try:
            header_line = handle.readline()
            dialect = csv.Sniffer().sniff(header_line)
        except csv.Error as exc:
            raise ValueError(f"Unable to detect delimiter in {path}") from exc

        # Set up CSV reader
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=dialect.delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")

        # Check if sequence and count columns exist
        if (
            sequence_column not in reader.fieldnames
            or count_column not in reader.fieldnames
        ):
            raise ValueError(
                f"Missing required columns {sequence_column}, {count_column} in {path}; "
                f"available: {reader.fieldnames}"
            )

        # Read all (unique) sequences
        for row in reader:
            sequence = row[sequence_column].strip()
            count_str = row[count_column].strip()
            try:
                count = int(count_str)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid count value {count_str!r} in {path}"
                ) from exc
            sequences.append(UniqueSequence(sequence=sequence, count=count))

    return sequences


def split_by_length(
    sequences: Sequence[UniqueSequence],
    output_store: Path,
    chunk_size: int,
) -> None:
    """Write per-length tables into a Zarr store with a group per length."""
    output_store.parent.mkdir(parents=True, exist_ok=True)

    grouped: dict[int, list[UniqueSequence]] = defaultdict(list)
    for record in sequences:
        grouped[len(record.sequence)].append(record)

    sorted_grouped = dict(sorted(grouped.items()))
    root = zarr.open_group(output_store.absolute(), mode="w")

    total_sequences = len(sequences)
    total_reads = sum(record.count for record in sequences)
    root.attrs["total_sequences"] = total_sequences
    root.attrs["total_reads"] = total_reads

    for length, records in sorted_grouped.items():
        length_reads = sum(r.count for r in records)
        print(f"Length {length}: {len(records):,} sequences, {length_reads:,} reads")

        group = root.create_group(f"length_{length}", overwrite=True)
        group.attrs["unique_sequences"] = len(records)
        group.attrs["total_reads"] = length_reads

        dtype = f"<U{length}"
        sequences_arr = np.array([r.sequence for r in records], dtype=dtype)
        counts_arr = np.array([r.count for r in records], dtype=np.int64)

        chunk = min(chunk_size, len(records))
        group.create_dataset(
            "sequence",
            data=sequences_arr,
            chunks=(chunk,),
        )
        group.create_dataset(
            "count",
            data=counts_arr,
            chunks=(chunk,),
        )


def connect_sequences_same_length(
    sequences: Sequence[UniqueSequence], n_edits: int, offset: int
) -> list[tuple[int, int]]:
    """Find all pairs of sequences within an edit distance for same-length sequences."""
    if not sequences:
        return []
    partitions = generate_partitions(len(sequences[0].sequence), n_edits + 1)
    edges: list[tuple[int, int]] = []
    for start, end in partitions:
        seed_to_bucket = fill_buckets(sequences, start, end)
        for bucket in seed_to_bucket.values():
            if len(bucket) < 2:
                continue
            bucket_a = list(bucket)
            bucket_b = list(bucket)
            compare_buckets(bucket_a, bucket_b, sequences, sequences, n_edits, edges)

    # Adjust indices by offset
    edges = [(i + offset, j + offset) for i, j in edges]

    return edges


def connect_sequences_different_length(
    sequences_a: Sequence[UniqueSequence],
    sequences_b: Sequence[UniqueSequence],
    n_edits: int,
    offset_a: int,
    offset_b: int,
) -> list[tuple[int, int]]:
    """Find all pairs of sequences within an edit distance for different-length sequences."""
    if not sequences_a or not sequences_b:
        return []
    len_a = len(sequences_a[0].sequence)
    len_b = len(sequences_b[0].sequence)

    if abs(len_a - len_b) > n_edits:
        return []

    partitions = generate_partitions(len_a, n_edits + 1)
    edges: list[tuple[int, int]] = []
    max_shift = min(n_edits, len_b - len_a)

    for start, end in partitions:
        seed_to_bucket_a = fill_buckets(sequences_a, start, end)
        seed_to_bucket_b: dict[str, list[int]] = {}
        for shift in range(max_shift + 1):
            shifted = fill_buckets(sequences_b, start + shift, end + shift)
            for key, value in shifted.items():
                seed_to_bucket_b.setdefault(key, []).extend(value)

        for seed, bucket_a in seed_to_bucket_a.items():
            bucket_b = seed_to_bucket_b.get(seed)
            if not bucket_b:
                continue
            compare_buckets(
                list(bucket_a), list(bucket_b), sequences_a, sequences_b, n_edits, edges
            )

    # Adjust indices by offsets
    edges = [(i + offset_a, j + offset_b) for i, j in edges]

    return edges


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


def run_split(args) -> None:
    """Split unique sequences into per-length tables."""
    start = time.time()
    input_path = Path(args.input)
    output_store = Path(args.output)
    sequences = read_sequences_table_with_columns(
        input_path,
        sequence_column=args.sequence_column,
        count_column=args.count_column,
    )
    split_by_length(sequences, output_store, args.chunk_size)
    print(f"Split {len(sequences):,} sequences into Zarr groups under {output_store}")
    print(f"Time elapsed: {time.time() - start:.2g} seconds")


def generate_length_pairs(
    lengths: Sequence[int], length_to_count: dict[int, int], max_distance: int
) -> list[tuple[int, int]]:
    """Return all length pairs (a <= b) within the given distance."""
    pairs: list[tuple[int, int]] = []
    lengths = list(sorted(lengths))
    for i, a in enumerate(lengths):
        for b in lengths[i:]:
            if abs(a - b) <= max_distance:
                pairs.append((a, b))

    # Sort them by the expected number of comparisons (product of counts)
    pairs.sort(
        key=lambda ab: length_to_count[ab[0]] * length_to_count[ab[1]],
        reverse=True
    )
    return pairs


def load_length_groups(
    length_store: Path,
    sequence_to_index: dict[str, int],
) -> tuple[dict[int, list[UniqueSequence]], dict[int, list[int]], dict[int, int]]:
    """Load per-length sequences and map them back to global indices."""
    store = zarr.open_group(length_store.absolute(), mode="r")
    pattern = re.compile(r"length_(\d+)")

    sequences_by_length: dict[int, list[UniqueSequence]] = {}
    indices_by_length: dict[int, list[int]] = {}
    counts_by_length: dict[int, int] = {}

    for name, group in store.groups():
        match = pattern.fullmatch(name)
        if not match:
            continue
        length = int(match.group(1))
        seq_array = group["sequence"][:]
        count_array = group["count"][:]
        seq_list = seq_array.tolist()
        count_list = count_array.tolist()

        sequences_list = [
            UniqueSequence(sequence=str(seq), count=int(cnt))
            for seq, cnt in zip(seq_list, count_list, strict=True)
        ]
        if not sequences_list:
            continue

        indices: list[int] = []
        for seq in seq_list:
            idx = sequence_to_index.get(seq)
            if idx is None:
                raise ValueError(
                    f"Sequence {seq!r} in group {name!r} not found in unique table."
                )
            indices.append(idx)

        sequences_by_length[length] = sequences_list
        indices_by_length[length] = indices
        counts_by_length[length] = len(sequences_list)

    return sequences_by_length, indices_by_length, counts_by_length


def compute_edges_for_pair(
    sequences_a: Sequence[UniqueSequence],
    sequences_b: Sequence[UniqueSequence],
    n_edits: int,
    same_length: bool,
) -> list[tuple[int, int]]:
    """Return local index pairs within edit distance between two length buckets."""
    if same_length:
        edges = connect_sequences_same_length(sequences_a, n_edits, 0)
    else:
        edges = connect_sequences_different_length(sequences_a, sequences_b, n_edits, 0, 0)
    return list(set(edges))


def run_cluster(args) -> None:
    """Build clusters by computing edges with Dask and unioning them locally."""
    start = time.time()
    unique_path = Path(args.unique)
    length_store = Path(args.length_store)
    output_path = Path(args.output)
    n_edits = args.distance

    sequences = read_sequences_table(unique_path)
    if not sequences:
        raise ValueError(f"No sequences found in {unique_path}")

    sequence_to_index = {
        record.sequence: idx for idx, record in enumerate(sequences)
    }

    (
        sequences_by_length,
        indices_by_length,
        counts_by_length,
    ) = load_length_groups(length_store, sequence_to_index)

    if not sequences_by_length:
        raise ValueError(
            f"No per-length groups found in {length_store}. "
            "Run the split command first."
        )

    pairs = generate_length_pairs(
        sequences_by_length.keys(),
        counts_by_length,
        n_edits,
    )

    dsu = DisjointSetUnion(len(sequences))
    n_edges = 0

    client: Client | None = None
    cluster: LocalCluster | None = None
    futures = {}

    try:
        if pairs:
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

            scattered_sequences = {
                length: client.scatter(seqs, broadcast=True)
                for length, seqs in sequences_by_length.items()
            }

            for length_a, length_b in pairs:
                future = client.submit(
                    compute_edges_for_pair,
                    scattered_sequences[length_a],
                    scattered_sequences[length_b],
                    n_edits,
                    length_a == length_b,
                    pure=False,
                )
                futures[future] = (length_a, length_b)

            for future in as_completed(futures):
                length_a, length_b = futures[future]
                same_length = length_a == length_b
                try:
                    local_edges = future.result()
                except Exception as exc:  # noqa: BLE001
                    raise RuntimeError(
                        f"Failed to compute edges for lengths ({length_a}, {length_b})"
                    ) from exc

                indices_a = indices_by_length[length_a]
                indices_b = indices_by_length[length_b]

                for local_a, local_b in local_edges:
                    global_a = indices_a[local_a]
                    global_b = indices_a[local_b] if same_length else indices_b[local_b]
                    dsu.union(global_a, global_b)

                n_edges += len(local_edges)
                print(
                    f"Length pair ({length_a}, {length_b}) produced {len(local_edges):,} edges."
                )

        else:
            print("No length pairs within the requested distance. Skipping edge computation.")

    finally:
        for future in list(futures):
            future.release()
        futures.clear()
        if client is not None:
            client.close()
        if cluster is not None:
            cluster.close()

    components = dsu.get_components()
    clusters: list[tuple[str, int, int]] = []
    for component in components:
        total_count = sum(sequences[idx].count for idx in component)
        representative_idx = max(component, key=lambda idx: sequences[idx].count)
        representative = sequences[representative_idx].sequence
        clusters.append((representative, len(component), total_count))

    clusters.sort(key=lambda item: item[2], reverse=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "cluster_size", "total_count"])
        for representative, cluster_size, total_count in clusters:
            writer.writerow([representative, str(cluster_size), str(total_count)])

    elapsed = time.time() - start
    print(
        f"Processed {len(sequences):,} sequences with {n_edges:,} edges "
        f"into {len(clusters):,} clusters."
    )
    print(f"Wrote cluster representatives to {output_path}")
    print(f"Time elapsed: {elapsed:.2g} seconds")
