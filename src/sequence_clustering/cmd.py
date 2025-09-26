from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from .dsu import DisjointSetUnion
from .types import UniqueSequence
from .io import read_sequences_table, write_sequences_table, FastQReader
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


def split_by_length(
    sequences: Sequence[UniqueSequence],
    output_dir: Path,
    min_length: int | None,
    max_length: int | None,
) -> None:
    """Write per-length tables with summary headers."""
    # Group sequences by length
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[int, list[UniqueSequence]] = defaultdict(list)
    for record in sequences:
        grouped[len(record.sequence)].append(record)
    grouped = dict(sorted(grouped.items()))

    # Filter lengths outside the specified range
    min_length = min_length or min(grouped.keys(), default=0)
    max_length = max_length or max(grouped.keys(), default=0)

    # Write per-length files
    for length, records in grouped.items():
        length_reads = sum(r.count for r in records)
        if length < min_length or length > max_length:
            print(
                f"Skipping length {length} "
                f"({len(records):,} sequences, {length_reads:,} reads)"
            )
            continue

        length_path = output_dir / f"length_{length}.csv"
        print(f"Length {length}: {len(records):,} sequences, {length_reads:,} reads")
        write_sequences_table(records, length_path)


def load_sequences_for_length(
    length_file: Path, sequence_map: dict[str, UniqueSequence]
) -> list[UniqueSequence]:
    """Return the UniqueSequence records referenced in a per-length file."""
    records: list[UniqueSequence] = []
    with length_file.open("r", encoding="ascii") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if row[0] == "sequence":
                continue
            sequence = row[0]
            try:
                records.append(sequence_map[sequence])
            except KeyError as exc:
                raise KeyError(
                    f"Sequence {sequence!r} from {length_file} missing in unique table"
                ) from exc
    return records


def connect_sequences_same_length(
    sequences: Sequence[UniqueSequence], n_edits: int
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
    return edges


def connect_sequences_different_length(
    sequences_a: Sequence[UniqueSequence],
    sequences_b: Sequence[UniqueSequence],
    n_edits: int,
) -> list[tuple[int, int]]:
    """Find all pairs of sequences within an edit distance for different-length sequences."""
    if not sequences_a or not sequences_b:
        return []
    len_a = len(sequences_a[0].sequence)
    len_b = len(sequences_b[0].sequence)
    if len_a > len_b:
        edges = connect_sequences_different_length(sequences_b, sequences_a, n_edits)
        return [(j, i) for i, j in edges]

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

    return edges


def deduplicate_edges(
    edges: Iterable[tuple[int, int]],
    sequences_a: Sequence[UniqueSequence],
    sequences_b: Sequence[UniqueSequence],
    same_length: bool,
) -> list[tuple[int, int]]:
    seen: set[tuple[int, int]] = set()
    result: list[tuple[int, int]] = []
    for i, j in edges:
        idx_a = sequences_a[i].index
        idx_b = sequences_b[j].index
        if same_length:
            if idx_a == idx_b:
                continue
            pair = (idx_a, idx_b) if idx_a < idx_b else (idx_b, idx_a)
        else:
            pair = (idx_a, idx_b)
        if pair in seen:
            continue
        seen.add(pair)
        result.append(pair)
    result.sort()
    return result


def write_edge_file(
    output_dir: Path,
    length_a: int,
    length_b: int,
    distance: int,
    edges: Sequence[tuple[int, int]],
) -> Path:
    """Write an edge list to a CSV file."""
    output_dir.mkdir(parents=True, exist_ok=True)
    name = f"pairs_len{min(length_a, length_b)}_len{max(length_a, length_b)}_d{distance}.csv"
    path = output_dir / name
    with path.open("w", newline="", encoding="ascii") as handle:
        handle.write(
            f"# length_a={length_a}\tlength_b={length_b}\tdistance={distance}\n"
        )
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["index_a", "index_b"])
        for index_a, index_b in edges:
            writer.writerow([str(index_a), str(index_b)])
    return path


def read_edge_file(path: Path) -> list[tuple[int, int]]:
    """Read an edge list from a CSV file."""
    edges: list[tuple[int, int]] = []
    with path.open("r", encoding="ascii") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "index_a":
                continue
            edges.append((int(row[0]), int(row[1])))
    return edges


def run_unique(args) -> None:
    """Extract unique sequences from a FASTQ file."""
    fastq_path = Path(args.fastq)
    output_path = Path(args.output)
    sequences, total_reads, skipped = collect_unique_sequences(fastq_path)
    write_sequences_table(sequences, output_path)
    print(
        f"Found {len(sequences):,} unique sequences "
        f"({total_reads:,} total reads, {skipped:,} skipped)."
    )
    print(f"Wrote unique sequence table to {output_path}")


def run_split(args) -> None:
    """Split unique sequences into per-length tables."""
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    sequences = read_sequences_table(input_path)
    split_by_length(
        sequences,
        output_dir,
        args.min_length,
        args.max_length,
    )
    print(
        f"Split {len(sequences):,} sequences into per-length tables under {output_dir}"
    )


def run_pairs(args) -> None:
    """Find and write sequence pairs within a certain edit distance."""
    unique_path = Path(args.unique)
    length_dir = Path(args.length_dir)
    length_a = args.length_a
    length_b = args.length_b
    distance = args.distance

    sequences = read_sequences_table(unique_path)
    sequence_map = {record.sequence: record for record in sequences}

    file_a = length_dir / f"length_{length_a}.csv"
    file_b = length_dir / f"length_{length_b}.csv"

    if not file_a.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_a}")
    if not file_b.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_b}")

    seqs_a = load_sequences_for_length(file_a, sequence_map)
    seqs_b = load_sequences_for_length(file_b, sequence_map)

    same_length = length_a == length_b
    if same_length:
        raw_edges = connect_sequences_same_length(seqs_a, distance)
        sequences_b = seqs_a
    else:
        raw_edges = connect_sequences_different_length(seqs_a, seqs_b, distance)
        sequences_b = seqs_b
    edges = deduplicate_edges(raw_edges, seqs_a, sequences_b, same_length)

    output_dir = Path(args.output_dir)
    edge_path = write_edge_file(output_dir, length_a, length_b, distance, edges)
    print(
        f"Found {len(edges):,} pairs within distance {distance} "
        f"for lengths ({length_a}, {length_b})."
    )
    print(f"Wrote edge list to {edge_path}")


def run_cluster(args) -> None:
    """Assemble clusters from edge lists and write representatives."""
    unique_path = Path(args.unique)
    output_path = Path(args.output)
    sequences = read_sequences_table(unique_path)
    total_reads = sum(record.count for record in sequences)
    dsu = DisjointSetUnion(len(sequences))

    edges: list[tuple[int, int]] = []
    for edge_file in args.edges:
        edge_path = Path(edge_file)
        file_edges = read_edge_file(edge_path)
        for u, v in file_edges:
            dsu.union(u, v)
        edges.extend(file_edges)

    components = dsu.get_components()
    clusters: list[tuple[UniqueSequence, int, int]] = []
    for component in components:
        total_count = sum(sequences[idx].count for idx in component)
        representative_idx = max(component, key=lambda idx: sequences[idx].count)
        representative = sequences[representative_idx]
        clusters.append((representative, representative_idx, total_count))

    clusters.sort(key=lambda item: item[2], reverse=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "count", "length", "frequency"])
        for representative, _, total_count in clusters:
            frequency = total_count / total_reads if total_reads else 0.0
            writer.writerow(
                [
                    representative.sequence,
                    str(total_count),
                    str(representative.length),
                    f"{frequency:.12g}",
                ]
            )

    print(
        f"Processed {len(sequences):,} sequences with {len(edges):,} edges into "
        f"{len(clusters):,} clusters."
    )
    print(f"Wrote cluster representatives to {output_path}")
