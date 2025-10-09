import concurrent.futures
import subprocess
import sys
import csv
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Sequence

from .dsu import DisjointSetUnion
from .types import UniqueSequence
from .io import (
    read_sequences_table,
    write_sequences_table,
    read_edge_file,
    write_edge_file,
    extract_counts,
    FastQReader
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


def split_by_length(
    sequences: Sequence[UniqueSequence],
    output_dir: Path,
) -> None:
    """Write per-length tables with summary headers."""
    # Group sequences by length
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[int, list[UniqueSequence]] = defaultdict(list)
    for record in sequences:
        grouped[len(record.sequence)].append(record)
    grouped = dict(sorted(grouped.items()))

    # Write per-length files
    for length, records in grouped.items():
        length_reads = sum(r.count for r in records)
        length_path = output_dir / f"length_{length}.csv"
        print(f"Length {length}: {len(records):,} sequences, {length_reads:,} reads")
        write_sequences_table(records, length_path)


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
    output_dir = Path(args.output_dir)
    sequences = read_sequences_table(input_path)
    split_by_length(sequences, output_dir)
    print(
        f"Split {len(sequences):,} sequences into per-length tables under {output_dir}"
    )
    print(f"Time elapsed: {time.time() - start:.2g} seconds")


def discover_length_files(length_dir: Path) -> dict[int, Path]:
    """Return mapping of sequence length to per-length CSV path."""
    pattern = re.compile(r"length_(\d+)\.csv")
    result: dict[int, Path] = {}
    for candidate in sorted(length_dir.glob("length_*.csv")):
        match = pattern.fullmatch(candidate.name)
        if not match:
            raise ValueError(f"Unexpected per-length filename: {candidate}")
        length = int(match.group(1))
        result[length] = candidate
    return result


def run_pairs(args) -> None:
    """Find and write sequence pairs within a certain edit distance."""
    start = time.time()
    length_dir = Path(args.length_dir)
    length_a = min(args.length_a, args.length_b)
    length_b = max(args.length_a, args.length_b)
    distance = args.distance

    # Read headers of all files to compute offsets
    length_to_n = {}
    for file in length_dir.glob("length_*.csv"):
        length = int(re.match(r"length_(\d+)\.csv", file.name).group(1))
        n_sequences, _ = extract_counts(file)
        length_to_n[length] = n_sequences
    offset_a = sum(
        n for length, n in length_to_n.items() if length < length_a
    )
    offset_b = sum(
        n for length, n in length_to_n.items() if length < length_b
    )

    # Load all sequences with given lengths
    file_a = length_dir / f"length_{length_a}.csv"
    file_b = length_dir / f"length_{length_b}.csv"

    if not file_a.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_a}")
    if not file_b.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_b}")

    seqs_a = read_sequences_table(file_a)
    seqs_b = read_sequences_table(file_b)

    same_length = length_a == length_b
    if same_length:
        raw_edges = connect_sequences_same_length(seqs_a, distance, offset_a)
    else:
        raw_edges = connect_sequences_different_length(
            seqs_a, seqs_b, distance, offset_a, offset_b
        )
    edges = list(set(raw_edges))  # Remove duplicates

    output_path = Path(args.output_dir) / f"pairs_len{length_a}_len{length_b}_d{distance}.csv"
    write_edge_file(output_path, edges)
    print(
        f"({length_a}, {length_b}): Found {len(edges):,} pairs within distance {distance} "
    )
    print(f"({length_a}, {length_b}): Wrote edge list to {output_path}")
    print(f"({length_a}, {length_b}): Time elapsed: {time.time() - start:.2g} seconds")


def generate_length_pairs(lengths: list[int], max_distance: int) -> list[tuple[int, int]]:
    """Return all length pairs (a <= b) within the given distance."""
    pairs: list[tuple[int, int]] = []
    for i, a in enumerate(lengths):
        for b in lengths[i:]:
            if abs(a - b) <= max_distance:
                pairs.append((a, b))
    return pairs


def run_all_pairs(args) -> None:
    """Run pairs for all length combinations within an edit distance."""
    start = time.time()
    length_dir = Path(args.length_dir)
    output_dir = Path(args.output_dir)
    max_workers = max(args.workers, 1)

    length_files = discover_length_files(length_dir)
    lengths = sorted(length_files)
    pairs = generate_length_pairs(lengths, args.distance)
    if not pairs:
        print("No eligible length pairs found.")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    def launch_pair(length_a: int, length_b: int) -> None:
        cmd = [
            sys.executable, "-m",
            "sequence_clustering", "pairs",
            "--length-dir", str(length_dir),
            "--length-a", str(length_a),
            "--length-b", str(length_b),
            "--distance", str(args.distance),
            "--output-dir", str(output_dir),
        ]
        subprocess.run(cmd, check=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_pair = {
            executor.submit(launch_pair, a, b): (a, b) for a, b in pairs
        }
        for future in concurrent.futures.as_completed(future_to_pair):
            a, b = future_to_pair[future]
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001
                raise RuntimeError(
                    f"pairs command failed for lengths ({a}, {b})"
                ) from exc

    print(f"Processed {len(pairs):,} length pairs with {max_workers} workers.")
    print(f"Time elapsed: {time.time() - start:.2g} seconds")


def run_cluster(args) -> None:
    """Assemble clusters from edge lists and write representatives."""
    start = time.time()
    unique_path = Path(args.unique)
    edges_path = Path(args.edges_dir)
    output_path = Path(args.output)
    sequences = read_sequences_table(unique_path)
    dsu = DisjointSetUnion(len(sequences))

    # Read all edge files and union connected sequences
    n_edges = 0
    for edge_file in edges_path.glob("*.csv"):
        edge_path = Path(edge_file)
        file_edges = read_edge_file(edge_path)
        for u, v in file_edges:
            dsu.union(u, v)
        n_edges += len(file_edges)

    # Generate clusters (connected components of the graph)
    # Choose most abundant sequence as representative
    components = dsu.get_components()
    clusters: list[tuple[str, int, int]] = []
    for component in components:
        total_count = sum(sequences[idx].count for idx in component)
        representative_idx = max(component, key=lambda idx: sequences[idx].count)
        representative = sequences[representative_idx].sequence
        clusters.append((representative, len(component), total_count))

    # Sort clusters by total count (descending)
    clusters.sort(key=lambda item: item[2], reverse=True)

    # Write cluster representatives to output file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "cluster_size", "total_count"])
        for representative, cluster_size, total_count in clusters:
            writer.writerow([representative, str(cluster_size), str(total_count)])

    print(
        f"Processed {len(sequences):,} sequences with {n_edges:,} edges into "
        f"{len(clusters):,} clusters."
    )
    print(f"Wrote cluster representatives to {output_path}")
    print(f"Time elapsed: {time.time() - start:.2g} seconds")
