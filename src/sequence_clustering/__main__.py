import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from .dsu import DisjointSetUnion
from .types import UniqueSequence
from .utils import (
    compare_buckets,
    fill_buckets,
    generate_partitions,
    is_valid_sequence,
)


def collect_unique_sequences(fastq_path: Path) -> tuple[list[UniqueSequence], int, int]:
    """Return unique sequences, total reads, and skipped reads from a FASTQ file."""
    counts: dict[str, int] = {}
    total_reads = 0
    skipped_reads = 0

    with fastq_path.open("r", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise ValueError("Unexpected FASTQ structure: incomplete record detected")

            sequence = sequence.strip()
            if not is_valid_sequence(sequence):
                skipped_reads += 1
                continue

            total_reads += 1
            counts[sequence] = counts.get(sequence, 0) + 1

    ordered = sorted(counts.items(), key=lambda item: item[1], reverse=True)
    uniques = [
        UniqueSequence(sequence=seq, count=count, length=len(seq), index=idx)
        for idx, (seq, count) in enumerate(ordered)
    ]
    return uniques, total_reads, skipped_reads


def write_unique_sequences_table(
    sequences: Sequence[UniqueSequence], total_reads: int, output_path: Path
) -> None:
    """Persist unique sequences to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "count", "length", "frequency"])
        for record in sequences:
            frequency = record.count / total_reads if total_reads else 0.0
            writer.writerow(
                [
                    record.sequence,
                    str(record.count),
                    str(record.length),
                    f"{frequency:.12g}",
                ]
            )


def read_unique_sequences_table(path: Path) -> tuple[list[UniqueSequence], int]:
    """Load unique sequences written by step 1 and return records plus total reads."""
    sequences: list[UniqueSequence] = []
    total_reads = 0
    with path.open("r", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")
        expected = {"sequence", "count", "length", "frequency"}
        if set(reader.fieldnames) != expected:
            raise ValueError(
                f"Unexpected columns in {path}: {reader.fieldnames}"
            )
        for index, row in enumerate(reader):
            sequence = row["sequence"].strip()
            count = int(row["count"])
            length = int(row["length"])
            total_reads += count
            sequences.append(
                UniqueSequence(sequence=sequence, count=count, length=length, index=index)
            )
    return sequences, total_reads


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def split_by_length(
    sequences: Sequence[UniqueSequence],
    total_reads: int,
    output_dir: Path,
    min_length: int | None,
    max_length: int | None,
) -> None:
    """Write per-length tables with summary headers."""
    ensure_directory(output_dir)
    grouped: dict[int, list[UniqueSequence]] = defaultdict(list)
    for record in sequences:
        if min_length is not None and record.length < min_length:
            continue
        if max_length is not None and record.length > max_length:
            continue
        grouped[record.length].append(record)

    for length, records in grouped.items():
        length_path = output_dir / f"length_{length}.tsv"
        length_reads = sum(r.count for r in records)
        with length_path.open("w", newline="", encoding="ascii") as handle:
            handle.write(
                f"# unique_sequences={len(records)}\treads={length_reads}\n"
            )
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sequence", "count", "length", "frequency"])
            for record in records:
                frequency = record.count / total_reads if total_reads else 0.0
                writer.writerow(
                    [
                        record.sequence,
                        str(record.count),
                        str(record.length),
                        f"{frequency:.12g}",
                    ]
                )


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
            compare_buckets(list(bucket_a), list(bucket_b), sequences_a, sequences_b, n_edits, edges)

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
    ensure_directory(output_dir)
    name = f"pairs_len{min(length_a, length_b)}_len{max(length_a, length_b)}_d{distance}.tsv"
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


def run_unique(args: argparse.Namespace) -> None:
    fastq_path = Path(args.fastq)
    output_path = Path(args.output)
    sequences, total_reads, skipped = collect_unique_sequences(fastq_path)
    write_unique_sequences_table(sequences, total_reads, output_path)
    print(
        f"Found {len(sequences):,} unique sequences "
        f"({total_reads:,} total reads, {skipped:,} skipped)."
    )
    print(f"Wrote unique sequence table to {output_path}")


def run_split(args: argparse.Namespace) -> None:
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    sequences, total_reads = read_unique_sequences_table(input_path)
    split_by_length(
        sequences,
        total_reads,
        output_dir,
        args.min_length,
        args.max_length,
    )
    print(
        f"Split {len(sequences):,} sequences into per-length tables under {output_dir}"
    )


def run_pairs(args: argparse.Namespace) -> None:
    unique_path = Path(args.unique)
    length_dir = Path(args.length_dir)
    length_a = args.length_a
    length_b = args.length_b
    distance = args.distance

    sequences, _ = read_unique_sequences_table(unique_path)
    sequence_map = {record.sequence: record for record in sequences}

    file_a = length_dir / f"length_{length_a}.tsv"
    file_b = length_dir / f"length_{length_b}.tsv"

    if not file_a.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_a}")
    if not file_b.exists():
        raise FileNotFoundError(f"Missing per-length file: {file_b}")

    seqs_a = load_sequences_for_length(file_a, sequence_map)
    seqs_b = load_sequences_for_length(file_b, sequence_map)

    same_length = length_a == length_b
    if same_length:
        raw_edges = connect_sequences_same_length(seqs_a, distance)
    else:
        raw_edges = connect_sequences_different_length(seqs_a, seqs_b, distance)
    edges = deduplicate_edges(raw_edges, seqs_a, seqs_b if not same_length else seqs_a, same_length)

    output_dir = Path(args.output_dir)
    edge_path = write_edge_file(output_dir, length_a, length_b, distance, edges)
    print(
        f"Found {len(edges):,} pairs within distance {distance} "
        f"for lengths ({length_a}, {length_b})."
    )
    print(f"Wrote edge list to {edge_path}")


def run_cluster(args: argparse.Namespace) -> None:
    unique_path = Path(args.unique)
    output_path = Path(args.output)
    sequences, total_reads = read_unique_sequences_table(unique_path)
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sequence_clustering",
        description="Modular sequence clustering pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    unique_parser = subparsers.add_parser(
        "unique", help="Extract unique sequences and counts from FASTQ"
    )
    unique_parser.add_argument("--fastq", "-i", required=True, help="Input FASTQ file")
    unique_parser.add_argument(
        "--output", "-o", required=True, help="Output TSV for unique sequences"
    )
    unique_parser.set_defaults(func=run_unique)

    split_parser = subparsers.add_parser(
        "split", help="Split unique table into per-length TSVs"
    )
    split_parser.add_argument("--input", "-i", required=True, help="Unique TSV from step 1")
    split_parser.add_argument(
        "--output-dir", "-o", required=True, help="Directory for per-length TSV files"
    )
    split_parser.add_argument(
        "--min-length", type=int, default=None, help="Discard sequences shorter than this"
    )
    split_parser.add_argument(
        "--max-length", type=int, default=None, help="Discard sequences longer than this"
    )
    split_parser.set_defaults(func=run_split)

    pairs_parser = subparsers.add_parser(
        "pairs", help="Find sequence pairs within an edit distance"
    )
    pairs_parser.add_argument("--unique", required=True, help="Unique TSV from step 1")
    pairs_parser.add_argument(
        "--length-dir",
        required=True,
        help="Directory containing per-length TSV files from step 2",
    )
    pairs_parser.add_argument("--length-a", type=int, required=True, help="First sequence length")
    pairs_parser.add_argument("--length-b", type=int, required=True, help="Second sequence length")
    pairs_parser.add_argument("--distance", type=int, required=True, help="Maximum edit distance")
    pairs_parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write the edge list",
    )
    pairs_parser.set_defaults(func=run_pairs)

    cluster_parser = subparsers.add_parser(
        "cluster", help="Assemble clusters from edge lists"
    )
    cluster_parser.add_argument("--unique", required=True, help="Unique TSV from step 1")
    cluster_parser.add_argument(
        "--edges",
        nargs="*",
        default=[],
        help="Edge list files produced by the pairs subcommand",
    )
    cluster_parser.add_argument(
        "--output", "-o", required=True, help="Output TSV for cluster representatives"
    )
    cluster_parser.set_defaults(func=run_cluster)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
