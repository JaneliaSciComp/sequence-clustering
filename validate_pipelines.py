"""Validate overlap and reactivity agreement between two clustering pipelines.

This script compares the top N sequences produced by two FASTA outputs and
reports:
- the fraction of sequences shared among the top N, and
- the mean absolute error between their normalized reactivities.

The FASTA headers are expected to contain a trailing token of the form
``N:<value>`` which is interpreted as the (unnormalized) reactivity/count for the
sequence.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, Tuple

HEADER_COUNT_PATTERN = re.compile(r"N:(-?\d+(?:\.\d+)?)")


def parse_fasta_counts(path: Path) -> Dict[str, float]:
    """Return a mapping of sequence -> raw count/reactivity extracted from FASTA."""
    counts: Dict[str, float] = {}
    with path.open("r", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            if not sequence:
                raise ValueError(f"Sequence missing after header in {path}")

            header = header.strip()
            sequence = sequence.strip()

            match = HEADER_COUNT_PATTERN.search(header)
            if not match:
                raise ValueError(
                    f"Could not find 'N:<value>' token in header: {header}"
                )
            counts[sequence] = float(match.group(1))
    return counts


def normalize_counts(counts: Dict[str, float]) -> Dict[str, float]:
    """Normalize counts so they sum to 1 and clip the result to [0, 1]."""
    total = sum(counts.values())
    if total <= 0:
        return {sequence: 0.0 for sequence in counts}

    normalized = {
        sequence: clip_between_0_and_1(value / total)
        for sequence, value in counts.items()
    }
    return normalized


def clip_between_0_and_1(value: float) -> float:
    """Clip a floating point value to the inclusive range [0, 1]."""
    if value < 0:
        return 0.0
    if value > 1:
        return 1.0
    return value


def top_sequences(counts: Dict[str, float], top_n: int) -> Iterable[str]:
    """Return the identifiers of the top_n sequences ordered by descending count."""
    return [seq for seq, _ in sorted(counts.items(), key=lambda item: item[1], reverse=True)[:top_n]]


def compute_metrics(
    counts_a: Dict[str, float],
    counts_b: Dict[str, float],
    top_n: int,
) -> Tuple[float, float]:
    """
    Compute fraction of shared top sequences and MAE of normalized reactivities.

    Returns a tuple ``(fraction_shared, mae)``. The MAE is ``float('nan')`` if
    there are no shared sequences in the top N.
    """
    top_a = top_sequences(counts_a, top_n)
    top_b = top_sequences(counts_b, top_n)

    set_a = set(top_a)
    set_b = set(top_b)
    shared_sequences = sorted(set_a & set_b)

    if not shared_sequences:
        return 0.0, float("nan")

    denominator = min(top_n, len(top_a), len(top_b))
    fraction_shared = len(shared_sequences) / denominator if denominator else 0.0

    norm_a = normalize_counts(counts_a)
    norm_b = normalize_counts(counts_b)

    absolute_errors = [
        abs(norm_a[sequence] - norm_b[sequence]) for sequence in shared_sequences
    ]
    mae = sum(absolute_errors) / len(absolute_errors)

    return fraction_shared, mae


def format_report(
    fraction_shared: float,
    mae: float,
    top_n: int,
    shared_count: int,
) -> str:
    """Format the final report for stdout."""
    lines = [
        f"Top {top_n} overlap fraction: {fraction_shared:.3f}",
    ]
    if shared_count:
        lines.append(
            f"Mean absolute error of normalized reactivities (shared {shared_count}): {mae:.3f}"
        )
    else:
        lines.append("Mean absolute error of normalized reactivities: n/a (no overlap)")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate pipeline agreement")
    parser.add_argument("old", type=Path, help="FASTA output from the legacy pipeline")
    parser.add_argument("new", type=Path, help="FASTA output from the new pipeline")
    parser.add_argument(
        "--top",
        type=int,
        default=100,
        help="Number of top sequences to compare (default: 100)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts_old = parse_fasta_counts(args.old)
    counts_new = parse_fasta_counts(args.new)

    fraction_shared, mae = compute_metrics(counts_old, counts_new, args.top)

    top_a = top_sequences(counts_old, args.top)
    top_b = top_sequences(counts_new, args.top)
    shared_count = len(set(top_a) & set(top_b))

    print(format_report(fraction_shared, mae, args.top, shared_count))


if __name__ == "__main__":
    main()
