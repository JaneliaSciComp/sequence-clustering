"""Validate overlap and reactivity agreement between two clustering pipelines.

This script compares the top N sequences produced by two FASTA outputs and
reports:
- the fraction of sequences shared among the top N, and
- the mean absolute error between their normalized reactivities.

The FASTA headers are expected to contain a trailing token of the form
``N:<value>`` which is interpreted as the (unnormalized) reactivity/count for the
sequence.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

HEADER_COUNT_PATTERN = re.compile(r"N:(-?\d+(?:\.\d+)?)(?!.*N:)")


def parse_fasta_counts(path: Path) -> dict[str, float]:
    """Return a mapping of sequence -> raw count/reactivity extracted from FASTA."""
    counts: dict[str, float] = {}
    with path.open("r", encoding="ascii") as handle:
        while True:
            header = handle.readline().strip()
            if not header:
                break
            sequence = handle.readline().strip()

            match = HEADER_COUNT_PATTERN.search(header)
            if not match:
                raise ValueError(
                    f"Could not find 'N:<value>' token in header: {header}"
                )
            counts[sequence] = int(match.group(1))

    return counts

def compute_metrics(
    counts_a: dict[str, float],
    counts_b: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute fraction of shared top sequences and cdfs of counts."""
    counts_a = dict(sorted(counts_a.items(), key=lambda item: item[1], reverse=True))
    counts_b = dict(sorted(counts_b.items(), key=lambda item: item[1], reverse=True))

    # Compute cumulative distribution functions of counts
    cdf_a = np.array(list(counts_a.values()))
    cdf_a = np.cumsum(cdf_a) / np.sum(cdf_a)
    cdf_b = np.array(list(counts_b.values()))
    cdf_b = np.cumsum(cdf_b) / np.sum(cdf_b)

    # Compute fraction of shared sequences in top N
    unique_sequences = set()
    fraction_shared = []
    for i, (seq_a, seq_b) in enumerate(zip(counts_a.keys(), counts_b.keys()), start=1):
        unique_sequences.add(seq_a)
        unique_sequences.add(seq_b)
        fraction_shared.append(2 - len(unique_sequences) / i)
    fraction_shared = np.array(fraction_shared)

    return fraction_shared, cdf_a, cdf_b


def format_report(
    fraction_shared: float,
) -> str:
    """Format the final report for stdout."""
    line_idx = 10 ** np.arange(0, np.log10(len(fraction_shared)), dtype=int) - 1
    lines = []
    for idx in line_idx:
        lines.append(f"Top {idx + 1:>7,}: {fraction_shared[idx] * 100:>4.2f}% shared")

    return "\n".join(lines)


def plot_histograms(
    fraction_shared: np.ndarray,
    cdf_a: np.ndarray,
    cdf_b: np.ndarray,
) -> None:
    """Plot histograms of fraction shared and cdfs."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Pipeline Agreement Validation", fontsize=16)

    # Left plot: Fraction shared
    ax1.plot(fraction_shared * 100, color="tab:blue")
    ax1.set_xscale("log")
    ax1.set_xlabel("Top N sequences")
    ax1.set_ylabel("Fraction shared (%)")
    ax1.tick_params(axis="y")
    ax1.grid(True, which="both", linestyle="--", linewidth=0.5)

    # Right plot: Cumulative distributions
    ax2.plot(cdf_a, label="Old pipeline", color="tab:orange", linestyle="--")
    ax2.plot(cdf_b, label="New pipeline", color="tab:green", linestyle="-")
    ax2.set_xscale("log")
    ax2.set_xlabel("Top N sequences")
    ax2.set_ylabel("Cumulative fraction of counts")
    ax2.legend(loc="lower right")
    ax2.grid(True, which="both", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.show()


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Validate pipeline agreement")
    parser.add_argument("old", type=Path, help="FASTA output from the legacy pipeline")
    parser.add_argument("new", type=Path, help="FASTA output from the new pipeline")
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    counts_old = parse_fasta_counts(args.old)
    counts_new = parse_fasta_counts(args.new)

    print(f"Parsed {len(counts_old):,} clusters from old pipeline")
    print(f"Parsed {len(counts_new):,} clusters from new pipeline")

    fraction_shared, new_cdf, old_cdf = compute_metrics(counts_old, counts_new)

    print(format_report(fraction_shared))
    plot_histograms(fraction_shared, new_cdf, old_cdf)


if __name__ == "__main__":
    main()
