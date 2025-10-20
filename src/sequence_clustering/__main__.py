import argparse

from .cmd import run_cluster, run_unique


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the sequence clustering pipeline."""
    parser = argparse.ArgumentParser(
        description="Sequence clustering pipeline",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: unique
    unique_parser = subparsers.add_parser(
        "unique", help="Extract unique sequences and counts from FASTQ"
    )
    unique_parser.add_argument(
        "--fastq", "-i", required=True, help="Input FASTQ file"
    )
    unique_parser.add_argument(
        "--output", "-o", required=True, help="Output CSV for unique sequences"
    )
    unique_parser.set_defaults(func=run_unique)

    # Subcommand: cluster
    cluster_parser = subparsers.add_parser(
        "cluster", help="Compute edges with Dask and assemble clusters"
    )
    cluster_parser.add_argument(
        "--input", "-i", required=True, help="CSV file with global unique sequences"
    )
    cluster_parser.add_argument(
        "--length-store",
        help="Optional override for the per-length Zarr store (default: by_length.zarr)",
    )
    cluster_parser.add_argument(
        "--distance", "-d",
        type=int,
        required=True,
        help="Maximum edit distance",
    )
    cluster_parser.add_argument(
        "--output", "-o", required=True, help="Output CSV for cluster representatives"
    )
    cluster_parser.add_argument(
        "--sequence-column",
        default="sequence",
        help="Column name containing sequence strings (default: sequence)",
    )
    cluster_parser.add_argument(
        "--count-column",
        default="count",
        help="Column name containing sequence counts (default: count)",
    )
    cluster_parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Chunk size to use when writing Zarr arrays (default: 10000)",
    )
    cluster_parser.add_argument(
        "--tile-size",
        type=int,
        default=10000,
        help="Tile size to do pairwise comparisons on (default: 10000)",
    )
    cluster_parser.add_argument(
        "--workers", "-w",
        type=int,
        default=0,
        help="Number of Dask workers to launch (default: auto)",
    )
    cluster_parser.add_argument(
        "--threads-per-worker",
        type=int,
        default=0,
        help="Threads per Dask worker (default: auto)",
    )
    cluster_parser.set_defaults(func=run_cluster)

    return parser


def main() -> None:
    """Entry point for the sequence clustering pipeline."""
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
