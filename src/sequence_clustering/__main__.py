import argparse

from .cmd import run_cluster, run_pairs, run_split, run_unique


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the sequence clustering pipeline."""
    parser = argparse.ArgumentParser(
        description="Modular sequence clustering pipeline",
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

    # Subcommand: split
    split_parser = subparsers.add_parser(
        "split", help="Split unique table into per-length CSVs"
    )
    split_parser.add_argument(
        "--input", "-i", required=True, help="Unique CSV from step 1"
    )
    split_parser.add_argument(
        "--output-dir", "-o", required=True, help="Directory for per-length CSV files"
    )
    split_parser.set_defaults(func=run_split)

    # Subcommand: pairs
    pairs_parser = subparsers.add_parser(
        "pairs", help="Find sequence pairs within an edit distance"
    )
    pairs_parser.add_argument(
        "--length-dir",
        required=True,
        help="Directory containing per-length CSV files",
    )
    pairs_parser.add_argument(
        "--length-a", type=int, required=True, help="First sequence length"
    )
    pairs_parser.add_argument(
        "--length-b", type=int, required=True, help="Second sequence length"
    )
    pairs_parser.add_argument(
        "--distance", type=int, required=True, help="Maximum edit distance"
    )
    pairs_parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write the edge list",
    )
    pairs_parser.set_defaults(func=run_pairs)

    # Subcommand: cluster
    cluster_parser = subparsers.add_parser(
        "cluster", help="Assemble clusters from edge lists"
    )
    cluster_parser.add_argument(
        "--unique", required=True, help="Unique CSV from step 1"
    )
    cluster_parser.add_argument(
        "--edges",
        nargs="*",
        default=[],
        help="Edge list files produced by the pairs subcommand",
    )
    cluster_parser.add_argument(
        "--output", "-o", required=True, help="Output CSV for cluster representatives"
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
