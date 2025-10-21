import csv
import random
import sys
import logging
from typing import Any, Sequence
from pathlib import Path
from collections import defaultdict

import zarr
import numpy as np

from .utils import is_valid_sequence
from .types import UniqueSequence


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


class FastQReader:
    """
    A simple FASTQ file reader that yields sequences.
    """
    def __init__(self, fastq_file: str):
        self.fastq_file = fastq_file
        self.total = 0
        self.skipped = 0

    def __iter__(self):
        with open(self.fastq_file, 'r', encoding='ascii') as f:
            while True:
                header_line = f.readline().strip()
                if not header_line:
                    break
                sequence = f.readline().strip()
                f.readline()  # Skip '+' line
                f.readline()  # Skip quality line

                if not is_valid_sequence(sequence):
                    self.skipped += 1
                    continue
                self.total += 1

                yield sequence

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class FastAReader:
    """
    A simple FASTA file reader that yields sequences.
    """
    def __init__(self, fastq_file: str):
        self.fastq_file = fastq_file
        self.total = 0
        self.skipped = 0

    def __iter__(self):
        with open(self.fastq_file, 'r', encoding='ascii') as f:
            while True:
                header_line = f.readline().strip()
                if not header_line:
                    break
                sequence = f.readline().strip()

                if not is_valid_sequence(sequence):
                    self.skipped += 1
                    continue
                self.total += 1

                yield sequence

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class FastAWriter:
    """
    A simple FASTA file writer.
    """
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        self.handle = open(fasta_file, 'w', encoding='ascii')
        self._idx = 0

    def write(self, sequence: str, total_count: int):
        """
        Write a sequence to the FASTA file with a header including the total count.
        """
        self._idx += 1
        self.handle.write(f">{self._idx} N:{total_count} len:{len(sequence)}\n{sequence}\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.handle.close()


class ZarrStoreByLength:
    """
    A class to manage Zarr stores for sequences grouped by their lengths.
    """
    def __init__(self, base_path: Path):
        self.base_path = base_path

    @staticmethod
    def write(
        sequences: Sequence[UniqueSequence],
        base_path: Path,
        chunk_size: int
    ) -> None:
        """Persist unique sequences to a Zarr store, grouped by sequence length."""
        base_path.parent.mkdir(parents=True, exist_ok=True)
        root = zarr.open_group(str(base_path), mode="w")

        # Collect sequences by length
        grouped: dict[int, list[UniqueSequence]] = defaultdict(list)
        for record in sequences:
            grouped[len(record.sequence)].append(record)
        sorted_grouped = dict(sorted(grouped.items()))

        # Write overall stats
        total_sequences = len(sequences)
        total_reads = sum(record.count for record in sequences)
        root.attrs["total_sequences"] = total_sequences
        root.attrs["total_reads"] = total_reads

        for length, records in sorted_grouped.items():
            group = root.create_group(f"length_{length}", overwrite=True)

            # Randomize order to avoid similarity clusters
            # (-> better load balancing in pairwise comparisons later)
            random.shuffle(records)

            # Write stats for this length
            length_reads = sum(r.count for r in records)
            group.attrs["unique_sequences"] = len(records)
            group.attrs["total_reads"] = length_reads
            group.attrs["sequence_length"] = length
            logger.info(
                "Length %d: %s sequences, %s reads",
                length,
                format(len(records), ","),
                format(length_reads, ","),
            )

            # Write sequences and counts as zarr arrays
            str_type = f"<U{length}"
            sequences_arr = np.array([r.sequence for r in records], dtype=str_type)
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

    def load_sequences(
        self,
        length: int,
        idx: Any = slice(None)
    ) -> list[str]:
        """Load unique sequences of a specific length from the Zarr store."""
        root = zarr.open_group(str(self.base_path), mode="r")
        group_name = f"length_{length}"
        if group_name not in root:
            raise ValueError(f"No sequences of length {length} found in store {self.base_path}")

        group = root[group_name]
        sequences_arr = group["sequence"][idx]

        return [str(s) for s in sequences_arr]

    def load_counts(
        self,
        length: int,
        idx: Any = slice(None)
    ) -> list[int]:
        """Load counts of unique sequences of a specific length from the Zarr store."""
        root = zarr.open_group(str(self.base_path), mode="r")
        group_name = f"length_{length}"
        if group_name not in root:
            raise ValueError(f"No sequences of length {length} found in store {self.base_path}")

        group = root[group_name]
        return group["count"][idx]


def write_sequences_table(
    sequences: Sequence[UniqueSequence], output_path: Path
) -> None:
    """Persist unique sequences to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        # Write the header
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "count"])

        # Write data rows
        for record in sequences:
            writer.writerow([record.sequence, str(record.count)])


def read_sequences_table(
    path: Path,
    sequence_column: str,
    count_column: str,
) -> list[UniqueSequence]:
    """Load unique sequences from a CSV file."""
    sequences: list[UniqueSequence] = []
    with path.open("r", encoding="ascii") as handle:
        # Skip comment lines
        while True:
            pos = handle.tell()
            line = handle.readline()
            if not line.startswith("#"):
                handle.seek(pos)
                break

        # Detect dialect (in particular, the delimiter)
        try:
            header_line = handle.readline()
            dialect = csv.Sniffer().sniff(header_line)
        except csv.Error as exc:
            raise ValueError(f"Unable to detect delimiter in {path}") from exc

        # Read the header
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")
        expected = {sequence_column, count_column}
        if not expected.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Missing required columns {sequence_column}, {count_column} in {path}; "
                f"available: {reader.fieldnames}"
            )

        # Read data rows
        for row in reader:
            sequence = row["sequence"].strip()
            count = int(row["count"].strip())
            sequences.append(
                UniqueSequence(sequence=sequence, count=count)
            )

    return sequences
