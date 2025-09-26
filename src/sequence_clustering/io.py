import csv
from typing import Sequence
from pathlib import Path

from .utils import is_valid_sequence
from .types import UniqueSequence


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


def write_sequences_table(
    sequences: Sequence[UniqueSequence], output_path: Path
) -> None:
    """Persist unique sequences to a tab-delimited file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="ascii") as handle:
        # Write a comment line with summary statistics
        unique_sequences = len(sequences)
        total_reads = sum(record.count for record in sequences)
        handle.write(f"# {unique_sequences=}, {total_reads=}\n")

        # Write the header
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence", "count"])

        # Write data rows
        for record in sequences:
            writer.writerow([record.sequence, str(record.count)])


def read_sequences_table(path: Path) -> tuple[list[UniqueSequence], int]:
    """Load unique sequences written by step 1 and return records plus total reads."""
    sequences: list[UniqueSequence] = []
    total_reads = 0
    with path.open("r", encoding="ascii") as handle:
        # Skip comment lines
        handle.readline()

        # Read the header
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")
        expected = {"sequence", "count"}
        if set(reader.fieldnames) != expected:
            raise ValueError(
                f"Unexpected columns in {path}: {reader.fieldnames}"
            )

        # Read data rows
        for row in reader:
            sequence = row["sequence"].strip()
            count = int(row["count"])
            total_reads += count
            sequences.append(
                UniqueSequence(sequence=sequence, count=count)
            )

    return sequences, total_reads
