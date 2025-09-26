from .utils import is_valid_sequence


class FastQReader:
    """
    A simple FASTQ file reader that yields (header, sequence) tuples.
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

                yield header_line, sequence
                header = header_line[1:] if header_line.startswith('@') else header_line
                yield header, sequence

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class FastAReader:
    """
    A simple FASTA file reader that yields (header, sequence) tuples.
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

                yield header_line, sequence
                header = header_line[1:] if header_line.startswith('>') else header_line
                yield header, sequence

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

    def write(self, header: str, sequence: str, total_count: int):
        """
        Write a sequence to the FASTA file with a header including the total count.
        """
        self.handle.write(f">{header} N:{total_count}\n{sequence}\n")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.handle.close()
