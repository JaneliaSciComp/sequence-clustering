# 25N Variable Region Extraction Script

This script (`fastp_extract_25n.sh`) trims and filters paired-end FASTQ data
using **fastp**, then extracts highly enriched ~25-nucleotide variable regions
and produces sequence count statistics.

---

## Usage
```bash
./fastp_extract_25n.sh <input_R1.fastq.gz> <input_R2.fastq.gz> <output_prefix>
```

| Argument | Description |
|----------|-------------|
| `input_R1.fastq.gz` | Read 1 FASTQ file (gzipped) |
| `input_R2.fastq.gz` | Read 2 FASTQ file (gzipped) |
| `output_prefix`     | Prefix for all output files |

Example:
```bash
./fastp_extract_25n.sh sample_R1.fq.gz sample_R2.fq.gz sample_out
```

---

## Requirements
Install or have available in your `$PATH`:
- **fastp** (adapter detection, quality trimming)
- Standard UNIX utilities: `bash`, `gzip`/`zcat`, `awk`, `sort`, `uniq`, `bc`,
  `sed`, `rev`, `tr`

You can install fastp via conda:
```bash
conda install -c bioconda fastp
```

---

## Key Output Files
All outputs share the `<output_prefix>` you provide:

| File | Description |
|------|-------------|
| `<prefix>_fastp.html` | Interactive QC report from fastp |
| `<prefix>_fastp.json` | Machine-readable QC summary |
| `<prefix>_counts.txt` | Tab-delimited table of unique sequences: `sequence`, `count`, `length`, `frequency` |
| `<prefix>_sequences.fasta` | FASTA of unique sequences (headers include counts and length) |
| `<prefix>_summary.txt` | Top 20 most abundant sequences |
| `<prefix>_length_distribution.txt` | Length vs. unique count and total read count |

Temporary intermediate FASTQ and raw count files are removed automatically.

---

## Pipeline Overview
1. **Adapter trimming & QC** – `fastp` automatically detects adapters, corrects
   bases, trims low-quality regions, and filters by length and complexity.
2. **Length analysis** – Counts reads in the 20–30 bp range to identify which
   mate (R1 or R2) carries the 25 N variable region.
3. **Sequence extraction** – Chooses the higher-signal mate (or both), reverse-complements R2 if used, and collapses identical sequences.
4. **Counting & reporting** – Generates count tables, FASTA, summary, and length-distribution statistics.

---

## Notes
- Inspect the `*_fastp.html` report and `*_length_distribution.txt`
  to confirm expected ~25 bp enrichment before downstream analysis.
- The script automatically cleans up large intermediate files to save space.
