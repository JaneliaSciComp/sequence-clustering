#!/bin/bash

# 25N Variable Region Extraction Pipeline using fastp
# Usage: ./fastp_extract_25n.sh input_R1.fq.gz input_R2.fq.gz output_prefix

set -euo pipefail

INPUT_R1=$1
INPUT_R2=$2
OUTPUT_PREFIX=$3

echo "Starting 25N extraction pipeline using fastp..."
echo "Input R1: $INPUT_R1"
echo "Input R2: $INPUT_R2"
echo "Output prefix: $OUTPUT_PREFIX"

# Step 1: Use fastp for automatic adapter detection and trimming
echo "Step 1: Running fastp for adapter detection and trimming..."

fastp \
    --in1 "$INPUT_R1" \
    --in2 "$INPUT_R2" \
    --out1 "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" \
    --out2 "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --correction \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 15 \
    --length_limit 50 \
    --low_complexity_filter \
    --complexity_threshold 30 \
    --thread 16 \
    --json "${OUTPUT_PREFIX}_fastp.json" \
    --html "${OUTPUT_PREFIX}_fastp.html" \
    --report_title "${OUTPUT_PREFIX} 25N extraction"

echo "Step 2: Analyzing fastp results..."

# Check if fastp found our expected structure by looking at length distribution
echo "Length distribution after fastp trimming:"

# Get length distribution from both R1 and R2
{
    zcat "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" | awk 'NR%4==2 {print length($0)}'
    zcat "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" | awk 'NR%4==2 {print length($0)}'
} | sort -n | uniq -c | sort -nr | head -10

# Step 3: Extract sequences from the file with more ~25bp reads
echo "Step 3: Determining which file contains the 25N variable regions..."

r1_25bp=$(zcat "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" | awk 'NR%4==2 {if(length($0) >= 20 && length($0) <= 30) count++} END {print count+0}')
r2_25bp=$(zcat "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" | awk 'NR%4==2 {if(length($0) >= 20 && length($0) <= 30) count++} END {print count+0}')

total_r1=$(zcat "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" | awk 'NR%4==2' | wc -l)
total_r2=$(zcat "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" | awk 'NR%4==2' | wc -l)

r1_percent=$(awk -v a=$r1_25bp -v b=$total_r1 'BEGIN{printf "%.1f", 100*a/b}')
r2_percent=$(awk -v a=$r2_25bp -v b=$total_r2 'BEGIN{printf "%.1f", 100*a/b}')

printf "R1: %'d/%'d reads in 20-30bp range (%.1f%%)\n" $r1_25bp $total_r1 $r1_percent
printf "R2: %'d/%'d reads in 20-30bp range (%.1f%%)\n" $r2_25bp $total_r2 $r2_percent

# Decide which file(s) to use for counting
if (( $(echo "$r1_percent > $r2_percent + 5" | bc -l) )); then
    echo "✅ Using R1 for sequence extraction (higher 25bp percentage: ${r1_percent}%)"
    SOURCE_FILE="${OUTPUT_PREFIX}_R1_trimmed.fastq.gz"
    SOURCE_NAME="R1"
elif (( $(echo "$r2_percent > $r1_percent + 5" | bc -l) )); then
    echo "✅ Using R2 for sequence extraction (higher 25bp percentage: ${r2_percent}%)"
    echo "⚠️  Note: R2 will be reverse-complemented to get biological sequence orientation"
    SOURCE_FILE="${OUTPUT_PREFIX}_R2_trimmed.fastq.gz"  
    SOURCE_NAME="R2-RC"
else
    echo "⚠️  Similar percentages in both files (R1:${r1_percent}%, R2:${r2_percent}%)"
    echo "✅ Using both R1 and R2 for sequence extraction"
    echo "📝 Note: R2 sequences will be reverse-complemented before merging to avoid duplicates"
    SOURCE_FILE="both"
    SOURCE_NAME="R1+R2-RC"
fi

# Step 4: Count sequences
echo "Step 4: Extracting and counting sequences..."

if [ "$SOURCE_FILE" = "both" ]; then
    # Combine both files - R2 needs reverse complement!
    echo "Converting R2 to reverse complement before merging..."
    {
        zcat "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" | awk 'NR%4==2'
        zcat "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz" | awk 'NR%4==2' | \
            rev | tr 'ATCG' 'TAGC'  # Reverse complement R2
    } | LC_ALL=C sort | uniq -c | LC_ALL=C sort -nr > "${OUTPUT_PREFIX}_raw_counts.txt"
elif [[ "$SOURCE_FILE" == *"R2"* ]]; then
    # R2 only - needs reverse complement
    echo "Converting R2 sequences to reverse complement..."
    zcat "$SOURCE_FILE" | awk 'NR%4==2' | rev | tr 'ATCG' 'TAGC' | \
        LC_ALL=C sort | uniq -c | LC_ALL=C sort -nr > "${OUTPUT_PREFIX}_raw_counts.txt"
else
    # R1 only - use as-is
    zcat "$SOURCE_FILE" | awk 'NR%4==2' | LC_ALL=C sort | uniq -c | LC_ALL=C sort -nr > "${OUTPUT_PREFIX}_raw_counts.txt"
fi

# Step 5: Create output files
echo "Step 5: Creating output files..."

# Filter sequences in reasonable length range (15-40bp to be inclusive)
awk '$1 >= 1 && length($2) >= 15 && length($2) <= 40' "${OUTPUT_PREFIX}_raw_counts.txt" > "${OUTPUT_PREFIX}_filtered_counts.txt"

# Create FASTA
awk '{print ">"NR"_count_"$1"_len_"length($2)"\n"$2}' "${OUTPUT_PREFIX}_filtered_counts.txt" > "${OUTPUT_PREFIX}_sequences.fasta"

# Create count table
total_reads=$(awk '{sum+=$1} END {print sum}' "${OUTPUT_PREFIX}_filtered_counts.txt")
awk -v total=$total_reads '{
    printf "%s\t%d\t%d\t%.6f\n", $2, $1, length($2), $1/total
}' "${OUTPUT_PREFIX}_filtered_counts.txt" > "${OUTPUT_PREFIX}_counts.txt"

# Add header
sed -i '1i sequence\tcount\tlength\tfrequency' "${OUTPUT_PREFIX}_counts.txt"

# Create summary
echo "Top 20 most abundant sequences:" > "${OUTPUT_PREFIX}_summary.txt"
echo "Rank\tSequence\tLength\tCount\tFrequency" >> "${OUTPUT_PREFIX}_summary.txt"
head -20 "${OUTPUT_PREFIX}_filtered_counts.txt" | awk -v total=$total_reads '{
    printf "%d\t%s\t%d\t%d\t%.4f\n", NR, $2, length($2), $1, $1/total
}' >> "${OUTPUT_PREFIX}_summary.txt"

# Create length distribution
echo "length\tnum_unique_sequences\ttotal_reads" > "${OUTPUT_PREFIX}_length_distribution.txt"
awk '{
    len=length($2); 
    count[len]+=$1; 
    unique[len]++
} END {
    for(l in count) printf "%d\t%d\t%d\n", l, unique[l], count[l]
}' "${OUTPUT_PREFIX}_filtered_counts.txt" | sort -n >> "${OUTPUT_PREFIX}_length_distribution.txt"

# Step 6: Generate statistics
echo "Step 6: Generating statistics..."

unique_sequences=$(wc -l < "${OUTPUT_PREFIX}_filtered_counts.txt")
sequences_25bp=$(awk 'length($2) >= 23 && length($2) <= 27' "${OUTPUT_PREFIX}_filtered_counts.txt" | wc -l)
reads_25bp=$(awk 'length($2) >= 23 && length($2) <= 27 {sum+=$1} END {print sum+0}' "${OUTPUT_PREFIX}_filtered_counts.txt")

echo ""
echo "=== EXTRACTION RESULTS ==="
printf "Source used for counting: %s\n" "$SOURCE_NAME"
printf "Total unique sequences (15-40bp): %'d\n" $unique_sequences
printf "Unique sequences in 23-27bp range: %'d\n" $sequences_25bp
printf "Total reads in 23-27bp range: %'d\n" $reads_25bp
printf "Total reads after filtering: %'d\n" $total_reads

if [ $reads_25bp -gt 0 ]; then
    echo "✅ Successfully extracted sequences in expected 25bp range!"
else
    echo "⚠️  Few sequences in 23-27bp range. Check length distribution file."
fi

# Cleanup intermediate files
rm -f "${OUTPUT_PREFIX}_R1_trimmed.fastq.gz" "${OUTPUT_PREFIX}_R2_trimmed.fastq.gz"
rm -f "${OUTPUT_PREFIX}_raw_counts.txt" "${OUTPUT_PREFIX}_filtered_counts.txt"

echo ""
echo "Files created:"
echo "  ${OUTPUT_PREFIX}_sequences.fasta - FASTA with abundance"
echo "  ${OUTPUT_PREFIX}_counts.txt - Count table with length info"
echo "  ${OUTPUT_PREFIX}_summary.txt - Top 20 sequences"
echo "  ${OUTPUT_PREFIX}_length_distribution.txt - Length distribution"
echo "  ${OUTPUT_PREFIX}_fastp.html - Quality control report"
echo "  ${OUTPUT_PREFIX}_fastp.json - Detailed statistics"

echo ""
echo "Next steps:"
echo "1. Check ${OUTPUT_PREFIX}_fastp.html for quality control"
echo "2. Review length distribution to confirm 25bp enrichment" 
echo "3. Proceed with clustering analysis if results look good"

echo ""
echo "Pipeline completed successfully!"
