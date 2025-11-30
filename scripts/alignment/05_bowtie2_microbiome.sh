#!/bin/bash
# =============================================================================
# Microbial Gene Quantification with Bowtie2 + featureCounts
# =============================================================================
# This script aligns reads to the HOMD (Human Oral Microbiome Database)
# reference genomes and quantifies gene expression using featureCounts
# with PROKKA annotations.
#
# The two-step approach:
#   1. Bowtie2: Align reads to HOMD reference genomes
#   2. featureCounts: Count reads per gene using PROKKA GFF annotations
#
# Reference: HOMD V10.1 with PROKKA gene annotations
#
# Input:  mRNA-enriched FASTQ files (from rRNA filtering)
# Output: Gene count matrix for differential expression analysis
#
# Author: Maria J. Rus
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - Modify these paths for your system
# -----------------------------------------------------------------------------
INPUT_DIR="preprocessing/04_rrna_filtered/fastq"
OUTPUT_DIR="alignment/microbiome_homd"
THREADS=8

# HOMD reference database
# Build Bowtie2 index with: bowtie2-build HOMD_genomes.fna HOMD_v10.1
HOMD_INDEX="/path/to/databases/HOMD_v10.1/HOMD_v10.1"

# PROKKA annotation file (GFF format)
# Contains gene coordinates and annotations for all HOMD genomes
PROKKA_GFF="/path/to/databases/HOMD_v10.1/PROKKA/ALL_genomes.gff"

# Bowtie2 parameters
# --very-sensitive: Maximum sensitivity for accurate alignment
# --no-unal: Do not output unaligned reads to SAM
BOWTIE2_PARAMS="--very-sensitive --no-unal"

# featureCounts parameters
# -p: Paired-end mode
# -B: Both reads must be aligned
# -C: Do not count chimeric fragments
# -t CDS: Count reads overlapping CDS features
# -g ID: Use 'ID' attribute as gene identifier
FEATURECOUNTS_PARAMS="-p -B -C -t CDS -g ID"

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Microbial Gene Quantification Pipeline"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  - HOMD index: ${HOMD_INDEX}"
echo "  - PROKKA annotation: ${PROKKA_GFF}"
echo "  - Threads: ${THREADS}"
echo ""

mkdir -p "${OUTPUT_DIR}/bam"
mkdir -p "${OUTPUT_DIR}/counts"
mkdir -p "${OUTPUT_DIR}/logs"

# -----------------------------------------------------------------------------
# Step 1: Align reads to HOMD with Bowtie2
# -----------------------------------------------------------------------------
echo "Step 1: Aligning reads to HOMD reference..."
echo ""

# Initialize alignment summary
ALIGN_SUMMARY="${OUTPUT_DIR}/bowtie2_alignment_summary.tsv"
echo -e "Sample\tTotal_Reads\tAligned\tAlignment_Rate" > "${ALIGN_SUMMARY}"

# Process each filtered sample
for R1_FILE in "${INPUT_DIR}"/*_mRNA_R1.fastq.gz "${INPUT_DIR}"/*_R1*.fastq.gz 2>/dev/null; do
    # Skip if no files found
    [[ -e "${R1_FILE}" ]] || continue

    # Determine R2 file and sample name
    SAMPLE_NAME=$(basename "${R1_FILE}" | sed 's/_mRNA_R1.*//g; s/_R1.*//g')
    R2_FILE="${R1_FILE/_R1/_R2}"

    # Check if R2 exists
    if [[ ! -f "${R2_FILE}" ]]; then
        echo "Warning: R2 file not found for ${SAMPLE_NAME}, skipping..."
        continue
    fi

    echo "  Aligning: ${SAMPLE_NAME}"

    # Output files
    BAM_FILE="${OUTPUT_DIR}/bam/${SAMPLE_NAME}_homd.bam"
    LOG_FILE="${OUTPUT_DIR}/logs/${SAMPLE_NAME}_bowtie2.log"

    # Run Bowtie2
    bowtie2 \
        -x "${HOMD_INDEX}" \
        -1 "${R1_FILE}" \
        -2 "${R2_FILE}" \
        ${BOWTIE2_PARAMS} \
        --threads "${THREADS}" \
        2> "${LOG_FILE}" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" -o "${BAM_FILE}" -

    # Index BAM
    samtools index "${BAM_FILE}"

    # Extract alignment statistics
    TOTAL_READS=$(grep "reads; of these:" "${LOG_FILE}" | awk '{print $1}')
    ALIGNMENT_RATE=$(grep "overall alignment rate" "${LOG_FILE}" | awk '{print $1}')

    echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${ALIGNMENT_RATE}" >> "${ALIGN_SUMMARY}"
    echo "    Alignment rate: ${ALIGNMENT_RATE}"
done

echo ""
echo "Bowtie2 alignment completed."
echo ""

# -----------------------------------------------------------------------------
# Step 2: Count reads per gene with featureCounts
# -----------------------------------------------------------------------------
echo "Step 2: Counting reads per gene with featureCounts..."
echo ""

# Collect all BAM files
BAM_FILES=$(ls "${OUTPUT_DIR}/bam/"*_homd.bam | tr '\n' ' ')

# Output file
COUNT_FILE="${OUTPUT_DIR}/counts/all_samples_prokka_counts.txt"
FEATURECOUNTS_LOG="${OUTPUT_DIR}/logs/featureCounts.log"

# Run featureCounts
featureCounts \
    ${FEATURECOUNTS_PARAMS} \
    -T "${THREADS}" \
    -a "${PROKKA_GFF}" \
    -o "${COUNT_FILE}" \
    ${BAM_FILES} \
    2>&1 | tee "${FEATURECOUNTS_LOG}"

# featureCounts outputs a summary file automatically
mv "${COUNT_FILE}.summary" "${OUTPUT_DIR}/counts/featureCounts_summary.txt"

# -----------------------------------------------------------------------------
# Step 3: Clean up count matrix
# -----------------------------------------------------------------------------
echo ""
echo "Step 3: Cleaning count matrix..."

# featureCounts output has metadata columns (Chr, Start, End, Strand, Length)
# Remove these for downstream analysis
CLEAN_COUNT_FILE="${OUTPUT_DIR}/counts/gene_counts_matrix.tsv"

# Extract header (sample names) - column 1 is Geneid, columns 2-6 are metadata, rest are samples
head -2 "${COUNT_FILE}" | tail -1 | cut -f1,7- | sed 's|'${OUTPUT_DIR}'/bam/||g; s|_homd.bam||g' > "${CLEAN_COUNT_FILE}"

# Extract counts (skip header lines and metadata columns)
tail -n +3 "${COUNT_FILE}" | cut -f1,7- >> "${CLEAN_COUNT_FILE}"

echo "Cleaned count matrix saved to: ${CLEAN_COUNT_FILE}"

# -----------------------------------------------------------------------------
# Step 4: Generate count statistics
# -----------------------------------------------------------------------------
echo ""
echo "Step 4: Generating statistics..."

# Create per-sample statistics
STATS_FILE="${OUTPUT_DIR}/counts/count_statistics.tsv"
echo -e "Sample\tTotal_Assigned\tTotal_Genes_Detected" > "${STATS_FILE}"

# Read featureCounts summary for assigned reads
SUMMARY_FILE="${OUTPUT_DIR}/counts/featureCounts_summary.txt"

# Parse summary to get assigned counts per sample
# Format: Status | Sample1 | Sample2 | ...
awk 'NR==1 {
    for(i=2; i<=NF; i++) {
        gsub(/.*\//, "", $i);
        gsub(/_homd.bam/, "", $i);
        samples[i] = $i
    }
}
NR>1 && $1=="Assigned" {
    for(i=2; i<=NF; i++) assigned[i] = $i
}
END {
    for(i=2; i<=length(samples)+1; i++) {
        print samples[i] "\t" assigned[i]
    }
}' "${SUMMARY_FILE}" >> "${STATS_FILE}"

# Count genes detected per sample (genes with >0 counts)
echo "Counting detected genes per sample..."

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Microbial Gene Quantification Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - BAM files: ${OUTPUT_DIR}/bam/"
echo "  - Raw count table: ${COUNT_FILE}"
echo "  - Clean count matrix: ${CLEAN_COUNT_FILE}"
echo "  - featureCounts summary: ${OUTPUT_DIR}/counts/featureCounts_summary.txt"
echo "  - Alignment summary: ${ALIGN_SUMMARY}"
echo "  - Count statistics: ${STATS_FILE}"
echo ""
echo "Alignment summary:"
cat "${ALIGN_SUMMARY}"
echo ""
echo "The clean count matrix (${CLEAN_COUNT_FILE}) is ready for:"
echo "  - Differential expression analysis with DESeq2"
echo "  - Functional enrichment with eggNOG-mapper annotations"
echo ""
echo "Next step: Run functional annotation (functional_annotation/02_eggnog_mapper.sh)"
echo ""

# -----------------------------------------------------------------------------
# Notes
# -----------------------------------------------------------------------------
# Expected alignment rates for oral metatranscriptomics:
# - Alignment to HOMD: 20-60% depending on sample type and quality
# - Lower rates may indicate non-oral bacteria or sequencing issues
# - Very high rates (>80%) suggest good sample quality and database coverage
# -----------------------------------------------------------------------------
