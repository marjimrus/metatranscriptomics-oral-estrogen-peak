#!/bin/bash
# =============================================================================
# Host (Human) Removal with Bowtie2
# =============================================================================
# This script removes human-derived reads from metatranscriptomic samples
# using Bowtie2 alignment against the GRCh38 human reference genome.
#
# For oral samples, human RNA can constitute a significant proportion of reads.
# Removing host reads before downstream analysis:
#   - Reduces computational burden
#   - Prevents false positive classifications
#   - Increases signal-to-noise ratio for microbial analysis
#
# Input:  Trimmed FASTQ files (from fastp)
# Output: Non-host (microbial) FASTQ files
#
# Author: Maria J. Rus
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - Modify these paths for your system
# -----------------------------------------------------------------------------
INPUT_DIR="preprocessing/02_trimmed/fastq"
OUTPUT_DIR="preprocessing/03_host_removed"
THREADS=8

# Human genome reference (Bowtie2 index prefix)
# Build with: bowtie2-build GRCh38.fna GRCh38
HUMAN_INDEX="/path/to/databases/GRCh38/GRCh38"

# Bowtie2 parameters
# --very-sensitive: Maximum sensitivity for accurate host detection
# --un-conc-gz: Output non-concordantly aligned (non-host) reads compressed
BOWTIE2_PARAMS="--very-sensitive"

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Host Removal Pipeline - Bowtie2"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  - Human reference: ${HUMAN_INDEX}"
echo "  - Bowtie2 mode: ${BOWTIE2_PARAMS}"
echo "  - Threads: ${THREADS}"
echo ""

mkdir -p "${OUTPUT_DIR}/fastq"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/stats"

# -----------------------------------------------------------------------------
# Process each sample
# -----------------------------------------------------------------------------
echo "Processing samples..."
echo ""

# Initialize summary file
SUMMARY_FILE="${OUTPUT_DIR}/host_removal_summary.tsv"
echo -e "Sample\tTotal_Reads\tAligned_Host\tUnaligned_NonHost\tHost_Percent\tRetention_Percent" > "${SUMMARY_FILE}"

# Process each trimmed sample
for R1_FILE in "${INPUT_DIR}"/*_trimmed_R1.fastq.gz "${INPUT_DIR}"/*_R1*.fastq.gz 2>/dev/null; do
    # Skip if no files found
    [[ -e "${R1_FILE}" ]] || continue

    # Determine R2 file and sample name
    SAMPLE_NAME=$(basename "${R1_FILE}" | sed 's/_trimmed_R1.*//g; s/_R1.*//g')
    R2_FILE="${R1_FILE/_R1/_R2}"

    # Check if R2 exists
    if [[ ! -f "${R2_FILE}" ]]; then
        echo "Warning: R2 file not found for ${SAMPLE_NAME}, skipping..."
        continue
    fi

    echo "Processing: ${SAMPLE_NAME}"
    echo "  Input R1: ${R1_FILE}"
    echo "  Input R2: ${R2_FILE}"

    # Output files
    OUT_R1="${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R2.fastq.gz"
    LOG_FILE="${OUTPUT_DIR}/logs/${SAMPLE_NAME}_bowtie2.log"
    SAM_FILE="${OUTPUT_DIR}/stats/${SAMPLE_NAME}_aligned.sam"

    # Run Bowtie2
    # --un-conc-gz outputs reads that do NOT align concordantly (non-host reads)
    bowtie2 \
        -x "${HUMAN_INDEX}" \
        -1 "${R1_FILE}" \
        -2 "${R2_FILE}" \
        ${BOWTIE2_PARAMS} \
        --threads "${THREADS}" \
        --un-conc-gz "${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R%.fastq.gz" \
        -S "${SAM_FILE}" \
        2>&1 | tee "${LOG_FILE}"

    # Rename output files to standard naming convention
    if [[ -f "${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R1.fastq.gz" ]]; then
        mv "${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R1.fastq.gz" "${OUT_R1}" 2>/dev/null || true
        mv "${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_nonhost_R2.fastq.gz" "${OUT_R2}" 2>/dev/null || true
    fi

    # Extract statistics from Bowtie2 log
    TOTAL_READS=$(grep "reads; of these:" "${LOG_FILE}" | awk '{print $1}')
    ALIGNED=$(grep "aligned concordantly" "${LOG_FILE}" | head -1 | awk '{print $1}')
    if [[ -z "${ALIGNED}" ]]; then
        ALIGNED=$(grep "overall alignment rate" "${LOG_FILE}" | awk '{print $1}' | sed 's/%//')
    fi

    # Calculate retention (non-host reads)
    OVERALL_RATE=$(grep "overall alignment rate" "${LOG_FILE}" | awk '{print $1}' | sed 's/%//')
    RETENTION=$(echo "scale=2; 100 - ${OVERALL_RATE}" | bc)

    echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${OVERALL_RATE}%\t${RETENTION}%\t${OVERALL_RATE}\t${RETENTION}" >> "${SUMMARY_FILE}"

    # Remove SAM file to save space (alignment not needed downstream)
    rm -f "${SAM_FILE}"

    echo "  Host alignment rate: ${OVERALL_RATE}%"
    echo "  Retained (non-host): ${RETENTION}%"
    echo "  Completed: ${SAMPLE_NAME}"
    echo ""
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Host Removal Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - Non-host FASTQ files: ${OUTPUT_DIR}/fastq/"
echo "  - Bowtie2 logs: ${OUTPUT_DIR}/logs/"
echo "  - Summary statistics: ${SUMMARY_FILE}"
echo ""
echo "Summary statistics:"
cat "${SUMMARY_FILE}"
echo ""
echo "Next step: Run rRNA filtering (04_rrna_filtering.sh)"
echo ""

# -----------------------------------------------------------------------------
# Notes on expected results
# -----------------------------------------------------------------------------
# For oral metatranscriptomic samples:
# - Host content can vary from 50-90% depending on sample type
# - Subgingival plaque samples typically have lower host content than saliva
# - Higher bacterial content (lower host %) indicates successful sampling
# -----------------------------------------------------------------------------
