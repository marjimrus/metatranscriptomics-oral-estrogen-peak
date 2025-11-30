#!/bin/bash
# =============================================================================
# Quality Trimming with fastp
# =============================================================================
# This script performs quality trimming and adapter removal using fastp.
# Parameters optimized for metatranscriptomic data:
#   - Quality threshold: Q20 (Phred score >= 20)
#   - Minimum read length: 70 bp
#   - Cut front: 10 bp (remove low-quality bases at read start)
#   - Adapter detection: Automatic
#
# Input:  Raw FASTQ files (paired-end)
# Output: Trimmed FASTQ files and QC reports
#
# Author: Maria J. Rus
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - Modify these paths for your system
# -----------------------------------------------------------------------------
RAW_DATA_DIR="raw_data"
OUTPUT_DIR="preprocessing/02_trimmed"
THREADS=8

# Trimming parameters
QUALITY_THRESHOLD=20    # Minimum Phred quality score
MIN_LENGTH=70           # Minimum read length after trimming
CUT_FRONT=10            # Remove first N bases from each read

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Quality Trimming Pipeline - fastp"
echo "=============================================="
echo ""
echo "Parameters:"
echo "  - Quality threshold: Q${QUALITY_THRESHOLD}"
echo "  - Minimum length: ${MIN_LENGTH} bp"
echo "  - Cut front: ${CUT_FRONT} bp"
echo "  - Adapter detection: Automatic"
echo ""
echo "Creating output directories..."

mkdir -p "${OUTPUT_DIR}/fastq"
mkdir -p "${OUTPUT_DIR}/reports"
mkdir -p "${OUTPUT_DIR}/json"

# -----------------------------------------------------------------------------
# Process each sample
# -----------------------------------------------------------------------------
echo ""
echo "Processing samples..."
echo ""

# Initialize summary file
SUMMARY_FILE="${OUTPUT_DIR}/fastp_summary.tsv"
echo -e "Sample\tRaw_Reads\tTrimmed_Reads\tRetention_Rate\tQ30_Before\tQ30_After" > "${SUMMARY_FILE}"

# Find unique sample prefixes (assuming naming pattern: SAMPLE_R1.fastq.gz, SAMPLE_R2.fastq.gz)
for R1_FILE in "${RAW_DATA_DIR}"/*_R1*.fastq.gz "${RAW_DATA_DIR}"/*_1.fastq.gz 2>/dev/null; do
    # Skip if no files found
    [[ -e "${R1_FILE}" ]] || continue

    # Determine R2 file name
    if [[ "${R1_FILE}" == *"_R1"* ]]; then
        R2_FILE="${R1_FILE/_R1/_R2}"
        SAMPLE_NAME=$(basename "${R1_FILE}" | sed 's/_R1.*//g')
    else
        R2_FILE="${R1_FILE/_1.fastq/_2.fastq}"
        SAMPLE_NAME=$(basename "${R1_FILE}" | sed 's/_1\.fastq.*//g')
    fi

    # Check if R2 exists
    if [[ ! -f "${R2_FILE}" ]]; then
        echo "Warning: R2 file not found for ${SAMPLE_NAME}, skipping..."
        continue
    fi

    echo "Processing: ${SAMPLE_NAME}"

    # Output file names
    OUT_R1="${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/fastq/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
    HTML_REPORT="${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.html"
    JSON_REPORT="${OUTPUT_DIR}/json/${SAMPLE_NAME}_fastp.json"

    # Run fastp
    fastp \
        --in1 "${R1_FILE}" \
        --in2 "${R2_FILE}" \
        --out1 "${OUT_R1}" \
        --out2 "${OUT_R2}" \
        --html "${HTML_REPORT}" \
        --json "${JSON_REPORT}" \
        --thread "${THREADS}" \
        --qualified_quality_phred "${QUALITY_THRESHOLD}" \
        --length_required "${MIN_LENGTH}" \
        --cut_front \
        --cut_front_window_size 4 \
        --cut_front_mean_quality "${QUALITY_THRESHOLD}" \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality "${QUALITY_THRESHOLD}" \
        --detect_adapter_for_pe \
        --correction \
        --overrepresentation_analysis \
        --report_title "${SAMPLE_NAME} - fastp Report" \
        2>&1 | tee "${OUTPUT_DIR}/reports/${SAMPLE_NAME}_fastp.log"

    # Extract statistics from JSON for summary
    if command -v jq &> /dev/null; then
        RAW_READS=$(jq '.summary.before_filtering.total_reads' "${JSON_REPORT}")
        TRIMMED_READS=$(jq '.summary.after_filtering.total_reads' "${JSON_REPORT}")
        Q30_BEFORE=$(jq '.summary.before_filtering.q30_rate' "${JSON_REPORT}")
        Q30_AFTER=$(jq '.summary.after_filtering.q30_rate' "${JSON_REPORT}")
        RETENTION=$(echo "scale=4; ${TRIMMED_READS}/${RAW_READS}*100" | bc)
        echo -e "${SAMPLE_NAME}\t${RAW_READS}\t${TRIMMED_READS}\t${RETENTION}%\t${Q30_BEFORE}\t${Q30_AFTER}" >> "${SUMMARY_FILE}"
    fi

    echo "  Completed: ${SAMPLE_NAME}"
    echo ""
done

# -----------------------------------------------------------------------------
# Generate MultiQC report for trimmed data
# -----------------------------------------------------------------------------
echo "Generating MultiQC report for fastp results..."

multiqc \
    --outdir "${OUTPUT_DIR}" \
    --filename "multiqc_fastp" \
    --title "fastp Trimming Report - Oral Metatranscriptomics" \
    --force \
    "${OUTPUT_DIR}/json"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Quality Trimming Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - Trimmed FASTQ files: ${OUTPUT_DIR}/fastq/"
echo "  - Individual reports: ${OUTPUT_DIR}/reports/"
echo "  - Summary statistics: ${SUMMARY_FILE}"
echo "  - MultiQC report: ${OUTPUT_DIR}/multiqc_fastp.html"
echo ""
echo "Next step: Run host removal (03_host_removal.sh)"
echo ""
