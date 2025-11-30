#!/bin/bash
# =============================================================================
# Quality Control with FastQC and MultiQC
# =============================================================================
# This script performs initial quality assessment of raw RNA-seq reads
# using FastQC for individual sample reports and MultiQC for aggregated
# visualization.
#
# Input:  Raw FASTQ files (paired-end)
# Output: FastQC reports and MultiQC summary
#
# Author: Maria J. Rus
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - Modify these paths for your system
# -----------------------------------------------------------------------------
RAW_DATA_DIR="raw_data"
OUTPUT_DIR="preprocessing/01_fastqc"
THREADS=8

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
echo "=============================================="
echo "Quality Control Pipeline - FastQC + MultiQC"
echo "=============================================="
echo ""
echo "Creating output directories..."

mkdir -p "${OUTPUT_DIR}/raw"
mkdir -p "${OUTPUT_DIR}/multiqc"

# -----------------------------------------------------------------------------
# Run FastQC on all raw FASTQ files
# -----------------------------------------------------------------------------
echo ""
echo "Running FastQC on raw data..."
echo "Input directory: ${RAW_DATA_DIR}"
echo "Output directory: ${OUTPUT_DIR}/raw"
echo ""

# Find all FASTQ files and run FastQC
fastqc \
    --outdir "${OUTPUT_DIR}/raw" \
    --threads "${THREADS}" \
    --quiet \
    "${RAW_DATA_DIR}"/*.fastq.gz

echo "FastQC completed for all samples."

# -----------------------------------------------------------------------------
# Aggregate reports with MultiQC
# -----------------------------------------------------------------------------
echo ""
echo "Generating MultiQC report..."

multiqc \
    --outdir "${OUTPUT_DIR}/multiqc" \
    --filename "multiqc_raw_reads" \
    --title "Raw Reads Quality Report - Oral Metatranscriptomics" \
    --force \
    "${OUTPUT_DIR}/raw"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Quality Control Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - Individual FastQC reports: ${OUTPUT_DIR}/raw/"
echo "  - Aggregated MultiQC report: ${OUTPUT_DIR}/multiqc/multiqc_raw_reads.html"
echo ""
echo "Next step: Review quality metrics and proceed to trimming (02_fastp_trimming.sh)"
echo ""
