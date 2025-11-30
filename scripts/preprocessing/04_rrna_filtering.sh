#!/bin/bash
# =============================================================================
# Ribosomal RNA Filtering with SortMeRNA
# =============================================================================
# This script removes ribosomal RNA (rRNA) sequences from metatranscriptomic
# samples using SortMeRNA. rRNA can constitute >90% of total RNA in bacterial
# samples and provides no useful information for gene expression analysis.
#
# Databases used:
#   - SILVA 138.1: 16S/18S and 23S/28S rRNA
#   - RFAM 14.9: 5S and 5.8S rRNA
#
# Input:  Non-host FASTQ files (from host removal step)
# Output: mRNA-enriched FASTQ files (rRNA-depleted)
#
# Author: Maria J. Rus
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configuration - Modify these paths for your system
# -----------------------------------------------------------------------------
INPUT_DIR="preprocessing/03_host_removed/fastq"
OUTPUT_DIR="preprocessing/04_rrna_filtered"
THREADS=8

# SortMeRNA database directory
# Download from: https://github.com/biocore/sortmerna/releases
SORTMERNA_DB="/path/to/databases/sortmerna"

# rRNA reference databases (SILVA + RFAM)
# These should be indexed before first use with: sortmerna --index
RRNA_DBS=(
    "${SORTMERNA_DB}/smr_v4.3_default_db/silva-bac-16s-id90.fasta"
    "${SORTMERNA_DB}/smr_v4.3_default_db/silva-bac-23s-id98.fasta"
    "${SORTMERNA_DB}/smr_v4.3_default_db/silva-euk-18s-id95.fasta"
    "${SORTMERNA_DB}/smr_v4.3_default_db/silva-euk-28s-id98.fasta"
    "${SORTMERNA_DB}/smr_v4.3_default_db/rfam-5s-database-id98.fasta"
    "${SORTMERNA_DB}/smr_v4.3_default_db/rfam-5.8s-database-id98.fasta"
)

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------
echo "=============================================="
echo "rRNA Filtering Pipeline - SortMeRNA"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  - Database directory: ${SORTMERNA_DB}"
echo "  - Threads: ${THREADS}"
echo ""

mkdir -p "${OUTPUT_DIR}/fastq"
mkdir -p "${OUTPUT_DIR}/rrna"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/workdir"

# Build --ref arguments for SortMeRNA
REF_ARGS=""
for db in "${RRNA_DBS[@]}"; do
    if [[ -f "${db}" ]]; then
        REF_ARGS="${REF_ARGS} --ref ${db}"
    else
        echo "Warning: Database not found: ${db}"
    fi
done

if [[ -z "${REF_ARGS}" ]]; then
    echo "Error: No rRNA databases found. Please check SORTMERNA_DB path."
    exit 1
fi

# -----------------------------------------------------------------------------
# Process each sample
# -----------------------------------------------------------------------------
echo "Processing samples..."
echo ""

# Initialize summary file
SUMMARY_FILE="${OUTPUT_DIR}/sortmerna_summary.tsv"
echo -e "Sample\tTotal_Reads\trRNA_Reads\tmRNA_Reads\trRNA_Percent\tmRNA_Percent" > "${SUMMARY_FILE}"

# Process each non-host sample
for R1_FILE in "${INPUT_DIR}"/*_nonhost_R1.fastq.gz "${INPUT_DIR}"/*_R1*.fastq.gz 2>/dev/null; do
    # Skip if no files found
    [[ -e "${R1_FILE}" ]] || continue

    # Determine R2 file and sample name
    SAMPLE_NAME=$(basename "${R1_FILE}" | sed 's/_nonhost_R1.*//g; s/_R1.*//g')
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
    OUT_PREFIX="${OUTPUT_DIR}/fastq/${SAMPLE_NAME}"
    RRNA_PREFIX="${OUTPUT_DIR}/rrna/${SAMPLE_NAME}"
    WORKDIR="${OUTPUT_DIR}/workdir/${SAMPLE_NAME}"
    LOG_FILE="${OUTPUT_DIR}/logs/${SAMPLE_NAME}_sortmerna.log"

    # Create sample-specific workdir
    mkdir -p "${WORKDIR}"

    # Run SortMeRNA
    # --paired_in: Input is paired-end
    # --paired_out: Keep paired-end structure in output
    # --out2: Output to two files (R1 and R2)
    # --other: Non-rRNA (mRNA) reads
    # --aligned: rRNA reads (for QC)
    # --fastx: Output in FASTQ format
    sortmerna \
        ${REF_ARGS} \
        --reads "${R1_FILE}" \
        --reads "${R2_FILE}" \
        --workdir "${WORKDIR}" \
        --threads "${THREADS}" \
        --paired_in \
        --paired_out \
        --out2 \
        --other "${OUT_PREFIX}" \
        --aligned "${RRNA_PREFIX}" \
        --fastx \
        2>&1 | tee "${LOG_FILE}"

    # SortMeRNA outputs files with specific naming
    # Rename to standard convention
    if [[ -f "${OUT_PREFIX}_fwd.fq.gz" ]]; then
        mv "${OUT_PREFIX}_fwd.fq.gz" "${OUT_PREFIX}_mRNA_R1.fastq.gz"
        mv "${OUT_PREFIX}_rev.fq.gz" "${OUT_PREFIX}_mRNA_R2.fastq.gz"
    elif [[ -f "${OUT_PREFIX}_fwd.fq" ]]; then
        gzip -c "${OUT_PREFIX}_fwd.fq" > "${OUT_PREFIX}_mRNA_R1.fastq.gz"
        gzip -c "${OUT_PREFIX}_rev.fq" > "${OUT_PREFIX}_mRNA_R2.fastq.gz"
        rm -f "${OUT_PREFIX}_fwd.fq" "${OUT_PREFIX}_rev.fq"
    fi

    # Compress rRNA files if not already
    if [[ -f "${RRNA_PREFIX}_fwd.fq" ]]; then
        gzip -c "${RRNA_PREFIX}_fwd.fq" > "${RRNA_PREFIX}_rRNA_R1.fastq.gz"
        gzip -c "${RRNA_PREFIX}_rev.fq" > "${RRNA_PREFIX}_rRNA_R2.fastq.gz"
        rm -f "${RRNA_PREFIX}_fwd.fq" "${RRNA_PREFIX}_rev.fq"
    fi

    # Extract statistics from SortMeRNA log
    TOTAL_READS=$(grep "Total reads" "${LOG_FILE}" | awk -F'=' '{print $2}' | tr -d ' ' || echo "NA")
    RRNA_READS=$(grep "rRNA reads" "${LOG_FILE}" | awk -F'=' '{print $2}' | tr -d ' ' || echo "NA")
    MRNA_READS=$(grep "other reads" "${LOG_FILE}" | awk -F'=' '{print $2}' | tr -d ' ' || echo "NA")

    # Calculate percentages
    if [[ "${TOTAL_READS}" != "NA" && "${TOTAL_READS}" != "0" ]]; then
        RRNA_PCT=$(echo "scale=2; ${RRNA_READS}/${TOTAL_READS}*100" | bc 2>/dev/null || echo "NA")
        MRNA_PCT=$(echo "scale=2; ${MRNA_READS}/${TOTAL_READS}*100" | bc 2>/dev/null || echo "NA")
    else
        RRNA_PCT="NA"
        MRNA_PCT="NA"
    fi

    echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${RRNA_READS}\t${MRNA_READS}\t${RRNA_PCT}%\t${MRNA_PCT}%" >> "${SUMMARY_FILE}"

    # Clean up workdir to save space
    rm -rf "${WORKDIR}"

    echo "  rRNA content: ${RRNA_PCT}%"
    echo "  mRNA retained: ${MRNA_PCT}%"
    echo "  Completed: ${SAMPLE_NAME}"
    echo ""
done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "rRNA Filtering Complete"
echo "=============================================="
echo ""
echo "Output files:"
echo "  - mRNA-enriched FASTQ files: ${OUTPUT_DIR}/fastq/"
echo "  - rRNA reads (for QC): ${OUTPUT_DIR}/rrna/"
echo "  - SortMeRNA logs: ${OUTPUT_DIR}/logs/"
echo "  - Summary statistics: ${SUMMARY_FILE}"
echo ""
echo "Summary statistics:"
cat "${SUMMARY_FILE}"
echo ""
echo "Next step: Proceed to alignment (alignment/01_star_human_alignment.sh)"
echo ""

# -----------------------------------------------------------------------------
# Notes on expected results
# -----------------------------------------------------------------------------
# For metatranscriptomic samples:
# - rRNA content typically ranges from 40-90% depending on sample type
# - Higher mRNA retention indicates better library prep or rRNA depletion
# - Some rRNA is always present and normal
# - Very low rRNA (<20%) might indicate degraded samples
# -----------------------------------------------------------------------------
