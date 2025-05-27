#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: ./gen_upstream_metrics.sh <benchmark_vcf> <input_vcf> <window_size> <debugfile> <output_prefix>"
    exit 1
fi

# Input arguments
BENCHMARK_VCF=$1
INPUT_VCF=$2
SCORE_THRESHOLD=$3
DEBUGFILE=$4
OUTPUT_PREFIX=$5

# Enable strict error handling
set -e -o pipefail

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to validate if a file exists
validate_file() {
    local file=$1
    if [ ! -f "${file}" ]; then
        echo "Error: File '${file}' not found."
        exit 1
    fi
}

# Function to validate VCF file format
validate_vcf() {
    local file=$1
    if ! bcftools view -H "${file}" &>/dev/null; then
        echo "Error: '${file}' is not a valid VCF file."
        exit 1
    fi
}

# Function to bgzip and index a VCF file if needed
prepare_vcf() {
    local vcf_file=$1

    # Check if the file is already compressed
    if [[ ${vcf_file} != *.gz ]]; then
        log "Compressing ${vcf_file} with bgzip..."
        bgzip -c "${vcf_file}" > "${vcf_file}.gz"
        vcf_file="${vcf_file}.gz"
    fi

    # Index the VCF file
    log "Indexing ${vcf_file}..."
    if ! bcftools index -f "${vcf_file}"; then
        echo "Error indexing ${vcf_file}"
        exit 1
    fi

    # Return the processed file path
    echo "${vcf_file}"
}

# Validate input files
log "Validating input files..."
validate_file "${BENCHMARK_VCF}"
validate_file "${INPUT_VCF}"

# Prepare benchmark and input VCF files
log "Preparing benchmark VCF..."
BENCHMARK_VCF=$(prepare_vcf "${BENCHMARK_VCF}")
log "Preparing input VCF..."
INPUT_VCF=$(prepare_vcf "${INPUT_VCF}")

# Validate that the files are still valid after preparation
validate_vcf "${BENCHMARK_VCF}"
validate_vcf "${INPUT_VCF}"

# Merge, normalize, and calculate metrics
log "Starting bcftools merge and score computation..."
if ! bcftools merge --force-samples -Ou "${BENCHMARK_VCF}" "${INPUT_VCF}" | \
   bcftools norm -m -any -Ou - | \
   ./build/score_matched_vcf - "${SCORE_THRESHOLD}" "${DEBUGFILE}" > "${OUTPUT_PREFIX}.json"; then
    echo "Error during bcftools merge or score computation."
    exit 1
fi

log "Output stats written to ${OUTPUT_PREFIX}.json"

