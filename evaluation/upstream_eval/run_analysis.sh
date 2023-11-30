#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: ./run_analysis.sh <benchmark_vcf> <input_vcf> <score_threshold> <output_prefix>"
    exit 1
fi

set -e -o pipefail
benchmark_vcf=$1
input_vcf=$2
score_threshold=$3  # Score threshold, 0 for call accuracy, 200 for window accuracy
output_prefix=$4

# Ensure VCF files are block compressed and indexed
bcftools index ${benchmark_vcf}
bcftools index ${input_vcf}

# Merging and normalization 
bcftools merge --force-samples -Ou ${benchmark_vcf} ${input_vcf} | \
bcftools norm -m -any -Ou - | \

# score_matched_vcf is the binary compiled from score_matched_vcf.cpp
# <adjust the path for the binary o nbuilding>
#BPATH=<path to the binary>
${BPATH}/score_matched_vcf - ${score_threshold} ${output_prefix}
