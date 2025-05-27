#!/bin/bash

# Usage: run_glimpse_impute.sh <input_vcf> <ref_panel> <map_file> <chromosome>
# Example: run_glimpse_impute.sh input.vcf.gz ref_panel.vcf.gz map_file.map chr1

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_vcf> <ref_panel> <map_file> <chromosome>"
    exit 1
fi

input_vcf=$1
ref_panel=$2
map_file=$3
chromosome=$4

# Paths can be set here or passed as environment variables
ipath="/path/to/input_vcfs"  # Replace with actual path or use environment variable
refpath="/path/to/ref_panel_files"  # Replace with actual path or use environment variable
mappath="/path/to/glimpse_map_files"  # Replace with actual path or use environment variable

# Run Glimpse Chunk
echo "Starting glimpse_chunk for ${chromosome}"
./glimpse_chunk.sh ${refpath}/${ref_panel} ${chromosome}

# Run Glimpse Phase
echo "Starting glimpse_phase for ${chromosome}"
./glimpse_phase.sh ${ipath}/${input_vcf} ${refpath}/${ref_panel} ${mappath}/${map_file} GLIMPSE_impute_${chromosome} chunks.${chromosome}.txt ${chromosome}_imputing

# Run Glimpse Ligate
echo "Starting glimpse_ligate for ${chromosome}"
./glimpse_ligate.sh GLIMPSE_ligate_${chromosome} GLIMPSE_impute_${chromosome} ${chromosome}_imputing 

# Run Glimpse Sample
echo "Starting glimpse_sample for ${chromosome}"
./glimpse_sample.sh GLIMPSE_phase_${chromosome} GLIMPSE_ligate_${chromosome} ${chromosome}_imputing 

echo "Glimpse imputation process complete for chromosome ${chromosome}."

