#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <phase_dir> <ligate_dir> <output_file>"
    exit 1
fi

PHASE_DIR=$1
LIG_DIR=$2
outputfile=$3
pathh="glimpse_static_bins"
PROC=32  # Default value, can be overridden by passing as an argument

mkdir -p ${PHASE_DIR}
VCF=${LIG_DIR}/${outputfile}.merged.bcf
OUT=${PHASE_DIR}/${outputfile}.phased.bcf

echo "Starting the phasing step..."
${pathh}/GLIMPSE_sample_static --thread ${PROC} --input ${VCF} --solve --output ${OUT}
echo "Phasing step complete."

