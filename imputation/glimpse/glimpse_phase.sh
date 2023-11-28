#!/bin/bash

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_vcf> <reference_panel> <map_file> <output_dir> <chunk_file> <output_file_tag>"
    exit 1
fi

VCF=$1
REF=$2
MAP=$3
DIR_NAME=$4
chunk_file=$5
outfile_tag=$6
pathh="glimpse_static_bins"
PROC=32  # Default value, can be overridden by passing as an argument

mkdir -p ${DIR_NAME}

echo "Started processing ..."

while IFS="" read -r LINE || [ -n "$LINE" ]; do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    OUT=${DIR_NAME}/${outfile_tag}.${ID}.bcf

    echo "Processing chunk ${ID}..."
    ${pathh}/GLIMPSE_phase_static --thread ${PROC} --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}  
    bcftools index -f ${OUT}
done < ${chunk_file}

echo "Processing complete."

