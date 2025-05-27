#!/bin/bash
# All default values, can be overridden based on your requirement. 

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference_vcf> <region> <chunk_size>"
    exit 1
fi

echo "Started chunking regions"

pathh="glimpse_static_bins"
PROC=32 #<threads, adjust it as per your needs> 

GLIMPSE_CHUNK=${pathh}/GLIMPSE_chunk_static
REF_VCF=$1
REGION=$2
CHUNK_SIZE=$3

${GLIMPSE_CHUNK} --thread ${PROC} --input ${REF_VCF} --region ${REGION} --window-size 2000000 --buffer-size 200000 --output chunks.${REGION}.c${CHUNK_SIZE}.txt

echo "Done chunking regions"

