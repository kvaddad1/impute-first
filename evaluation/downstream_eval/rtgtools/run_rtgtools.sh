#!/bin/bash
# Usage: ./run_rtgtools.sh <truth_vcf> <input_vcf_directory> <reference_sdf> <bed_file_directory> <output_directory>

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <truth_vcf> <input_vcf_directory> <reference_sdf> <bed_file_directory> <output_directory>"
    exit 1
fi

TRUTH_VCF=$1
INPUT_VCF_DIR=$2
REFERENCE_SDF=$3
BED_FILE_DIR=$4
OUTPUT_DIR=$5

# Path to RTG Tools jar file
RTG_JAR="<path_to_jar_file>/RTG.jar"

# Ensure VCF files are block compressed and tabix indexed
echo "Indexing VCF files..."
for VCF in ${INPUT_VCF_DIR}/*.vcf.gz; do
    tabix -p vcf "$VCF"
done

# RTG vcfeval command
echo "Running RTG vcfeval..."
for VCF in ${INPUT_VCF_DIR}/*.vcf.gz; do
    for BED in ${BED_FILE_DIR}/*.bed; do
        VCF_NAME=$(basename ${VCF})
        BED_NAME=$(basename ${BED})
        OUTPUT_PATH="${OUTPUT_DIR}/${VCF_NAME}_vs_${BED_NAME}"

        echo "Running vcfeval for ${VCF_NAME} against ${BED_NAME}"
        java -jar $RTG_JAR vcfeval -b "$TRUTH_VCF" -c "$VCF" -t "$REFERENCE_SDF" -e "$BED" -o "$OUTPUT_PATH"
        echo "Completed vcfeval for ${VCF_NAME}"
    done
done

echo "RTG Tools vcfeval process completed."

