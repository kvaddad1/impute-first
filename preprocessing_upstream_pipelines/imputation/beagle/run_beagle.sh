#!/bin/bash
# Script: run_beagle.sh
# This script runs Beagle for genotype imputation chromomsome-wise. It requires the chromosome number,
# input VCF file, reference panel VCF file, output VCF file name, linkage map file, 
# and the number of threads to be used.

# Usage:
# ./run_beagle.sh <chromosome_number> <input_vcf_file> <ref_panel_vcf_file> <output_vcf_file> <linkage_map_file> <number_of_threads>

# Checking if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <chromosome_number> <input_vcf_file> <ref_panel_vcf_file> <output_vcf_file> <linkage_map_file> <number_of_threads>"
    exit 1
fi

# Assigning command-line arguments to variables
chromosome_number="$1"
input_vcf_file="$2"
ref_panel_vcf_file="$3"
output_vcf_file="$4"
linkage_map_file="$5"
number_of_threads="$6"

# Running Beagle for genotype imputation
echo "Running Beagle for genotype imputation..."
java -jar beagle.18May20.d20.jar chrom=chr${chromosome_number} impute=true gt=${input_vcf_file} ref=${ref_panel_vcf_file} out=${output_vcf_file} map=${linkage_map_file} nthreads=${number_of_threads}
echo "Beagle processing complete."

