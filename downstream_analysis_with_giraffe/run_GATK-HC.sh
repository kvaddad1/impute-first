#!/bin/bash

# Usage: run_GATK-HC.sh <bam_filename> <output_directory> <reference_file> <threads>
# Example: run_GATK-HC.sh sample.bam /path/to/output /path/to/reference.fasta 16

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <bam_filename> <output_directory> <reference_file> <threads>"
    exit 1
fi

# Assign arguments to variables
bam_filename=$1
output_directory=$2
ref=$3
threads=$4

# Define file paths
sorted_bam="${output_directory}/sorted_${bam_filename}.bam"
sorted_rg_bam="${output_directory}/sorted_rg_${bam_filename}.bam"

# Sorting the BAM file
echo "Started sorting the BAM file..."
date
samtools sort -@${threads} "${bam_filename}" > "$sorted_bam" || { echo "Sorting failed."; exit 1; }
date
echo "Done sorting the BAM file."

# Indexing the sorted BAM file
echo "Started indexing the sorted BAM file..."
date
samtools index -@${threads} "$sorted_bam" || { echo "Indexing failed."; exit 1; }
date
echo "Done indexing the sorted BAM file."

# Adding ReadGroup to the sorted BAM
echo "Started adding ReadGroup to the BAM..."
date
samtools addreplacerg -r '@RG\tID:grch38\tSM:sample\tPL:ILLUMINA\tDS:novaseq\tPU:novaseq' "$sorted_bam" -o "$sorted_rg_bam" || { echo "Adding ReadGroup failed."; exit 1; }
date
echo "Done adding ReadGroup to the BAM."

# Indexing the BAM file with ReadGroup
echo "Started indexing the BAM file with ReadGroup..."
date
samtools index -@${threads} "$sorted_rg_bam" || { echo "Indexing failed."; exit 1; }
date
echo "Done indexing the BAM file with ReadGroup."

# Run GATK HaplotypeCaller to call variants
echo "Started GATK HaplotypeCaller..."
date
/bin/time -o "${output_directory}/run_gatk.stats" --format='user=%U system=%S elapsed=%e CPU=%P MemMax=%M' \
    gatk HaplotypeCaller --native-pair-hmm-threads ${threads} -R "$ref" -I "$sorted_rg_bam" -O "${output_directory}/${bam_filename}_variants.vcf" || { echo "GATK HaplotypeCaller failed."; exit 1; }
date
echo "Done with GATK HaplotypeCaller."
