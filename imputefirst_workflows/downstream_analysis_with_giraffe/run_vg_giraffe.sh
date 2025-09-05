#!/bin/bash

# Usage: run_vg_giraffe.sh <index_directory> <reads_directory> <bam_filename> <output_directory> <threads> <prefix>
# Example: run_vg_giraffe.sh /path/to/index /path/to/reads sample.bam /path/to/output 16 sample_prefix

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <index_directory> <reads_directory> <bam_filename> <output_directory> <threads> <prefix>"
    exit 1
fi

# Assign arguments to variables
path_index=$1
path_reads=$2
bam_filename=$3
output_directory=$4
threads=$5
prefix=$6

# Mapping Step
echo "Started mapping with vg giraffe..."
date
/bin/time -o "${output_directory}/run_map_giraffe.stats" --format='user=%U system=%S elapsed=%e CPU=%P MemMax=%M' \
    vg giraffe -Z "${path_index}/${prefix}.giraffe.gbz" -m "${path_index}/${prefix}.min" -d "${path_index}/${prefix}.dist" \
    -f "${path_reads}/read1.fastq.gz" -f "${path_reads}/read2.fastq.gz" \
    -p -o BAM -t ${threads} > "${output_directory}/${bam_filename}.bam"
date
echo "Done mapping with vg giraffe."

