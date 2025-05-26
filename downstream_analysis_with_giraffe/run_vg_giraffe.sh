#!/bin/bash

# Usage: run_vg_giraffe.sh <index_directory> <reads_directory> <bam_filename> <output_directory> <threads> <prefix> <fragment_mean> <fragment_stdev>
# Example: run_vg_giraffe.sh /path/to/index /path/to/reads sample.bam /path/to/output 16 sample_prefix 350 50

# Check if the correct number of arguments is provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <index_directory> <reads_directory> <bam_filename> <output_directory> <threads> <prefix> <fragment_mean> <fragment_stdev>"
    exit 1
fi

# Assign arguments to variables
path_index=$1
path_reads=$2
bam_filename=$3
output_directory=$4
threads=$5
prefix=$6
fragment_mean=$7
fragment_stdev=$8

# Mapping Step
echo "Started mapping with vg giraffe..."
date
/bin/time -o "${output_directory}/run_map_giraffe.stats" --format='user=%U system=%S elapsed=%e CPU=%P MemMax=%M' \
    vg giraffe -Z "${path_index}/${prefix}.giraffe.gbz" -m "${path_index}/${prefix}.min" -d "${path_index}/${prefix}.dist" \
    -f "${path_reads}/read1.fastq.gz" -f "${path_reads}/read2.fastq.gz" \
    --fragment-mean ${fragment_mean} --fragment-stdev ${fragment_stdev} -p -o BAM -t ${threads} > "${output_directory}/${bam_filename}.bam"
date
echo "Done mapping with vg giraffe."
