#!/bin/bash
# Usage: ./bwa_index_and_align.sh <ref_fasta> <fq1> <fq2> <output_bam> <sample_id> <threads>

set -e

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <ref_fasta> <fq1> <fq2> <output_bam> <sample_id> <threads>"
    exit 1
fi

ref_fasta="$1"
fq1="$2"
fq2="$3"
output_bam="$4"
sample_id="$5"
threads="$6"

echo "Indexing reference"
bwa index "$ref_fasta"

echo "Running BWA-MEM alignment"
bwa mem -t "$threads" -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" "$ref_fasta" "$fq1" "$fq2" | samtools view -b > "$output_bam"

