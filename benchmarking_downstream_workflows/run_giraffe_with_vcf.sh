#!/bin/bash
# Script: run_giraffe_with_vcf.sh
# This script builds Giraffe index and maps reads using a VCF and reference genome.
# Note:
# - If input VCF is a "truth VCF", the workflow corresponds to Giraffe (Benchmark).
# - If input VCF is "1kGP phase VCF", the workflow corresponds to Pangenome Giraffe(1kGP_Pangenome).
#
# Usage:
# ./run_giraffe_with_vcf.sh -r <reference.fa> -v <variants.vcf.gz> -p <prefix> -i <index_dir> -d <reads_dir> -b <bam_prefix> -o <output_dir> -t <threads>

while getopts "r:v:p:i:d:b:o:t:" opt; do
  case "$opt" in
    r) ref=$OPTARG ;;
    v) vcf=$OPTARG ;;
    p) prefix=$OPTARG ;;
    i) index_dir=$OPTARG ;;
    d) reads_dir=$OPTARG ;;
    b) bam_prefix=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    t) threads=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

if [ -z "$ref" ] || [ -z "$vcf" ] || [ -z "$prefix" ] || [ -z "$index_dir" ] || [ -z "$reads_dir" ] || [ -z "$bam_prefix" ] || [ -z "$out_dir" ] || [ -z "$threads" ]; then
  echo "Usage: $0 -r <ref.fa> -v <vcf.gz> -p <prefix> -i <index_dir> -d <reads_dir> -b <bam_prefix> -o <out_dir> -t <threads>"
  exit 1
fi

mkdir -p "$index_dir" "$out_dir"

echo "Running vg autoindex..."
/bin/time -o "${index_dir}/${prefix}_autoindex.stats" --format='user= %U system= %S elapsed= %e CPU= %P MemMax= %M' \
  vg autoindex --workflow giraffe -r "$ref" -v "$vcf" -p "${index_dir}/${prefix}" --threads "$threads"

echo "Mapping reads with vg giraffe..."
/bin/time -o "${out_dir}/run_map_giraffe.stats" --format='user=%U system=%S elapsed=%e CPU=%P MemMax=%M' \
  vg giraffe -Z "${index_dir}/${prefix}.giraffe.gbz" -m "${index_dir}/${prefix}.min" -d "${index_dir}/${prefix}.dist" \
  -f "${reads_dir}/read1.fastq.gz" -f "${reads_dir}/read2.fastq.gz" \
  -p -o BAM -t "$threads" > "${out_dir}/${bam_prefix}.bam"

echo "done bam generation"
