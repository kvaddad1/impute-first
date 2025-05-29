#!/bin/bash

# Usage:
# ./run_happy.sh <truth_vcf> <query_vcf> <truth_bed> <reference_fasta> <output_dir> <hap_py_sif> <rtg_sdf_template> <threads> [stratification_tsv]

truth_vcf=$1
query_vcf=$2
truth_bed=$3
ref_fasta=$4
output_dir=$5
hap_py_sif=$6
rtg_sdf_template=$7
threads=$8
strat_tsv=$9  # Optional

mkdir -p "$output_dir"

cmd=(
  singularity exec "$hap_py_sif" /opt/hap.py/bin/hap.py
  "$truth_vcf" "$query_vcf"
  --engine vcfeval
  -o "${output_dir}/out"
  -r "$ref_fasta"
  -f "$truth_bed"
  --threads "$threads"
  --pass-only
  --engine-vcfeval-path /path/to/rtg-tools
  --engine-vcfeval-template "$rtg_sdf_template"
)

# Add stratification if provided
if [ -n "$strat_tsv" ]; then
  cmd+=(--stratification "$strat_tsv")
fi

"${cmd[@]}"

