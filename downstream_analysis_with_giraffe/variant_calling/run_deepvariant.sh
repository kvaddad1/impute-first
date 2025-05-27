#!/bin/bash

# Usage: run_deepvariant.sh <input_bam> <output_directory> <reference_fasta> <threads> <deepvariant_sif_path>
# Example: run_deepvariant.sh sample.bam /path/to/output hg38.no_alt.fa 32 /path/to/deepvariant_1.5.0.sif

# Check input arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input_bam> <output_directory> <reference_fasta> <threads> <deepvariant_sif_path>"
    exit 1
fi

# input arguments
input_bam="$1"
output_dir="$2"
reference_fasta="$3"
threads="$4"
deepvariant_sif="$5"

# output directory
mkdir -p "$output_dir"

# Run DeepVariant
echo "Running DeepVariant..."
/bin/time -o "${output_dir}/run_deepvariant.stats" --format='user=%U system=%S elapsed=%e CPU=%P MemMax=%M' \
singularity run "$deepvariant_sif" /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="$reference_fasta" \
    --reads="$input_bam" \
    --output_vcf="${output_dir}/output.vcf.gz" \
    --num_shards="$threads" \
    --dry_run=false \
    --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true"
echo "DeepVariant completed successfully."

