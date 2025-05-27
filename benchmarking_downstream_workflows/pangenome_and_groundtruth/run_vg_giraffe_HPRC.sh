#!/bin/bash
# Usage: ./run_giraffe_hprc.sh <sample_id> <fq1> <fq2> <graph_dir> <threads>

set -e

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <sample_id> <fq1> <fq2> <graph_dir> <threads>"
    exit 1
fi

sample="$1"
fq1="$2"
fq2="$3"
graph_dir="$4"
threads="$5"

vg="./vg"
graph="${graph_dir}/hprc-v1.1-mc-grch38.d9.gbz"
ref_dict="${graph_dir}/ref.dict"
gaf_out="${sample}_aligned.gaf.gz"
bam_out="${sample}_aligned_38_sorted.bam"
header_out="${sample}_new_header.sam"

# Giraffe mapping
${vg} giraffe -t ${threads} \
  --read-group "ID:1 LB:lib1 SM:${sample} PL:illumina PU:unit1" \
  --sample "${sample}" \
  --output-format gaf \
  -f "${fq1}" -f "${fq2}" \
  -Z "${graph}" | gzip > "${gaf_out}"

# Surject GAF to BAM
${vg} surject -x "${graph}" -t "${threads}" \
  --bam-output --gaf-input \
  --sample "${sample}" \
  --read-group "ID:1 LB:lib1 SM:${sample} PL:illumina PU:unit1" \
  --prune-low-cplx --interleaved --max-frag-len 3000 \
  "${gaf_out}" > aligned.bam

# fixing BAM header and sorting
{
  samtools view -H aligned.bam | grep ^@HD
  grep ^@SQ "${ref_dict}" | awk '{print $1 "\t" $2 "\t" $3}'
  samtools view -H aligned.bam | grep -v ^@HD | grep -v ^@SQ
} > "${header_out}"

prefix_strip="GRCh38#0#"
cat "${header_out}" <(samtools view aligned.bam) \
  | sed "s/${prefix_strip}//g" \
  | samtools sort --threads "${threads}" -O BAM > "${bam_out}"

echo "Done generating bam ${bam_out}"

