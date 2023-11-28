#!/bin/bash
# Usage: $0 <input_bam> <benchmark_vcf> <reference_fasta> <output_directory> <balance_value> <output_pdf>

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <input_bam> <benchmark_vcf> <reference_fasta> <output_directory> <balance_value> <output_pdf>"
    exit 1
fi

INPUT_BAM=$1
BENCHMARK_VCF=$2
REFERENCE_FASTA=$3
OUTPUT_DIR=$4
BALANCE_VALUE=$5
OUTPUT_PDF=$6

# Splitting BAM by chromosome and running Biastools
for i in $(seq 1 22); do
    samtools view -h "$INPUT_BAM" "chr$i" > "${OUTPUT_DIR}/chr${i}_sorted.sam"
    echo "Processing chromosome ${i}"
    python3 ref_bi_adaptive_wgs.py -v "$BENCHMARK_VCF" -s "${OUTPUT_DIR}/chr${i}_sorted.sam" -f "$REFERENCE_FASTA" -r -o "${OUTPUT_DIR}/out_${i}"
done

# Aggregating results
for i in $(seq 1 22); do 
    cat "${OUTPUT_DIR}/out_${i}" >> "${OUTPUT_DIR}/agg_chr.txt"
done

# Plotting results
echo "Starting Plotting"
python3 plot_biastools.py -lr "${OUTPUT_DIR}/agg_chr.txt" -ln "YourLabelsHere" -vcf "$BENCHMARK_VCF" -bd "$BALANCE_VALUE" -real -out "$OUTPUT_PDF"
echo "Plotting completed"

echo "Biastools analysis and plotting done."
