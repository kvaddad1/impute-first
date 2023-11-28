#!/bin/bash
# Script: gen_diploid_ref.sh
# This script generates a personalized diploid reference genome from a VCF file.
# It creates two haplotypes based on the input VCF and a reference genome.

# Usage:
# ./gen_diploid_ref.sh <reference_genome_fasta> <personalized_vcf_file> <output_file_prefix> <output_fasta_directory> <number_of_threads>
#
# Arguments:
# reference_genome_fasta - Path to the reference genome fasta (such as GRCh38.fasta)
# personalized_vcf_file - Path to the personalized variant call set (VCF) from Impute-first workflow. 
# output_file_prefix - Prefix for naming the output files.
# output_fasta_directory - Directory name for storing output fasta files.
# number_of_threads - Number of threads to use for processing.

# Check for correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <reference_genome_fasta> <personalized_vcf_file> <output_file_prefix> <output_fasta_directory> <number_of_threads>"
    exit 1
fi

# Assigning command-line arguments to variables for clarity
reference_genome_fasta="$1"    # GRCh38 reference genome fasta file
personalized_vcf_file="$2"     # Personalized variant call set (VCF) file
output_file_prefix="$3"        # Prefix for output files
output_fasta_directory="$4"    # Directory for output fasta files
processing_threads="$5"        # Number of processing threads

echo "Starting construction of personalized diploid reference!"
date

# Indexing the reference genome fasta file
echo "Indexing the reference genome fasta"
samtools faidx "$reference_genome_fasta"

# Indexing the personalized VCF file
echo "Indexing the personalized VCF file"
bcftools index -f "$personalized_vcf_file" --threads "$processing_threads"

# Normalizing the personalized VCF file
echo "Normalizing the personalized VCF"
bcftools norm -f "$reference_genome_fasta" "$personalized_vcf_file" -m +any -Oz -o "normalized_$personalized_vcf_file" --threads "$processing_threads"

# Sorting the normalized VCF file
echo "Sorting the normalized VCF"
bcftools sort "normalized_$personalized_vcf_file" -Oz -o "sorted_normalized_$personalized_vcf_file"

# Indexing the normalized and sorted personalized VCF file
echo "Indexing the normalized and sorted VCF file"
bcftools index -f "sorted_normalized_$personalized_vcf_file" --threads "$processing_threads"

# Generating consensus for Haplotype_1
echo "Generating consensus for Haplotype_1"
mkdir -p "$output_fasta_directory"
bcftools consensus -f "$reference_genome_fasta" -o "${output_fasta_directory}/${output_file_prefix}.Haplotype_1.fa" -H 1 "sorted_normalized_$personalized_vcf_file" -c "${output_fasta_directory}/${output_file_prefix}.ref2Haplotype_1.chain"

# Generating consensus for Haplotype_2
echo "Generating consensus for Haplotype_2"
bcftools consensus -f "$reference_genome_fasta" -o "${output_fasta_directory}/${output_file_prefix}.Haplotype_2.fa" -H 2 "sorted_normalized_$personalized_vcf_file" -c "${output_fasta_directory}/${output_file_prefix}.ref2Haplotype_2.chain"

# Indexing consensus for Haplotype_1 fasta
echo "Indexing consensus for Haplotype_1 fasta"
samtools faidx "${output_fasta_directory}/${output_file_prefix}.Haplotype_1.fa"

# Indexing consensus for Haplotype_2 fasta
echo "Indexing consensus for Haplotype_2 fasta"
samtools faidx "${output_fasta_directory}/${output_file_prefix}.Haplotype_2.fa"

date
echo "Done constructing personalized diploid reference!"

