#!/bin/bash
# Script: gen_bwa_index.sh
# This script is used for indexing Haplotype_1 and Haplotype_2 fasta files of a personalized diploid reference 
# using BWA. The Haplotype fasta files can be generated by the script gen_diploid_ref.sh in this directory.   
#
# Usage:
# ./gen_bwa_index.sh <fasta_folder_path> <file_prefix>
#
# Arguments:
# fasta_folder_path - Path to the folder containing the output fasta files.
# file_prefix - Prefix used in naming the output files from diploid_ref.sh.

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_folder_path> <file_prefix>"
    exit 1
fi

# Assigning command-line arguments to variables for clarity
fasta_folder_path="$1" # Path to folder containing the output fasta files
file_prefix="$2"       # Prefix used in naming the output files

echo "Starting BWA indexing for personalized diploid reference!"
date

# Indexing consensus for Haplotype_1 fasta with BWA
echo "Indexing consensus for Haplotype_1 fasta with BWA"
bwa index "${fasta_folder_path}/${file_prefix}.Haplotype_1.fa"

# Indexing consensus for Haplotype_2 fasta with BWA
echo "Indexing consensus for Haplotype_2 fasta with BWA"
bwa index "${fasta_folder_path}/${file_prefix}.Haplotype_2.fa"

date
echo "Done with BWA indexing for personalized diploid reference!"