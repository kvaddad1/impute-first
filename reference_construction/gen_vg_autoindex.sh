#!/bin/bash
# Script Name: gen_vg_autoindex.sh
# This script is used to generate vg graph based index using giraffe pipeline using a VCF file and a reference genome.
#
# Usage:
# ./gen_vg_autoindex.sh -r <reference_file> -v <variant_file> -p <prefix_for_output_files> -t <number_of_threads>
#
# Arguments:
# -r <reference_file>       : Path to the reference genome file.
# -v <variant_file>         : Path to the variant file (VCF format).
# -p <prefix_for_output_files> : Prefix for naming the output files.
# -t <number_of_threads>    : Number of threads to use for processing.

# Initialize variables
ref=""
vars=""
prefix=""
threads=""

# Parse command-line arguments
while getopts "r:v:p:t:" flag; do
    case "${flag}" in
        r) ref=${OPTARG};;
        v) vars=${OPTARG};;
        p) prefix=${OPTARG};;
        t) threads=${OPTARG};;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

# Check if the necessary arguments are provided
if [ -z "$ref" ] || [ -z "$vars" ] || [ -z "$prefix" ] || [ -z "$threads" ]; then
    echo "Usage: $0 -r <reference_file> -v <variant_file> -p <prefix_for_output_files> -t <number_of_threads>"
    echo "Error: Missing arguments. Please provide all the required arguments."
    exit 1
fi

# VG Indexing Step
echo "Started vg autoindex"
date
/bin/time -o "${prefix}_autoindex.stats" --format='user= %U system= %S elapsed= %e CPU= %P MemMax= %M' vg autoindex --workflow giraffe -r "$ref" -v "$vars" -p "$prefix" --threads "$threads"
date
echo "Done with vg autoindex"

