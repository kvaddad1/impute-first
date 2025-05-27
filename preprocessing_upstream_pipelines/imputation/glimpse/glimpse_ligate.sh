#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <ligate_dir> <impute_dir> <output_file>"
    exit 1
fi

LIG_DIR=$1
IMP_DIR=$2
outputfile=$3
list_name=list.${outputfile}.txt
pathh="glimpse_static_bins"
PROC=32  # Default value, can be overridden by passing as an argument

mkdir -p ${LIG_DIR}
LST=${LIG_DIR}/${list_name}
ls ${IMP_DIR}/${outputfile}*bcf > ${LST}
OUT=${LIG_DIR}/${outputfile}.merged.bcf

echo "Starting ligation..."
${pathh}/GLIMPSE_ligate_static --thread ${PROC} --input ${LST} --output ${OUT}
bcftools index -f ${OUT}
bcftools view ${OUT} -Oz -o ${OUT}.vcf.gz 
bcftools index -f ${OUT}.vcf.gz
echo "Ligation complete."

