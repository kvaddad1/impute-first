#!/bin/bash

# Create ref_panel directory if it doesn't exist
mkdir -p ref_panel

# Download the VCF files into the ref_panel folder
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_snv_snv_alt.vcf.gz -O ref_panel/variants_freeze4_snv_snv_alt.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_indel_insdel_alt.vcf.gz -O ref_panel/variants_freeze4_indel_insdel_alt.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v2.0/integrated_callset/variants_freeze4_sv_insdel_alt.vcf.gz -O ref_panel/variants_freeze4_sv_insdel_alt.vcf.gz

bcftools index ref_panel/variants_freeze4_snv_snv_alt.vcf.gz
bcftools index ref_panel/variants_freeze4_indel_insdel_alt.vcf.gz 
bcftools index ref_panel/variants_freeze4_sv_insdel_alt.vcf.gz

# Concatenate, normalize, compress
bcftools concat -a -Ou ref_panel/variants_freeze4_indel_insdel_alt.vcf.gz ref_panel/variants_freeze4_snv_snv_alt.vcf.gz ref_panel/variants_freeze4_sv_insdel_alt.vcf.gz | \
bcftools norm -m+any -Ou - | \
bcftools view -Oz -o ref_panel/hsvc.all.vcf.gz

# Index the gzipped VCF file
bcftools index ref_panel/hsvc.all.vcf.gz

echo "done downloading and aggregating reference panel HGSVC2"
