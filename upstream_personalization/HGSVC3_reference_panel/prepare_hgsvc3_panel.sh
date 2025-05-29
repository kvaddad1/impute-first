#!/bin/bash
set -euo pipefail

# Download HGSVC3 VCF files (GRCh38, version 1.0)
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_indel_insdel_alt_HGSVC2024v1.0.vcf.gz
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_snv_snv_alt_HGSVC2024v1.0.vcf.gz
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz

# Index individual VCFs using bcftools
bcftools index --threads 48 variants_GRCh38_indel_insdel_alt_HGSVC2024v1.0.vcf.gz
bcftools index --threads 48 variants_GRCh38_snv_snv_alt_HGSVC2024v1.0.vcf.gz
bcftools index --threads 48 variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz

# Concatenate, sort, normalize, and compress the merged VCF
bcftools concat -Ou -a \
    variants_GRCh38_indel_insdel_alt_HGSVC2024v1.0.vcf.gz \
    variants_GRCh38_snv_snv_alt_HGSVC2024v1.0.vcf.gz \
    variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz | \
    bcftools sort -Ou | \
    bcftools norm -m+any -Oz -o hsvc.all.vcf.gz

# Index the merged VCF
bcftools index --threads 32 hsvc.all.vcf.gz

# Fill in AC (allele count) and AN (allele number) INFO tags
bcftools +fill-tags hsvc.all.vcf.gz -Oz -o hsvc3.all.vcf.gz -- -t AC,AN
bcftools index --threads 48 hsvc3.all.vcf.gz

# Fix the VCF (custom script must be provided separately)
python3 fix_VCF.py hsvc3.all.vcf.gz hsvc3.all.fixed.vcf.gz
bcftools index --threads 48 hsvc3.all.fixed.vcf.gz

# Create HG002-excluded panel (useful for trio-based analyses)
bcftools view -s ^NA24385 hsvc3.all.fixed.vcf.gz -Oz -o hsvc3.all.less1.fixed.vcf.gz
bcftools index --threads 48 hsvc3.all.less1.fixed.vcf.gz

