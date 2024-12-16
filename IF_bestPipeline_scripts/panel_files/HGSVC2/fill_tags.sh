bcftools +fill-tags chm13_cactus_filtered_ids.phase_fixed.vcf.gz -Oz -o chm13_cactus_filtered_ids.fixed.vcf.gz -- -t AC,AN 
bcftools index --threads 48 chm13_cactus_filtered_ids.fixed.vcf.gz
echo "done filling tags"
