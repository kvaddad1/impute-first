#downloading the genetic map files and the panel files - HGSVC2
./download_linkage_maps.sh  
./download_panel_files.sh

#optional for filling in AC, AN info - if missing 
#./fill_tags.sh

#fixes the missing GT or unphased info 
./fix_vcf.sh
echo "done wrap run for fixing the VCF"

# Remove sample NA12878 or NA24385 [based on the choice of Leave-One-Out analysis using the full panel file: ref_panel/hsvc.all.fixed.vcf.gz]

#Removing HG001 from panel 
bcftools view -s ^NA12878 ref_panel/hsvc.all.fixed.vcf.gz -Oz -o ref_panel/hsvc.all.less1.HG001.fixed.vcf.gz 

#Removing HG002 from panel 
bcftools view -s ^NA24385 ref_panel/hsvc.all.fixed.vcf.gz -Oz -o ref_panel/hsvc.all.less1.HG002.fixed.vcf.gz 

#indexing the panel-specific HGSVC2-LOO files for HG001/HG002 respectively. 
bcftools index ref_panel/hsvc.all.less1.HG001.fixed.vcf.gz 
bcftools index ref_panel/hsvc.all.less1.HG002.fixed.vcf.gz 

echo "Done preparing the HGSVC2-LOO panel for NA24385 (HG002)"
