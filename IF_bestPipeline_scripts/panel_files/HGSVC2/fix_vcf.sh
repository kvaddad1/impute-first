date;
python3 fix_VCF.py ref_panel/hsvc.all.vcf.gz ref_panel/hsvc.all.fixed.vcf.gz 
date;
echo "done fixing"
bcftools index ref_panel/hsvc.all.fixed.vcf.gz 
echo "done indexing the fixed vcf file"
