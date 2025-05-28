

picard="picard.jar"
graph="hprc-v1.1-mc-grch38"
threads=48  

# based on the scripts on https://github.com/vgteam/vg/wiki/Haplotype-Sampling and https://github.com/vgteam/vg_wdl/blob/master/workflows/haplotype_sampling.wdl


vg index --threads ${threads} -j ${graph}.dist ${graph}.gbz   # In practice this step uses 1 CPU.
vg gbwt -p --num-threads ${threads} -r ${graph}.ri -Z ${graph}.gbz
vg haplotypes  -v 2 -t ${threads} -r ${graph}.ri  -H ${graph}.hapl -d ${graph}.dist  ${graph}.gbz   



mkdir tmp_kmc
echo ${sample}.R1.fq.gz > list_read_files.txt 
echo ${sample}.R2.fq.gz >> list_read_files.txt
kmc -k29 -m128 -okff -t${threads} -hp @list_read_files.txt  tmp_kmc/sample_kmc tmp_kmc/



kff="tmp_kmc/sample_kmc.kff"



vg haplotypes --include-reference --diploid-sampling -k $kff -i ${graph}.hapl  -v 2 -t ${threads} -r ${graph}.ri  -g ${graph}_haplotype_sampled_graph.gbz  ${graph}.gbz   


vg giraffe  --progress -t ${threads}  --read-group "ID:1 LB:lib1 SM:sample PL:illumina PU:unit1" --sample "sample"  \
		--output-format gaf -f ${sample}.R1.fq.gz -f ${sample}.R2.fq.gz -Z ${graph}_haplotype_sampled_graph.gbz  | gzip > aligned.gaf.gz







vg paths -L -S GRCh38 -x ${graph}_haplotype_sampled_graph.gbz | sort > path_list.txt
grep -v _decoy path_list.txt | grep -v _random |  grep -v chrUn_ | grep -v chrEBV | grep -v chrM | grep -v chain_ > path_list.sub.txt

vg paths --extract-fasta -p path_list.sub.txt --xg ${graph}_haplotype_sampled_graph.gbz > ref.fa

in_prefix_to_strip="GRCh38#0#"
sed -e "s/>${in_prefix_to_strip}/>/g" ref.prefix.fa > ref.fa
samtools faidx ref.fa
samtools dict ref.fa > ref.dict


bam_file="aligned.bam"

vg surject  -F "path_list.sub.txt"  -x ${graph}_haplotype_sampled_graph.gbz  -t ${threads} --bam-output --gaf-input  --sample sample  \
		--read-group "ID:1 LB:lib1 SM:sample PL:illumina PU:unit1"  --prune-low-cplx --interleaved --max-frag-len 3000  aligned.gaf.gz > ${bam_file}



samtools view -H ${bam_file} | grep ^@HD > new_header.sam
grep ^@SQ ref.dict | awk '{print $1 "\t" $2 "\t" $3}' >> new_header.sam
samtools view -H ${bam_file}  | grep -v ^@HD | grep -v ^@SQ >> new_header.sam
cat <(cat new_header.sam) <(samtools view ${bam_file}) | sed -e "s/${in_prefix_to_strip}//g" |  samtools sort --threads ${threads} -O BAM > aligned_38_sorted.bam

samtools index aligned_38_sorted.bam

