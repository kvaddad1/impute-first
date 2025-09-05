mkdir -p resources
# download docker images
singularity pull ./resources/deepvariant_1.6.0.sif docker://google/deepvariant:1.6.0  
singularity pull ./resources/leviosam2_v0.5.0.sif docker://naechyun/leviosam2:v0.5.0 

# download reference genome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -P ./resources/
gunzip ./resources/chm13v2.0.fa.gz
wget --content-disposition --trust-server-names --no-check-certificate \
https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED
unzip ncbi_dataset.zip -d ./resources/

# download chain file
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain -P ./resources/

# make index files
samtools faidx ./resources/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
singularity exec ./resources/leviosam2_v0.5.0.sif leviosam2 index -c ./resources/chm13v2-grch38.chain -p ./resources/chm13v2-grch38 -F ./resources/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.fai
bwa index ./resources/chm13v2.0.fa
