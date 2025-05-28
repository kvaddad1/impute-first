This directory contains scripts to generate BAM files for a given sample reads (fastqs), in GRCh38 coordinates using different reference indexing and alignment workflows.

- `run_bwa.sh`  
  Linear alignment using BWA-MEM with GRCh38 reference.

- `run_vg_giraffe_HPRC.sh`  
  Downloads frequency-filtered HPRC v1.1 pangenome indexes and performs mapping with `vg giraffe`.

- `run_giraffe_with_vcf.sh`  
  Builds Giraffe index from input VCF and reference FASTA, then maps reads.  
  Used for both 1kGP pangenome and GIAB truth VCF workflows.

- `run_realignment_abra.sh`  
  Performs indel realignment using ABRA2 to improve alignment accuracy.

- `run_vg_personalized_pipeline.sh`  
  Builds a sample-specific personalized pangenome index by subsampling haplotypes from the HPRC v1.1 default pangenome indexes, and then maps reads.

