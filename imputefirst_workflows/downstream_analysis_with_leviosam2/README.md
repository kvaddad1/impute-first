# Downstream Analysis Pipeline with LevioSAM2

## Required software
- BWA MEM (0.7.17)
- samtools (1.21)
- singularity (3.11.4-1.el8)
- bcftools (1.22)

## Preparation
Run 
```
bash download_resource.sh
```
to download the prerequisite files listed as:
- docker image of deepvariant
- docker image of LevioSAM2
- reference genome GRCh38
- reference genome T2T-CHM13
- chain file from T2T-CHM13 to GRCh38

The script also generate the index files for the snakemake pipeline

## Running the snakemake pipeline
In the ```config.yaml```, at least three input path needs to be updated.
1. ```reads_r1```: read 1 sequence data file of the sample.
2. ```reads_r2```: read 2 sequence data file of the sample.
3. ```vcf```: the target personalized VCF file, i.e. the impute-first upstream generated VCF file. Though other VCF file can also be used.

4. ```vcf_name```: best practice to set as the prefix of the VCF file.

```chm13``` flag is set to True in default mode to include T2T-CHM13 in the analysis.
```invert.py``` is borrowed from [chaintools](https://github.com/milkschen/chaintools).

Run the full pipeline with
```
snakemake -j 32
```
