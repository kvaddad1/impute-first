## Upstream Personalization Workflow

This directory provides the Snakemake workflow for generating personalized diploid references using genotyping and imputation.

### Files and Subfolders

- `Snakefile`: Main Snakemake file that connects alignment, genotyping (`genotyping.smk`), and imputation (`beagle_imputation.smk` or `glimpse1_imputation.smk`) steps based on the config file.
- `genotyping.smk`: Performs pre-genotyping alignment with Bowtie2 or BWA and generates rough genotyping calls using `bcftools`.
- `beagle_imputation.smk`: Imputation workflow using Beagle.
- `glimpse1_imputation.smk`: Imputation workflow using GLIMPSE v1.
- `env.yml`: Conda environment file listing required software dependencies.

- `configs/`: Configuration directory.
  - `exampleData_config.yaml`: Sample config for running the workflow on chr21 of HG002.

- `download_exampleData.sh`: Downloads chr21 FASTQs and reference for HG002.
- `download_linkage_maps.sh`: Downloads genetic maps used by Beagle and GLIMPSE.
- `HGSVC3_reference_panel/`: HGSVC3-based reference panel VCF(s) used for imputation when running on full chromosomes or custom datasets.

### Example Configuration (`configs/exampleData_config.yaml`)

Below is a snippet of the configuration used for the provided example:

```yaml
sample_name: "HG002"
read_cov: 30
fq1: "exampleData/HG002_chr21_R1.fastq.gz"
fq2: "exampleData/HG002_chr21_R2.fastq.gz"
fa: "exampleData/GRCh38_no_alt_chr21.fasta"
covs: ["5"]
aligner: "bowtie2"
align_mode: "unpaired"
bcftools_mode: "filter"
bcftools_filter_QUAL: 20
bcftools_filter_DP: 100
impute_mode: "beagle"
ref_panel: "exampleData/hsvc.all.less1.chr21.fixed.vcf.gz"
```

### Notes

- For running instructions, see the top-level [README.md](../README.md).
- The configuration can be modified to analyze other samples, chromosomes, or reference panels.

