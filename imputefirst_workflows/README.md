## Repository Structure

- `upstream_personalization/`: Snakemake workflow to generate personalized diploid references from downsampled reads by performing genotyping and imputation using a reference panel.


- `downstream_analysis_with_giraffe/`: Downstream analysis using Impute-First personalized diploid references of a given sample, built into variation graphs, enabling alignment and lift-over with VG Giraffe (Figure 1, Section B.1) and to perform variant calling.

- `downstream_analysis_with_leviosam2/`: Downstream analysis using Impute-First personalized diploid FASTA references of a given sample, aligned with BWA-MEM and lifted with LevioSAM2 to perform variant calling (Figure 1, Section B.2).

## Running the Upstream Personalization Workflow

The upstream personalization workflow is implemented using Snakemake and defined in `upstream_personalization/`. It performs read downsampling, genotyping, and imputation to generate diploid personalized reference for a given sample.

All required software dependencies are listed in the Conda environment file `env.yml`. A compatible environment can be created with:

```bash
conda env create -f imputefirst_workflows/upstream_personalization/env.yml
conda activate genotyping_imputation
```

To test the workflow, a chromosome 21–based demonstration is provided in the `upstream_personalization/` directory. It includes the following files:
- `download_exampleData.sh`: Downloads a small test dataset (HG002 chr21 reads and reference files).
- `download_linkage_maps.sh`: Downloads linkage maps for use with Beagle and GLIMPSE.
- `configs/exampleData_config.yaml`: Configuration file for running the workflow on the example dataset.

Run the workflow with:
```bash
snakemake -j <threads> --configfile configs/exampleData_config.yaml
```

This config file parameters can be modified to enable using different input files.  

## Running the Downstream Analysis Workflow
[LevioSAM2 based workflow](downstream_analysis_with_leviosam2/)

[VG giraffe based workflow](downstream_analysis_with_giraffe/)
