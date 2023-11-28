## Imputation Module: Beagle & GLIMPSE
This module is designed to perform genotype imputation using two powerful tools: Beagle and GLIMPSE. It has been utilized to process rough genotyped calls from three methods: bowtie2+bcftools, rowbowt, and bayestyper in the Impute-First Alignment Workflow. 

### Beagle Directory (`./beagle`)
- `run_beagle.sh`: Script to execute Beagle imputation.

### GLIMPSE Directory (`./glimpse`)
- `run_glimpse_impute.sh`: Main script for the GLIMPSE imputation process.
- `glimpse_ligate.sh`: Handles the ligation step in GLIMPSE.
- `glimpse_phase.sh`: Manages the phasing step in GLIMPSE.
- `glimpse_sample.sh`: Conducts the sampling step in GLIMPSE.
- `glimpse_chunk.sh`: Performs the chunking step in GLIMPSE.
- `prepare_mapfiles.sh`: Prepares genetic map files for GLIMPSE.

## Beagle Imputation

### Requirements
- Beagle version 5.1 (`beagle.18May20.d20.jar`)
- Java version 8

### Genetic Map Files
Download GrCh38 genetic maps in PLINK format with cM units from [Beagle Genetic Maps](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

## GLIMPSE Imputation

### Building GLIMPSE from Scratch
Follow the installation guide at [GLIMPSE Installation](https://odelaneau.github.io/GLIMPSE/glimpse1/installation.html).

### Using Static Binaries
Download from [GLIMPSE GitHub Releases](https://github.com/odelaneau/GLIMPSE/releases/tag/v1.1.1).

### Tutorial
For a step-by-step guide, visit [GLIMPSE Tutorial](https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_hg19.html).
