## Context 
This module is the benchmarking component as described in the [rowbowt paper](https://pubmed.ncbi.nlm.nih.gov/36409181/). It is specifically designed for generating genotyping calls for `rowbowt`, `bayestyper`, and `bowtie2+bcftools`. Within the Impute-first Alignment framework, this module functions as the preprocessing module, integrating read sampling, alignment, and rough genotyping into a Snakemake workflow. This workflow utilizes HG001 read data and HGSVC2 reference panels, aligning with the methodologies discussed in the paper.

While this module is an integral part of the modular Impute-first Alignment framework, its use is optional. It is particularly useful to guide on generating genotype calls using their method of choice before plugging into the imputation module of the framework.

## Installation
1) download and install [rowbowt](https://github.com/alshai/rowbowt) somewhere (please use the `dev` branch!)

```
git clone --recursive git@github.com:alshai/rowbowt.git .
cd rowbowt
git checkout dev
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_DIRECTORY=<install dir>
cmake --build . && cmake --install .
export PATH=<install dir>:$PATH # optional, if <install dir> not already in PATH
```

2) clone this repo:

```git clone --recursive git@github.com:alshai/impute_first.git```

3) in your config file, modify `pfbwt_bins` and `rowbowt_bins` to point to their respective executable paths

4) modify the other parts of the config file accordingly

5) the Snakemake should handle most of the other files and executables needed

## File structure:

1) `Snakemake` - the "main" Snakemake file. contains common utility rules

2) `rowbowt.snk` - contains rules used to build rowbowt index

3) `marker_benchmark.snk` - contains rules used for benchmarking other tools (bowtie2_bcftools, bayestyper)

4) `configs/*.yaml` - an assortment of configurations for running the workflow.

    - ex. `configs/hg38.hsvc.yaml` is the config for building a rowbowt index using the variants provided by HSVC (SVs + short variants).

    - modify accordingly to suit your input reads/reference files. 

## Usage
Modify `rule all` in `Snakefile` to produce different results:

```
rule all:
    input:
        """ FNAME """
```

possible  substitutes for FNAME:
- `ROWBOWT_FNAME` - builds the rowbowt index on yaml file based data files 

Then, run (in this directory):
```
snakemake -j<nthreads> --configfile <your configfile> 
```
