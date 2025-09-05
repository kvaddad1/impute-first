## Reference Construction Module

Scripts to construct graph-based references.

### Scripts
- `gen_vg_autoindex_VCF.sh`: Utilizes vg autoindex for indexing genomic data, essential for constructing personalized reference graphs.

### Personalized Graph Reference
The `gen_vg_autoindex_VCF.sh` script is crucial for generating personalized reference graphs. This process involves using vg autoindex with the giraffe workflow, as detailed in the study. The script uses the GRCh38 reference FASTA and personalized diploid variant calls (VCF file), embedding unique variants into the graph. 

### Installation Prerequisites
- **VG Autoindex and VG Giraffe (v1.55.0)**:
  - For installation, visit the [VG GitHub Repository](https://github.com/vgteam/vg/releases).

### Usage
- Follow the installation guidelines for vg.
- Use `gen_vg_autoindex_VCF.sh` for indexing the references.

Refer to the comments within each script for detailed usage instructions. The indexes generated are prerequisites for the processes in the alignment and variant calling phases.


---

## Downstream Analysis Module

This module contains scripts used for mapping and variant calling in genomic data analysis, focusing on VG Giraffe-based mapping, followed by variant calling with GATK HaplotypeCaller and DeepVariant.

### Scripts
- `run_vg_giraffe.sh`: Executes the VG Giraffe mapping on paired-end reads using Giraffe indexes.
- `variant_calling/run_GATK-HC.sh`: Runs the GATK HaplotypeCaller on the BAM outputs from mapping pipelines to call variants.
- `variant_calling/run_deepvariant.sh`: Runs DeepVariant using a Singularity container to perform WGS variant calling from aligned BAM files.

### Installation Prerequisite for GATK HaplotypeCaller
- **GATK HaplotypeCaller**:
  - Installation and setup instructions can be found on the [GATK Getting Started Guide](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).
  - Ensure that Java is installed as it is a requirement for running GATK.
  - Download and configure GATK, following the detailed steps provided in the guide.

### Installation Prerequisite for DeepVariant
- **DeepVariant (via Singularity)**:
  - Install [Singularity](https://sylabs.io/guides/) on your system.
  - Pull the DeepVariant image (e.g., version 1.5.0) using:
    ```bash
    BIN_VERSION="1.5.0"
    singularity pull docker://google/deepvariant:"${BIN_VERSION}"
    ```
  - This will create a Singularity image file named `deepvariant_1.5.0.sif` in your working directory.

### Usage
- Make sure all prerequisite tools are installed and configured.
- Utilize `run_vg_giraffe.sh` for generating BAM files through different alignment strategies.
- Use `variant_calling/run_GATK-HC.sh` or `variant_calling/run_deepvariant.sh` to perform variant calling on BAM files using GATK HaplotypeCaller or DeepVariant, respectively.

