## Downstream Analysis Module

This module contains scripts used for mapping and variant calling in genomic data analysis, focusing on both linear mapping and VG Giraffe-based mapping, followed by variant calling with GATK HaplotypeCaller.

### Scripts
- `run_bwamem.sh`: Performs mapping on paired-end reads using BWA indexes.
- `run_vg_giraffe.sh`: Executes the VG Giraffe mapping on paired-end reads using Giraffe indexes.
- `run_GATK-HC.sh`: Runs the GATK HaplotypeCaller on the BAM outputs from linear (BWA) and VG Giraffe mappings.

### Installation Prerequisite for GATK HaplotypeCaller
- **GATK HaplotypeCaller**:
  - Installation and setup instructions can be found on the [GATK Getting Started Guide](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).
  - Ensure that Java is installed as it is a requirement for running GATK.
  - Download and configure GATK, following the detailed steps provided in the guide.

### Usage
- Make sure all prerequisite tools are installed and configured.
- Utilize `run_bwamem.sh` and `run_vg_giraffe.sh` for generating BAM files through different alignment strategies.
- Use `run_GATK-HC.sh` to perform variant calling on these BAM files using GATK HaplotypeCaller.

Refer to individual script comments for detailed usage instructions. Prerequisite indexes for these scripts can be generated using the scripts available in the `reference_construction/` module.

