## Overview
This README details the evaluation metrics of the Impute-First Alignment Workflow. 

### Upstream Analysis (Personalization phase)
- **Metrics**: Genotyping and Imputation Calls evaluation: Call Accuracy, Window Accuracy.

### Downstream Analysis
- **Metrics**: Allelic Balance at HETs, Variant Call Accuracy.

## Installation and Usage

1. **Biastools**:
   - Main scripts provided in `./downstream_eval/biastools`.
   - Additional usage information is available in the associated script directory or in the [Biastools GitHub Repository](https://github.com/maojanlin/biastools).

2. **RTG Tools**:
   - Available at: [RTG Tools 3.12.1](https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12.1/rtg-tools-3.12.1-linux-x64.zip).

3. **hap.py (via Singularity)**:
   - The hap.py image (v0.3.15) can be pulled using:
     ```bash
     singularity pull docker://pkrusche/hap.py:v0.3.15
     ```
   - This generates a Singularity image file named `hap.py_v0.3.15.sif` in the working directory.
   - When using `--engine=vcfeval`, the RTG tools path must be specified with:
     ```bash
     --engine-vcfeval-path /path/to/rtg-tools
     ```

Additional instructions and parameters are documented within each respective script directory.

