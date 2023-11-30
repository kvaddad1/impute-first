## Context - Metrics 
These scripts were used to generate the call accuracy and window accuracy metrics measured in the Personalization phase of the Impute-first Alignment workflow. 
- **Call Accuracy**: Generates the precision, recall and other details on genotyped and imputed vcf files.
- **Window Accuracy**: Evaluates the accuracy of genomic windows in the context of short-range phasing.
## Contents
- `CMakeLists.txt`: CMake configuration for building the C++ program.
- `run_analysis.sh`: Script to execute the analysis process.
- `score_matched_vcf.cpp`: C++ code for scoring VCF files according to accuracy metrics.

Build instructions as follows. 
```
        mkdir -p score_matched_vcf/build;
        cd score_matched_vcf/build;
        cmake ..;
        cmake --build
```
