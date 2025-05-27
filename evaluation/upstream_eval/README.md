## Impute-first Call Accuracy and Window Accuracy Metrics  
These scripts were used to generate the `call accuracy` and `window accuracy` metrics measured in the Personalization phase of the Impute-first workflow.  
- **Call Accuracy**: Evaluates precision, recall and other metrics by comparing genotyped/imputed calls against corresponding benchmark genotyped/imputed VCF files.
- **Window Accuracy**: Measures the concordance of variants within sliding genomic windows (default 200bp) to evaluate short-range phasing accuracy.

## Contents
- `CMakeLists.txt`: CMake configuration for building the `score_matched_vcf` executable.
- `gen_upstream_metrics.sh`: Script to handle VCF merging, normalization and prepare input for metrics calculation.
- `score_matched_vcf.cpp`: C++ code to calculate the metrics.
- `upstream_metrics.yml`: Conda environment file for setting up the environment.

## Metrics Calculation

### Call Accuracy Metrics
For alternative allele calls:
- Each allele in a genotype is compared between test and truth VCF
- True Positive (TP): Alternative allele present in both test and truth genotype
- True Negative (TN): Reference allele present in both genotypes
- False Positive (FP): Alternative allele in test but not in truth genotype
- False Negative (FN): Reference allele in test but alternative in truth genotype

For heterozygous calls:
- True Heterozygous (thets): Heterozygous variant present in both test and truth
- False Heterozygous (fhets): Heterozygous in test but homozygous in truth
- False Homozygous (fhoms): Homozygous in test but heterozygous in truth
- True Homozygous (thoms): Homozygous variant present in both test and truth

Precision and recall are calculated for both metrics:
- Precision = TP / (TP + FP)
- Recall = TP / (TP + FN)

### Window Accuracy Metrics
For each sliding window of 200bp, accuracy is measured as:
- Window Accuracy = Matching Windows / Total Windows

Windows are categorized by variant density:
- Low density: 1-5 variants per window
- Medium density: 6-10 variants per window  
- High density: >10 variants per window

For example:
```
perc_matching_200bp_windows_1_5 = matching_200bp_windows_1_5 / total_200bp_windows_1_5
```

## Build Instructions
```
conda env create -f upstream_metrics.yml
conda activate IF_upstream_metrics

# Build the tool
mkdir build
cd build
cmake ../
cmake --build .
```
## Usage Options

### Using the Wrapper Script 
The `gen_upstream_metrics.sh` script handles the preprocessing steps and calls the executable:
```
./gen_upstream_metrics.sh <benchmark_vcf> <input_vcf> <window_size> <debugfile> <output_prefix>
```
Required arguments:
- benchmark_vcf: Truth/reference VCF file
- input_vcf: VCF file to evaluate
- window_size: Window size in base pairs (use `0` for `call accuracy` metrics, `200` for `window accuracy`)
- debugfile: Output file for variant comparison details (optional, use /dev/null if not needed)
```
# Each line represents one variant position comparison (0-based positions)
# <sequence_name> <chromosome_internal_id> <position-1> <truth_genotype> <call_genotype>
chr1 0 10067 1/1 1/0    # Example: chromosome1, internal_id=0, position=10068 in VCF, truth=homozygous alt, call=heterozygous
```
- output_prefix: Prefix for output JSON file

Output stats will be written to `<output_prefix>.json`

The script will:
1. Merge the benchmark and input VCFs
2. Normalize variants
3. Feed the normalized VCF to score_matched_vcf executable

## Output Format
The output JSON file contains metrics in the following format:
```
{
  "var_type": "any",          # Variant type (any/snp/indel/sv)
  "total_vars": int,          # Total variants analyzed
  "fps": int,                 # False positives
  "fns": int,                 # False negatives  
  "tps": int,                 # True positives
  "tns": int,                 # True negatives
  "fhoms": int,              # False homozygous
  "fhets": int,              # False heterozygous
  "thoms": int,              # True homozygous
  "thets": int,              # True heterozygous
  "alt_precision": float,     # Precision for alternate alleles
  "alt_recall": float,        # Recall for alternate alleles
  "het_precision": float,     # Precision for heterozygous calls
  "het_recall": float,        # Recall for heterozygous calls
  
  # Window accuracy metrics (are available for 'any' var_type when window_size is set to 200)
  "matching_200bp_windows_1_5": int,     # Windows with 1-5 variants matching
  "matching_200bp_windows_6_10": int,    # Windows with 6-10 variants matching
  "matching_200bp_windows_11_n": int,    # Windows with 11+ variants matching
  "total_200bp_windows_1_5": int,        # Total windows with 1-5 variants
  "total_200bp_windows_6_10": int,       # Total windows with 6-10 variants  
  "total_200bp_windows_11_n": int,       # Total windows with 11+ variants
  "perc_matching_200bp_windows_1_5": float,  # Percentage matching 1-5
  "perc_matching_200bp_windows_6_10": float, # Percentage matching 6-10
  "perc_matching_200bp_windows_11_n": float  # Percentage matching 11+
}
```
