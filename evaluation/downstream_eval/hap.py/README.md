### Metrics Overview

Variant calls are benchmarked against a truth set using `hap.py` with the `vcfeval` engine.  
Evaluations are restricted to confident regions defined by a BED file. Optional regional stratification is also supported.  
The reference FASTA must have a corresponding SDF-formatted version required by `vcfeval`.

More details: https://github.com/Illumina/hap.py/blob/master/doc/happy.md

---

### Metric Extraction

From each evaluation output, per-variant-type metrics (e.g., SNP, INDEL) are extracted.  
The following columns are used:

- `TRUTH.TP` — True positives in truth representation  
- `QUERY.TP` — True positives in query representation  
- `QUERY.FP` — False positives  
- `TRUTH.FN` — False negatives  
- `METRIC.Recall`  
- `METRIC.Precision`  

These metrics reflect performance within the confident regions for each sample and variant type.

---

### Aggregated Metrics

Accuracy statistics are then aggregated across all variant types for each sample and pipeline.  
Computed metrics include:

- **Precision** = QUERY.TP / (QUERY.TP + QUERY.FP)  
- **Recall** = TRUTH.TP / (TRUTH.TP + TRUTH.FN)  
- **F1 Score** = 2 × Precision × Recall / (Precision + Recall)

This allows high-level comparison of overall variant calling performance across different workflows.

