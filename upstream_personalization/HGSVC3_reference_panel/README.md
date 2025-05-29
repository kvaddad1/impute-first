## HGSVC3 Panel Preparation

The HGSVC3 reference panel yielded the highest accuracy in our framework evaluations, and this directory provides the corresponding VCF files used in those analyses.

### Included Script

- `prepare_hgsvc3_panel.sh`: This script automates the complete preparation of the HGSVC3 panel. It:
  - Downloads SNV, indel, and SV VCFs from the official HGSVC3 GRCh38 release
  - Indexes them using `bcftools`
  - Concatenates and normalizes the files into a single unified panel
  - Fills allele count (AC) and number (AN) tags
  - Applies a custom fix script (`fix_VCF.py`)
  - Generates a version with HG002 (NA24385) excluded, for trio-based imputation

### Notes on Panel Use

- `hsvc3.all.fixed.vcf.gz`: **Full HGSVC3 panel** including all samples.  
  Used for imputing HG001 and HG005 (not present in HGSVC3).

- `hsvc3.all.less1.fixed.vcf.gz`: **HG002-excluded panel**, created by removing sample NA24385.  
  Used for imputing **HG002, HG003, and HG004**, which form an Ashkenazim trio.

### Files Generated

| File Name                        | Description                                  |
|----------------------------------|----------------------------------------------|
| `hsvc.all.vcf.gz`               | Concatenated and normalized raw panel        |
| `hsvc3.all.vcf.gz`              | Tag-filled version (AC, AN)                  |
| `hsvc3.all.fixed.vcf.gz`        | Post-processed version after `fix_VCF.py`    |
| `hsvc3.all.less1.fixed.vcf.gz`  | HG002-excluded version                       |

Make sure `fix_VCF.py` is available in the directory and properly handles your custom field corrections.

