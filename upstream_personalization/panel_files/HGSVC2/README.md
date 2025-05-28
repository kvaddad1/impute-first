# Genetic Map and VCF Fixing Scripts

This folder contains a set of scripts designed for downloading genetic map files, panel files, and fixing VCF files for HGSVC2 release files, to be used in imputation. 

## Contents
- `download_linkage_maps.sh`: Script to download genetic map files.
- `download_panel_files.sh`: Script to download panel files (e.g., HGSVC2).
- `fill_tags.sh`: (Optional) Script to fill in `AC` and `AN` info if missing.
- `fix_VCF.py`: Python script to fix missing `GT` or unphased genotype info.
- `fix_vcf.sh`: Shell script to run `fix_VCF.py` for fixing VCF issues.
- `wrap_run_fixing_vcf.sh`: Wrapper script to automate the process of fixing VCF files and preparing a Leave-One-Out (LOO) analysis panel.

## Usage
### Fixing the VCF and Preparing the LOO Panel
To fix VCF issues and prepare the LOO panel, run the wrapper script:

```bash
./wrap_run_fixing_vcf.sh
```
The `HGSVC2_LOO` panel files for `HG001/HG002` will be at `ref_panel/hsvc.all.less1.HG001.fixed.vcf.gz`
and `ref_panel/hsvc.all.less1.HG002.fixed.vcf.gz`

