import sys
import pysam

def fix_ref_panel(input_vcf_path, output_vcf_path):
    # Open input VCF file
    in_vcf = pysam.VariantFile(input_vcf_path)
    
    # Open output VCF file for writing
    out_vcf = pysam.VariantFile(output_vcf_path, "w", header=in_vcf.header)
    
    # Process each record in the input VCF
    for rec in in_vcf:
        gt_counts = {}
        # Count genotype occurrences
        for sample in in_vcf.header.samples:
            for i in rec.samples[sample]["GT"]:
                if i is not None:
                    if i not in gt_counts:
                        gt_counts[i] = 0
                    gt_counts[i] += 1
        
        # Find the most common genotype
        if gt_counts:
            common_gt = max(gt_counts.items(), key=lambda x: x[1])[0]
        else:
            common_gt = 0
        
        # Update genotypes and phase them
        for sample in in_vcf.header.samples:
            rec.samples[sample]["GT"] = [
                i if i is not None else common_gt
                for i in rec.samples[sample]["GT"]
            ]
            rec.samples[sample].phased = True
        
        # Write the modified record to the output VCF
        out_vcf.write(rec)

    # Close the VCF files
    in_vcf.close()
    out_vcf.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_ref_panel.py <input_vcf.gz> <output_vcf.gz>")
        sys.exit(1)
    
    input_vcf_path = sys.argv[1]
    output_vcf_path = sys.argv[2]

    fix_ref_panel(input_vcf_path, output_vcf_path)

