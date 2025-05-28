# glimpse1_imputation.smk

# Install GLIMPSE binaries (keep existing)
rule install_glimpse:
    output:
        "glimpse/GLIMPSE_chunk_static",
        "glimpse/GLIMPSE_phase_static", 
        "glimpse/GLIMPSE_ligate_static",
        "glimpse/GLIMPSE_sample_static"
    shell:
        """ 
        mkdir -p glimpse
        wget -O glimpse/GLIMPSE_chunk_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        wget -O glimpse/GLIMPSE_phase_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        wget -O glimpse/GLIMPSE_ligate_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_ligate_static
        wget -O glimpse/GLIMPSE_sample_static https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_sample_static
        chmod +x glimpse/GLIMPSE_*_static
        """

# Existing split_ref_panel and split_by_chr rules remain the same

rule glimpse_chunk:
    input:
        chunk="glimpse/GLIMPSE_chunk_static",
        vcf="chr_vcfs/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.vcf.gz",
        ref="chr_vcfs/ref_panel.chr{chr}.vcf.gz"
    output:
        chunks="chunks/{sample_name}.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.txt"
    params:
        region=lambda w: f"chr{w.chr}"
    threads: 32
    shell:
        """
        mkdir -p chunks
        {input.chunk} \
            --input {input.vcf} \
            --thread {threads} \
            --region {params.region} \
            --window-size 2000000 \
            --buffer-size 200000 \
            --output {output.chunks}
        """

rule glimpse_phase:
    input:
        phase="glimpse/GLIMPSE_phase_static",
        vcf="chr_vcfs/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.vcf.gz",
        ref="chr_vcfs/ref_panel.chr{chr}.vcf.gz",
        map=lambda w: config["glimpse_maps"][w.chr],
        chunks="chunks/{sample_name}.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.txt"
    output:
        phased=directory("phase/{sample_name}.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}")
    threads: 32
    shell:
        """
        mkdir -p {output.phased}
        while IFS=$'\t' read -r chunk_num chrom input_reg output_reg window_size variant_count; do
            # Pad chunk number to 2 digits
            chunk_id=$(printf "%02d" $chunk_num)
            
            {input.phase} \
                --thread {threads} \
                --input {input.vcf} \
                --reference {input.ref} \
                --map {input.map} \
                --input-region "$input_reg" \
                --output-region "$output_reg" \
                --output {output.phased}/phased_$chunk_id.bcf
                
            bcftools index -f {output.phased}/phased_$chunk_id.bcf
        done < <(tail -n +2 {input.chunks})  # Skip header line
        """

rule glimpse_ligate:
    input:
        ligate="glimpse/GLIMPSE_ligate_static",
        phased=directory("phase/{sample_name}.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}")
    output:
        vcf="imputed/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.glimpse1.vcf.gz"
    threads: 32
    shell:
        """
        ls {input.phased}/phased_*.bcf | sort -V > filelist.txt
        {input.ligate} \
            --thread {threads} \
            --input filelist.txt \
            --output {output.vcf}
        bcftools index -t {output.vcf}
        """
rule concat_imputed:
    input:
        expand("imputed/{{sample_name}}.r_real.c{{cov}}.{{aligner}}.{{align_mode}}.{{bcftools_mode}}.chr{chr}.glimpse1.vcf.gz",
               chr=config["glimpse_maps"].keys())
    output:
        vcf="imputed_final/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.glimpse1.vcf.gz"
    threads: 48
    shell:
        """
        bcftools concat --threads {threads} -Oz -o {output.vcf} {input}
        bcftools index -t {output.vcf}
        """
