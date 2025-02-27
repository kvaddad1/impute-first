# Main rule that defines the final output
rule beagle_download:
    input:
    output:
        "beagle.18May20.d20.jar"
    shell:
        "wget http://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar"

rule imputed_calls:
    input:
        expand("imputed_final/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.beagle.vcf.gz",
               sample_name=config['sample_name'],
               cov=config['covs'],
               aligner=config['aligner'],
               align_mode=['paired', 'unpaired'],
               bcftools_mode=['filter', 'unfilter'])
rule split_ref_panel:
    input:
        ref_panel=config["ref_panel"]
    output:
        vcf="chr_vcfs/ref_panel.chr{chr}.vcf.gz",
        index="chr_vcfs/ref_panel.chr{chr}.vcf.gz.tbi"
    log:
        "logs/split_ref/chr{chr}.split.log"
    shell:
        """
        set -x
        mkdir -p chr_vcfs
        bcftools index -t {input.ref_panel}
        bcftools view {input.ref_panel} --regions chr{wildcards.chr} -Oz -o {output.vcf} 2>> {log}
        bcftools index -t {output.vcf}
        """

rule split_by_chr:
    input:
        vcf="genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.calls.vcf.gz"
    output:
        vcf="chr_vcfs/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.vcf.gz",
        index="chr_vcfs/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.vcf.gz.tbi"
    log:
        "logs/split/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.split.log"
    benchmark:
        "benchmarks/split/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.split.benchmark"
    shell:
        """
        set -x
        mkdir -p chr_vcfs
        bcftools view {input.vcf} --regions chr{wildcards.chr} -Oz -o {output.vcf} 2>> {log}
        bcftools index -t {output.vcf}
        """

rule beagle_impute:
    input:
        vcf="chr_vcfs/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.vcf.gz",
        ref="chr_vcfs/ref_panel.chr{chr}.vcf.gz",
        map=lambda w: config["beagle_maps"][w.chr],
        #jar=config["beagle_jar"]
        jar="beagle.18May20.d20.jar"
    output:
        vcf="imputed/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.imputed.vcf.gz"
    params:
        outprefix="imputed/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.imputed"
    threads: 48
    log:
        "logs/impute/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.impute.log"
    benchmark:
        "benchmarks/impute/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.impute.benchmark"
    shell:
        """
        set -x
        mkdir -p imputed
        java -jar {input.jar} \
            gt={input.vcf} \
            ref={input.ref} \
            map={input.map} \
            out={params.outprefix} \
            nthreads={threads} 2> {log}
        """
rule index_imputed:
    input:
        vcf="imputed/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.imputed.vcf.gz"
    output:
        index="imputed/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.imputed.vcf.gz.tbi"
    log:
        "logs/index/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.chr{chr}.index.log"
    shell:
        """
        bcftools index -t {input.vcf} 2> {log}
        """

rule concat_imputed:
    input:
        vcfs=expand("imputed/{{sample_name}}.r_real.c{{cov}}.{{aligner}}.{{align_mode}}.{{bcftools_mode}}.chr{chr}.imputed.vcf.gz",
                   chr=config["beagle_maps"].keys()),
        # Ensure all input files are indexed
        indexes=expand("imputed/{{sample_name}}.r_real.c{{cov}}.{{aligner}}.{{align_mode}}.{{bcftools_mode}}.chr{chr}.imputed.vcf.gz.tbi",
                      chr=config["beagle_maps"].keys())
    output:
        vcf="imputed_final/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.beagle.vcf.gz",
        index="imputed_final/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.beagle.vcf.gz.tbi"
    threads: 48
    log:
        "logs/concat/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.concat.log"
    benchmark:
        "benchmarks/concat/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.concat.benchmark"
    shell:
        """
        set -x
        mkdir -p imputed_final
        bcftools concat \
            --threads {threads} \
            -Oz \
            -o {output.vcf} \
            {input.vcfs} 2> {log}
        bcftools index --threads {threads} -t {output.vcf}
        """
