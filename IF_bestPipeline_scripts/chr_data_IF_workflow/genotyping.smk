# Validate configurations upfront
if 'align_mode' not in config:
    raise ValueError("align_mode must be specified in config")
if config['align_mode'] not in ['paired', 'unpaired']:
    raise ValueError("align_mode in config must be one of: paired, unpaired")

if 'aligner' not in config:
    raise ValueError("aligner must be specified in config")
if config['aligner'] not in ['bowtie2', 'bwa']:
    raise ValueError("aligner in config must be one of: bowtie2, bwa")

if 'bcftools_mode' not in config:
    raise ValueError("bcftools_mode must be specified in config")
if config['bcftools_mode'] not in ['filter', 'unfilter']:
    raise ValueError("bcftools_mode in config must be one of: filter, unfilter")

if 'covs' not in config:
    raise ValueError("covs must be specified in config")

# Rule that defines the genotyping output
rule genotype_calls:
    input:
        "genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.{bcftools_mode}.calls.vcf.gz".format(
            sample_name=config['sample_name'],
            cov=config['covs'][0],  # Using first value since it's a list
            aligner=config['aligner'],
            align_mode=config['align_mode'],
            bcftools_mode=config['bcftools_mode'])

# Downsample reads based on desired coverage
rule downsample_reads:
    input:
        fq1=config['fq1'].format(sample_name=config['sample_name']),
        fq2=config['fq2'].format(sample_name=config['sample_name'])
    output:
        fq1="downsampled_reads/{sample_name}.r_real.c{cov}.R1.fq",
        fq2="downsampled_reads/{sample_name}.r_real.c{cov}.R2.fq"
    params:
        seed=config['seed'],
        original_cov=config['read_cov']
    benchmark:
        "benchmarks/downsample/{sample_name}.r_real.c{cov}.downsample_reads.benchmark"
    log:
        "logs/downsample/{sample_name}.r_real.c{cov}.downsample_reads.log"
    shell:
        """
        set -x 
        mkdir -p downsampled_reads
        if [ "{wildcards.cov}" = "{params.original_cov}" ] || [ "{wildcards.cov}" = "NA" ]; then
            # Create symlinks for either NA or when cov matches original coverage
            ln -sf {input.fq1} {output.fq1}
            ln -sf {input.fq2} {output.fq2}
        else
            if ! command -v bc &> /dev/null; then
                echo "bc could not be found, please install it."
                exit 1
            fi
            frac=$(echo "{wildcards.cov} / {params.original_cov}" | bc -l)
            seqtk sample -s{params.seed} {input.fq1} $frac | \
                awk '{{if(NR%4==1) print $0 "/1"; else print $0}}' > {output.fq1}
            seqtk sample -s{params.seed} {input.fq2} $frac | \
                awk '{{if(NR%4==1) print $0 "/2"; else print $0}}' > {output.fq2}
        fi
        """
# Handle read combination based on alignment mode
if config['align_mode'] == 'unpaired':
    rule combine_reads:
        input:
            r1="downsampled_reads/{sample_name}.r_real.c{cov}.R1.fq" if config['covs'] != 'NA' else config['fq1'].format(sample_name=config['sample_name']),
            r2="downsampled_reads/{sample_name}.r_real.c{cov}.R2.fq" if config['covs'] != 'NA' else config['fq2'].format(sample_name=config['sample_name'])
        output:
            combined="downsampled_reads/{sample_name}.r_real.c{cov}.combined.fq"
        log:
            "logs/combine/{sample_name}.r_real.c{cov}.combine_reads.log"
        shell:
            """
            set -x
            cat {input.r1} {input.r2} > {output.combined} 2> {log}
            """

# Handle aligner selection and configuration
if config['aligner'] == 'bowtie2':
    rule bowtie2_index:
        input:
            fa=config["fa"]
        output:
            expand("indexes/{aligner}/{basename}.{suffix}.bt2",
                  aligner="bowtie2",
                  basename=config["fa"].split('/')[-1].replace('.fna', ''),
                  suffix=["1", "2", "3", "4", "rev.1", "rev.2"])
        params:
            outprefix=lambda w, output: output[0].replace(".1.bt2", "")
        threads: 48
        benchmark:
            "benchmarks/index/bowtie2_{basename}_index.benchmark".format(
                basename=config["fa"].split('/')[-1].replace('.fna', ''))
        log:
            "logs/index/bowtie2_{basename}_index.log".format(
                basename=config["fa"].split('/')[-1].replace('.fna', ''))
        shell:
            """
            set -x
            mkdir -p indexes/bowtie2
            bowtie2-build --threads {threads} {input.fa} {params.outprefix} 2> {log}
            """

    if config['align_mode'] == 'paired':
        rule aligner_align:
            input:
                reads1="downsampled_reads/{sample_name}.r_real.c{cov}.R1.fq" if config['covs'] != 'NA' else config['fq1'].format(sample_name=config['sample_name']),
                reads2="downsampled_reads/{sample_name}.r_real.c{cov}.R2.fq" if config['covs'] != 'NA' else config['fq2'].format(sample_name=config['sample_name']),
                idx=expand("indexes/bowtie2/{basename}.{suffix}.bt2",
                          basename=config["fa"].split('/')[-1].replace('.fna', ''),
                          suffix=["1", "2", "3", "4", "rev.1", "rev.2"])
            output:
                bam="aligned/{sample_name}.r_real.c{cov}.bowtie2.paired.bam"
            params:
                fa_prefix="indexes/bowtie2/" + config["fa"].split('/')[-1].replace('.fna', '')
            threads: 48
            benchmark:
                "benchmarks/align/{sample_name}.r_real.c{cov}.bowtie2.paired.align.benchmark"
            log:
                "logs/align/{sample_name}.r_real.c{cov}.bowtie2.paired.align.log"
            shell:
                """
                set -x 
                mkdir -p aligned
                bowtie2 -p {threads} -x {params.fa_prefix} -1 {input.reads1} -2 {input.reads2} 2> {log} | \
                samtools sort -o {output.bam} -
                """
    elif config['align_mode'] == 'unpaired':
        rule aligner_align:
            input:
                reads_combined="downsampled_reads/{sample_name}.r_real.c{cov}.combined.fq",
                idx=expand("indexes/bowtie2/{basename}.{suffix}.bt2",
                          basename=config["fa"].split('/')[-1].replace('.fna', ''),
                          suffix=["1", "2", "3", "4", "rev.1", "rev.2"])
            output:
                bam="aligned/{sample_name}.r_real.c{cov}.bowtie2.unpaired.bam"
            params:
                fa_prefix="indexes/bowtie2/" + config["fa"].split('/')[-1].replace('.fna', '')
            threads: 48
            benchmark:
                "benchmarks/align/{sample_name}.r_real.c{cov}.bowtie2.unpaired.align.benchmark"
            log:
                "logs/align/{sample_name}.r_real.c{cov}.bowtie2.unpaired.align.log"
            shell:
                """
                set -x 
                mkdir -p aligned
                bowtie2 -p {threads} -x {params.fa_prefix} -U {input.reads_combined} 2> {log} | \
                samtools sort -o {output.bam} -
                """

elif config['aligner'] == 'bwa':
    rule bwa_mem_index:
        input:
            fa=config["fa"]
        output:
            expand("indexes/{aligner}/{basename}.{suffix}",
                  aligner="bwa",
                  basename=config["fa"].split('/')[-1].replace('.fna', ''),
                  suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
        params:
            outprefix=lambda w, output: output[0].replace(".amb", "")
        threads: 48
        benchmark:
            "benchmarks/index/bwa_{basename}_index.benchmark".format(
                basename=config["fa"].split('/')[-1].replace('.fna', ''))
        log:
            "logs/index/bwa_{basename}_index.log".format(
                basename=config["fa"].split('/')[-1].replace('.fna', ''))
        shell:
            """
            set -x
            mkdir -p indexes/bwa
            bwa index -p {params.outprefix} {input.fa} 2> {log}
            """

    if config['align_mode'] == 'paired':
        rule aligner_align:
            input:
                reads1="downsampled_reads/{sample_name}.r_real.c{cov}.R1.fq" if config['covs'] != 'NA' else config['fq1'].format(sample_name=config['sample_name']),
                reads2="downsampled_reads/{sample_name}.r_real.c{cov}.R2.fq" if config['covs'] != 'NA' else config['fq2'].format(sample_name=config['sample_name']),
                idx=expand("indexes/bwa/{basename}.{suffix}",
                          basename=config["fa"].split('/')[-1].replace('.fna', ''),
                          suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
            output:
                bam="aligned/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".paired.bam"
            params:
                fa_prefix="indexes/bwa/" + config["fa"].split('/')[-1].replace('.fna', ''),
                L_flag=config.get('bwa_L_flag', 'default')
            threads: 48
            benchmark:
                "benchmarks/align/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".paired.benchmark"
            log:
                "logs/align/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".paired.log"
            shell:
                """
                set -x
                mkdir -p aligned
                if [ "{params.L_flag}" != "default" ]; then
                    bwa mem -t {threads} -L {params.L_flag} {params.fa_prefix} {input.reads1} {input.reads2} 2> {log} | samtools sort -o {output.bam} -
                else
                    bwa mem -t {threads} {params.fa_prefix} {input.reads1} {input.reads2} 2> {log} | samtools sort -o {output.bam} -
                fi
                """
    elif config['align_mode'] == 'unpaired':
        rule aligner_align:
            input:
                reads_combined="downsampled_reads/{sample_name}.r_real.c{cov}.combined.fq",
                idx=expand("indexes/bwa/{basename}.{suffix}",
                          basename=config["fa"].split('/')[-1].replace('.fna', ''),
                          suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
            output:
                bam="aligned/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".unpaired.bam"
            params:
                fa_prefix="indexes/bwa/" + config["fa"].split('/')[-1].replace('.fna', ''),
                L_flag=config.get('bwa_L_flag', 'default')
            threads: 48
            benchmark:
                "benchmarks/align/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".unpaired.benchmark"
            log:
                "logs/align/{sample_name}.r_real.c{cov}.bwa" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" else "_default") + ".unpaired.log"
            shell:
                """
                set -x
                mkdir -p aligned
                if [ "{params.L_flag}" != "default" ]; then
                    bwa mem -t {threads} -L {params.L_flag} {params.fa_prefix} {input.reads_combined} 2> {log} | samtools sort -o {output.bam} -
                else
                    bwa mem -t {threads} {params.fa_prefix} {input.reads_combined} 2> {log} | samtools sort -o {output.bam} -
                fi
                """

else:
    raise ValueError("Invalid aligner in config. Must be one of: bowtie2, bwa")

# Handle bcftools mode selection
if config['bcftools_mode'] == 'filter':
    rule bcftools_call:
        input:
            fa=config["fa"],
            bam=lambda wildcards: "aligned/{sample_name}.r_real.c{cov}.{aligner}".format(**wildcards) + 
                 ("_L" + config['bwa_L_flag'] if wildcards.aligner == 'bwa' and config.get('bwa_L_flag') != "default" else "_default" if wildcards.aligner == 'bwa' else "") + 
                 ".{align_mode}.bam".format(**wildcards)
        output:
            vcf="genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" and "{aligner}" == "bwa" else "_default" if "{aligner}" == "bwa" else "") + ".{align_mode}.filter.calls.vcf.gz",
            csi="genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" and "{aligner}" == "bwa" else "_default" if "{aligner}" == "bwa" else "") + ".{align_mode}.filter.calls.vcf.gz.csi"
        params:
            sample=config['sample_name'],
            bcftools_filter_QUAL=config.get('bcftools_filter_QUAL', 20),
            bcftools_filter_DP=config.get('bcftools_filter_DP', 100)
        threads: 48
        benchmark:
            "benchmarks/calls/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.filter.calls.benchmark"
        log:
            "logs/calls/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.filter.calls.log"
        shell:
            """
            set -x 
            mkdir -p genotyped_vcf
            bcftools mpileup --threads {threads} -Ou -f {input.fa} {input.bam} 2> {log} | \
            bcftools call --threads {threads} -Ou -mv | \
            bcftools filter -s LowQual -e 'QUAL<{params.bcftools_filter_QUAL} || DP>{params.bcftools_filter_DP}' | \
            bcftools view --threads {threads} -Oz -o {output.vcf} && \
            bcftools index --threads {threads} {output.vcf} && \
            bcftools index --threads {threads} -t {output.vcf}
            """
elif config['bcftools_mode'] == 'unfilter':
    rule bcftools_call:
        input:
                fa=config["fa"],
                bam=lambda wildcards: "aligned/{sample_name}.r_real.c{cov}.{aligner}".format(**wildcards) + 
                     ("_L" + config['bwa_L_flag'] if wildcards.aligner == 'bwa' and config.get('bwa_L_flag') != "default" else "_default" if wildcards.aligner == 'bwa' else "") + 
                     ".{align_mode}.bam".format(**wildcards)
        output:
            vcf="genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" and "{aligner}" == "bwa" else "_default" if "{aligner}" == "bwa" else "") + ".{align_mode}.unfilter.calls.vcf.gz",
            csi="genotyped_vcf/{sample_name}.r_real.c{cov}.{aligner}" + ("_L" + config['bwa_L_flag'] if config.get('bwa_L_flag') != "default" and "{aligner}" == "bwa" else "_default" if "{aligner}" == "bwa" else "") + ".{align_mode}.unfilter.calls.vcf.gz.csi"
        params:
            sample=config['sample_name']
        threads: 48
        benchmark:
            "benchmarks/calls/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.unfilter.calls.benchmark"
        log:
            "logs/calls/{sample_name}.r_real.c{cov}.{aligner}.{align_mode}.unfilter.calls.log"
        shell:
            """
            set -x 
            mkdir -p genotyped_vcf
            bcftools mpileup --threads {threads} -Ou -f {input.fa} {input.bam} 2> {log} | \
            bcftools call --threads {threads} -Ou -mv | \
            bcftools view --threads {threads} -Oz -o {output.vcf} && \
            bcftools index --threads {threads} {output.vcf} && \
            bcftools index --threads {threads} -t {output.vcf}
            """
else:
    raise ValueError("Invalid bcftools_mode in config. Must be one of: filter, unfilter")
