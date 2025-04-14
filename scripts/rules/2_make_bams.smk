# Laurel Hiatt 04/10/2025

log_dir = out_dir + "/log/2_make_bams/"
bench_dir = out_dir + "/benchmark/2_make_bams/"

# Aligns cleaned fastq reads to the defined reference genome
# Mark duplicates and extract discordant and split reads from sam files
# Convert to bam (exclude unmapped reads with -F 4)

rule bwa_mem:
    input:
        r1_clean = rules.fastp.output.r1_clean,
        r2_clean = rules.fastp.output.r2_clean,
        ref = reference
    output:
        temp(out_dir + "/bam/{sample}.sam")
    threads: 16
    resources:
        mem_mb = mem_xlarge
    log:
        log_dir + "{sample}_bwa_mem.log"
    benchmark:
        bench_dir + "{sample}_bwa_mem.tsv"
    conda:
        "../../envs/make_bams.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1_clean} {input.r2_clean} > {output} 2> {log}
        """

rule samblaster:
    input:
        sam = rules.bwa_mem.output
    output:
        temp(out_dir + "/bam/{sample}.samblaster.sam")
    threads: 16
    resources:
        mem_mb = mem_xlarge
    log:
        log_dir + "{sample}_samblaster.log"
    benchmark:
        bench_dir + "{sample}_samblaster.tsv"
    conda:
        "../../envs/make_bams.yaml"
    shell:
        """
        samblaster {input.sam} | samtools view -b -@ {threads} {input.sam} > {output} 2> {log}
        """

rule samtools_sort:
    input:
        bam = rules.samblaster.output
    output:
        bam_sort = out_dir + "/bam/{sample}-sortednoRG.bam"
    threads: 16
    resources:
        mem_mb = mem_xlarge
    log:
        log_dir + "{sample}_samtools_sort.log"
    benchmark:
        bench_dir + "{sample}_samtools_sort.tsv"
    conda:
        "../../envs/make_bams.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam_sort} {input.bam} 2> {log}
        """

# read groups are necessary for joint calling so we're gonna make sure they're there
rule add_rg:
    input:
        bam_sort = rules.samtools_sort.output.bam_sort
    output:
        sort_bam_RG = out_dir + "/bam/{sample}-sorted.bam"
    params:
        donor = lambda wildcards: get_donor(wildcards.sample, matches)
    threads:
        8
    resources:
        mem_mb = mem_medium
    log:
        log_dir + "{sample}_rg.log"
    benchmark:
        bench_dir + "{sample}_rg.tsv"
    conda:
         "../../envs/make_bams.yaml"
    shell:
        """
        samtools addreplacerg -r "@RG\\tID:{params.donor}_{wildcards.sample}\\tSM:{params.donor}_{wildcards.sample}\\tLB:{params.donor}_{wildcards.sample}\\tPL:Illumina" \
        --threads {threads} \
        -o {output.sort_bam_RG} \
        {input.bam_sort} > {log} 2>&1
        """

# an index helps
rule index_bam:
    input:
        bam_sort = rules.add_rg.output.sort_bam_RG
    output:
        bai = out_dir + "/bam/{sample}-sorted.bam.bai"
    threads:
        2
    resources:
        mem_mb = mem_small
    log:
        log_dir + "{sample}_index.log"
    conda:
         "../../envs/make_bams.yaml"
    localrule: True
    shell:
        """
        samtools index --threads {threads} {input.bam_sort} > {log} 2>&1
        """
