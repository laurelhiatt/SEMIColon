# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/2_make_bams/"
bench_dir = out_dir + "/benchmark/2_make_bams/"

# Aligns cleaned fastq reads to the defined reference genome
# Mark duplicates and extract discordant and split reads from sam files
# Convert to bam (exclude unmapped reads with -F 4)

rule align_and_sort:
    input:
        r1_clean = rules.fastp.output.r1_clean,
        r2_clean = rules.fastp.output.r2_clean,
        ref = reference
    output:
        bam_sort = temp(out_dir + "/bam/{sample}-sortednoRG.bam")
    threads: 16
    resources:
        mem_mb = mem_large
    log:
        log_dir + "{sample}_align_sort.log"
    benchmark:
        bench_dir + "{sample}_align_sort.tsv"
    conda:
         "../../envs/make_bams.yaml"
    shell:
        """
        module load bwa
        module load samblaster
        module load samtools
        bwa mem -t {threads} {input.ref} {input.r1_clean} {input.r2_clean} | samblaster | samtools view -b | samtools sort -o {output.bam_sort}
        """

# read groups are necessary for joint calling so we're gonna make sure they're there
rule add_rg:
    input:
        bam_sort = rules.align_and_sort.output.bam_sort
    output:
        sort_bam_RG = out_dir + "/bam/{sample}-sorted.bam"
    params:
        donor = lambda wildcards: get_donor(wildcards.sample, matches)
    threads: 8
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
        samtools addreplacerg -r "@RG\\tID:{params.donor}_{wildcards.sample}\\tSM:{params.donor}_{wildcards.sample}\\tLB:{params.donor}_{wildcards.sample}\\tPL:Illumina" {input.bam_sort} | \
        samtools view -b > {output.sort_bam_RG}
        """

# an index helps
rule index_bam:
    input:
        bam_sort = rules.add_rg.output.sort_bam_RG
    output:
        bai = out_dir + "/bam/{sample}-sorted.bam.bai"
    threads: 4
    resources:
        mem_mb = mem_large
    log:
        log_dir + "{sample}_index.log"
    benchmark:
        bench_dir + "{sample}_index.tsv"
    conda:
         "../../envs/make_bams.yaml"
    shell:
        """
        samtools index {input.bam_sort}
        """
