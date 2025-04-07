# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/0_merge_fastq/"
bench_dir = out_dir + "/benchmark/0_merge_fastq/"

# Snakemake rule for merging initial input fastq.gz files to one R1 and one R2
# necessary for top off sequencing or samples split across flow cells
rule merge_fastq:
    input:
        R1= lambda wildcards: sample_files.get(wildcards.sample, {}).get('R1', []),
        R2= lambda wildcards: sample_files.get(wildcards.sample, {}).get('R2', [])
    output:
        R1_merged = in_dir + "/merged/{sample}_R1.fastq.gz",
        R2_merged = in_dir + "/merged/{sample}_R2.fastq.gz"
    shell:
        """
        zcat {input.R1} | gzip > {output.R1_merged}
        zcat {input.R2} | gzip > {output.R2_merged}
        """
    resources:
        mem_mb = mem_small
    threads:
        8
    log:
        log_dir + "{sample}.log"
    benchmark:
        bench_dir + "{sample}.tsv"



