# Laurel Hiatt 04/14/2025

log_dir = out_dir + "/log/3_check_bams"
bench_dir = out_dir + "/benchmark/3_check_bams"

# statistics from BAM files and outputs in a text format, summary below:
# CHK	Checksum
# SN	Summary numbers
# FFQ	First fragment qualities
# LFQ	Last fragment qualities
# GCF	GC content of first fragments
# GCL	GC content of last fragments
# GCC	ACGT content per cycle
# GCT	ACGT content per cycle, read oriented
# FBC	ACGT content per cycle for first fragments only
# FTC	ACGT raw counters for first fragments
# LBC	ACGT content per cycle for last fragments only
# LTC	ACGT raw counters for last fragments
# BCC	ACGT content per cycle for BC barcode
# CRC	ACGT content per cycle for CR barcode
# OXC	ACGT content per cycle for OX barcode
# RXC	ACGT content per cycle for RX barcode
# MPC	Mismatch distribution per cycle
# QTQ	Quality distribution for BC barcode
# CYQ	Quality distribution for CR barcode
# BZQ	Quality distribution for OX barcode
# QXQ	Quality distribution for RX barcode
# IS	Insert sizes
# RL	Read lengths
# FRL	Read lengths for first fragments only
# LRL	Read lengths for last fragments only
# MAPQ	Mapping qualities
# ID	Indel size distribution
# IC	Indels per cycle
# COV	Coverage (depth) distribution
# GCD	GC-depth
# output can be visualized graphically using plot-bamstats.

rule samtools_stats:
    input:
        bam_sort = rules.add_rg.output.sort_bam_RG
    output:
        stats = out_dir + "/bam/{sample}-sorted.stats",
    resources:
        mem_mb = mem_medium
    threads:
        4
    log:
        log_dir + "/{sample}_stats.log"
    conda:
         "../../envs/make_bams.yaml"
    shell:
        """
        samtools stats --threads {threads} {input.bam_sort} > {output.stats} 2> {log}
        """

# getting coverage across genome, chromosomes, etc
rule mosdepth:
    input:
        bam_sort = rules.add_rg.output.sort_bam_RG
    output:
        out_dir + "/mosdepth/{sample}.mosdepth.global.dist.txt"
    params:
        mos_dir = out_dir + "/mosdepth",
    resources:
        mem_mb = mem_medium
    threads:
        2
    envmodules:
        "mosdepth/0.3.1"
    log:
        log_dir + "/{sample}_mosdepth.log"
    shell:
        """
        mkdir -p {params.mos_dir}

        BASE_NAME=$(basename {input.bam_sort} -sorted.bam)

        mosdepth -n --fast-mode --by 500 \
           "${params.mos_dir}/${{BASE_NAME}}" \
           "${input.bam_sort}" 2>> {log}
        """

# plotting mosdepth results
rule plot_mosdepth:
    input:
        mosdepth = lambda wildcards: expand(
            out_dir + "/mosdepth/{sample}.mosdepth.global.dist.txt",
            sample=get_samples_for_donor(wildcards.donor, matches)
        )
    output:
        html = out_dir + "/mosdepth/{donor}_mosdepth_coverage.html"
    conda:
         "../../envs/plot_mosdepth.yaml"
    resources:
        mem_mb = mem_small
    threads:
        2
    log:
        log_dir + "/{donor}_plot_mosdepth.log"
    shell:
        """
        python ../quality_control/plot-dist.py {input.mosdepth} --output {output.html}
        """

# bam statistcs with downstream plotting of quality metrics
rule alfred_qc:
    input:
        bam = rules.add_rg.output.sort_bam_RG,  # Input is the sorted BAM from `add_rg`
        bai = rules.index_bam.output.bai,
        ref = reference
    output:
        temp(out_dir + "/alfred/{sample}.alfred.qc.json.gz")  # Per-sample Alfred QC report
    conda:
         "../../envs/alfred.yaml"
    resources:
        mem_mb = mem_small
    threads: 2
    log:
        log_dir + "/{sample}_alfred.log"
    benchmark:
        bench_dir + "/{sample}_alfred.tsv"
    shell:
        """
        alfred qc {input.bam} -r {input.ref} -o {output} > {log} 2>&1
        """

# plotting the alfred results
# This rule takes the output from the `alfred_qc` rule and generates a PDF report
rule alfred_summary:
    input:
        report = out_dir + "/alfred/{sample}.alfred.qc.json.gz"
    output:
        pdf = out_dir + "/alfred/{sample}.pdf"
    resources:
        mem_mb = mem_medium
    threads: 2
    log:
        log_dir + "/{sample}_plot_alfred.log"
    envmodules:
        "R/4.4.2"
    shell:
        """
        module load R
        Rscript quality_control/stats.R {input.report} {output.pdf}
        """