# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/1_fastp/"
bench_dir = out_dir + "/benchmark/1_fastp/"

# Rule 'fastp' preprocesses the reads and creates
# filtered fastq reads with preprocessing reports.
rule fastp:
    input:
        R1 = in_dir + "/merged/{sample}_R1.fastq.gz",
        R2 = in_dir + "/merged/{sample}_R2.fastq.gz"
    output:
        r1_clean = out_dir + "/fastq/{sample}_R1.clean.fastq.gz",
        r2_clean = out_dir + "/fastq/{sample}_R2.clean.fastq.gz",
        html_report = out_dir + "/reports/{sample}-fastp-report.html",
        json_report = out_dir + "/reports/{sample}-fastp-report.json"
    resources:
        mem_mb = mem_medium
    threads:
        8
    log:
        log_dir + "{sample}.log"
    benchmark:
        bench_dir + "{sample}.tsv"
    conda:
        "../../envs/fastp.yaml"
    shell:
        """
        fastp --in1 {input.R1} --in2 {input.R2} \
        --out1 {output.r1_clean} --out2 {output.r2_clean} \
        --disable_quality_filtering --disable_adapter_trimming --disable_trim_poly_g \
        --html {output.html_report} \
        --json {output.json_report} \
        --thread {threads}
        """
