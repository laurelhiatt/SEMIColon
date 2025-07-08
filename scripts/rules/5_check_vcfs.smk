# Laurel Hiatt 04/14/2025

log_dir = out_dir + "/log/5_check_vcfs"
bench_dir = out_dir + "/benchmark/5_check_vcfs"

# check relatedness
# check relatedness
rule somalier_extract:
    input:
        vcf = out_dir + "/vcf/{donor}-var.vcf.gz"
    output:
        somalier_dir = directory(out_dir + "/somalier/{donor}/extract/"),
    params:
        sites= "/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/sites.hg38.vcf.gz",
        fasta = reference
    resources:
        mem_mb = mem_medium
    log:
        log_dir + "/{donor}_somalier.log"
    benchmark:
        bench_dir + "/{donor}_somalier.tsv"
    threads:
        1
    shell:
        """
        mkdir -p {output.somalier_dir}
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier extract {input.vcf} --sites {params.sites} --fasta {params.fasta} -d {output.somalier_dir} > {log} 2>&1
        """

#finish somalier analysis
rule somalier_check:
    input:
        somalier_dir = rules.somalier_extract.output.somalier_dir
    output:
        html = out_dir + "/somalier/{donor}/relate.html",
        pairs = out_dir + "/somalier/{donor}/relate.pairs.tsv",
        groups = out_dir + "/somalier/{donor}/relate.groups.tsv",
        samples = out_dir + "/somalier/{donor}/relate.samples.tsv"
    resources:
        mem_mb = mem_medium
    threads: 1
    log:
        log_dir + "/{donor}_somalier_check.log"
    benchmark:
        bench_dir + "/{donor}_somalier_check.tsv"
    shell:
        """
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier relate {input.somalier_dir}/*.somalier -o {out_dir}/somalier/{wildcards.donor}/relate
        """

#code for the bcftools stats
rule bcftools_stats:
    input:
        vcf = out_dir + "/vcf/{donor}-annotated-var.vcf.gz"
    output:
        stats = temp(out_dir + "/vcf/{donor}-vcf_stats.txt")
    resources:
        mem_mb = mem_small
    threads: 4
    log:
        log_dir + "/{donor}_bcftools_stats.log"
    envmodules:
        "bcftools/1.16"
    shell:
        """
        module load bcftools/1.16
        bcftools stats -s - --verbose --threads {threads} {input.vcf} > {output.stats} 2> {log}
        """

# plot the bcftools stats output
rule plot_stats:
    input:
        stats = out_dir + "/vcf/{donor}-vcf_stats.txt"
    output:
         outdir = out_dir + "/vcf/{donor}",
         summary = out_dir + "/vcf/{donor}/summary.pdf"
    conda:
        "../../envs/vcfstats.yaml"
    resources:
        mem_mb = mem_medium
    threads: 2
    log:
        log_dir + "/{donor}_bcftools_plot.log"
    benchmark:
        bench_dir + "/{donor}_bcftools_plot.tsv"
    shell:
        """
        mkdir -p {output.outdir}
        tmpdir=$(mktemp -d)

        plot-vcfstats -p $tmpdir {input.stats}

        mv ${{tmpdir}}/summary.pdf {output.summary}
        rm -r $tmpdir
        """