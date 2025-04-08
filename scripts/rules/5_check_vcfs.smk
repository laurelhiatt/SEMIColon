# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/5_check_vcfs/"
bench_dir = out_dir + "/benchmark/5_check_vcfs/"

# check relatedness
rule somalier_extract:
    input:
        vcf = out_dir + "/vcf/{donor}-var.vcf.gz"
    output:
        somalier = out_dir + "/somalier/{donor}/extract/{donor}_{sample}.somalier"
    params:
        sites= "/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/sites.hg38.vcf.gz",
        somalier_dir = out_dir + "/somalier/{donor}/extract/",
        fasta = reference
    resources:
        mem_mb = mem_medium
    threads: 4
    log:
        log_dir + "{donor}_{sample}_somalier.log"
    benchmark:
        bench_dir + "{donor}_{sample}_somalier.tsv"
    shell:
        """
        mkdir -p {params.somalier_dir}
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier extract {input.vcf} --sites {params.sites} --fasta {params.fasta} -d {params.somalier_dir}
        """

#finish somalier analysis
rule somalier_check:
    input:
        somalier_files = lambda wildcards: expand(
            out_dir + "/somalier/{donor}/extract/{donor}_{sample}.somalier",
            donor=wildcards.donor,
            sample=matches[wildcards.donor]
        )
    output:
        html = out_dir + "/somalier/{donor}/relate.html",
        pairs = out_dir + "/somalier/{donor}/relate.pairs.tsv",
        groups = out_dir + "/somalier/{donor}/relate.groups.tsv",
        samples = out_dir + "/somalier/{donor}/relate.samples.tsv"
    params:
        donor = "{donor}"
    resources:
        mem_mb = mem_medium
    threads: 4
    log:
        log_dir + "{donor}_somalier_check.log"
    benchmark:
        bench_dir + "{donor}_somalier_check.tsv"
    shell:
        """
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier relate {input.somalier_dir}/*.somalier -o {out_dir}/somalier/{params.donor}/relate
        """

#code for the bcftools stats
rule bcftools_stats:
    input:
        vcf = out_dir + "/vcf/{donor}-annotated-var.vcf.gz"
    output:
        stats = out_dir + "/vcf/{donor}-vcf_stats.txt"

    resources:
        mem_mb = mem_small
    threads: 4
    log:
        log_dir + "{donor}_bcftools_stats.log"
    benchmark:
        bench_dir + "{donor}_bcftools_stat.tsv"
    shell:
        """
        module load bcftools
        bcftools stats -s - --verbose {input.vcf} > {output.stats}
        """

# plot the bcftools stats output
rule plot_stats:
    input:
        stats = out_dir + "/vcf/{donor}-vcf_stats.txt"
    output:
         outdir = out_dir + "/vcf/{donor}",
         summaries = out_dir + "/vcf/{donor}/summary.pdf"
    conda:
        "/uufs/chpc.utah.edu/common/HIPAA/u1264408/software/pkg/miniconda3/envs/vcfstats"
    resources:
        mem_mb = mem_medium
    threads: 4
    log:
        log_dir + "{donor}_bcftools_plot.log"
    benchmark:
        bench_dir + "{donor}_bcftools_plot.tsv"
    shell:
        """
        plot-vcfstats -p {output.outdir} {input.stats}
        """