rule somalier_extract:
    input:
        vcf = expand("{out_dir}/vcf/{donor}-var.vcf.gz", out_dir=out_dir, donor=donors)
    output:
        somalier_dir = expand("{out_dir}/somalier/{donor}/extract/", out_dir=out_dir, donor=donors),
        samples = expand("{out_dir}/somalier/{donor}/extract/{donor}_{sample}.somalier", out_dir = out_dir, donor = donors, sample = samples)
    params:
        sites= "/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/sites.hg38.vcf.gz",
        fasta= "/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"
    shell:
        """
        mkdir -p {output.somalier_dir}
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier extract {input.vcf} --sites {params.sites} --fasta {params.fasta} -d {output.somalier_dir}
        """

rule somalier_check:
    input:
        somalier_dir = expand("{out_dir}/somalier/{donor}/extract/", out_dir=out_dir, donor=donors)
    output:
        html= "{out_dir}/somalier/{donor}/relate.html",
        pairs= "{out_dir}/somalier/{donor}/relate.pairs.tsv",
        groups= "{out_dir}/somalier/{donor}/relate.groups.tsv",
        samples= "{out_dir}/somalier/{donor}/relate.samples.tsv"
    params:
        donor= "{donor}"
    shell:
        """
        /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier relate {input.somalier_dir}/*.somalier -o {out_dir}/somalier/{params.donor}/relate
        """

rule bcftools_stats:
    input:
        vcf= "{out_dir}/vcf/{donor}-annotated-var.vcf.gz"
    output:
        stats= "{out_dir}/vcf/{donor}-vcf_stats.txt"
    shell:
        "bcftools stats -s - --verbose {input.vcf} > {output.stats}"

rule plot_stats:
    input:
        stats= "{out_dir}/vcf/{donor}-vcf_stats.txt"
    output:
         outdir= "{out_dir}/vcf/{donor}",
         summaries = "{out_dir}/vcf/{donor}/summary.pdf"
    conda:
        "/uufs/chpc.utah.edu/common/HIPAA/u1264408/software/pkg/miniconda3/envs/vcfstats"
    shell:
        """
        plot-vcfstats -p {output.outdir} {input.stats}
        """