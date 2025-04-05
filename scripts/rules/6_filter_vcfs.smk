rule decompose_vcfs:
    input:
        vcf = "{out_dir}/vcf/{donor}-var.vcf.gz",
        fasta = reference
    output:
        clean_vcf = temp("{out_dir}/vcf/{donor}-clean-var.vcf.gz"),
    shell:
        """
        bcftools norm -m - {input.vcf} -w 10000 -f {input.fasta} -O b -o {output.clean_vcf}
        """

rule gnomad_VCFs:
    input:
        clean_vcf = "{out_dir}/vcf/{donor}-clean-var.vcf.gz",
    output:
        annotated_vcf = "{out_dir}/vcf/{donor}-annotated-var.vcf.gz",
    shell:
        """
        slivar expr -g /scratch/ucgd/lustre/common/data/Slivar/db/gnomad.hg38.genomes.v3.fix.zip -v {input.clean_vcf} -o {output.annotated_vcf}
        """

rule remove_lcr:
    input:
        annotated_vcf = "{out_dir}/vcf/{donor}-annotated-var.vcf.gz",
        lcr_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/LCR-hs38.bed.gz",
        simplerepeats_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/GRCh38.UCSC.SimpleRepeats.bed.gz"
    output:
        filtered_vcf = "{out_dir}/vcf/{donor}-annotated-var-noLCR.vcf.gz"
    shell:
        """
        bedtools intersect -header -v -a {input.annotated_vcf} -b {input.lcr_bed} | \
        bedtools intersect -header -v -a - -b {input.simplerepeats_bed} | bgzip -c > {output.filtered_vcf}
        tabix -p vcf {output.filtered_vcf}
        """