# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/6_filter_vcfs/"
bench_dir = out_dir + "/benchmark/6_filter_vcfs/"

# decompose the vcf for downstream
rule decompose_vcfs:
    input:
        vcf = out_dir + "/vcf/{donor}-var.vcf.gz",
        fasta = reference
    output:
        clean_vcf = temp(out_dir + "/vcf/{donor}-clean-var.vcf.gz")
   resources:
        mem_mb = mem_large
    threads: 8
    log:
        log_dir + "{donor}_decompose.log"
    benchmark:
        bench_dir + "{donor}_decompose.tsv"
    shell:
        """
        module load bcftools
        bcftools norm -m - {input.vcf} -w 10000 -f {input.fasta} -O b -o {output.clean_vcf}
        """

# annotate the vcf with gnomAD
rule gnomad_VCFs:
    input:
        clean_vcf = out_dir + "/vcf/{donor}-clean-var.vcf.gz",
    output:
        annotated_vcf = out_dir + "/vcf/{donor}-annotated-var.vcf.gz",
   resources:
        mem_mb = mem_large
    threads: 8
    log:
        log_dir + "{donor}_gnomad.log"
    benchmark:
        bench_dir + "{donor}_gnomad.tsv"
    shell:
        """
        module load slivar
        slivar expr -g /scratch/ucgd/lustre/common/data/Slivar/db/gnomad.hg38.genomes.v3.fix.zip -v {input.clean_vcf} -o {output.annotated_vcf}
        """

rule remove_lcr:
    input:
        annotated_vcf = out_dir + "/vcf/{donor}-annotated-var.vcf.gz",
        lcr_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/LCR-hs38.bed.gz",
        simplerepeats_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/GRCh38.UCSC.SimpleRepeats.bed.gz"
    output:
        filtered_vcf = out_dir + "/vcf/{donor}-annotated-var-noLCR.vcf.gz"
   resources:
        mem_mb = mem_xlarge
    threads: 8
    log:
        log_dir + "{donor}_noLCR.log"
    benchmark:
        bench_dir + "{donor}_noLCR.tsv"
    shell:
        """
        module load bedtools
        bedtools intersect -header -v -a {input.annotated_vcf} -b {input.lcr_bed} | \
        bedtools intersect -header -v -a - -b {input.simplerepeats_bed} | bgzip -c > {output.filtered_vcf}
        tabix -p vcf {output.filtered_vcf}
        """