# Laurel Hiatt 04/14/2025
log_dir = out_dir + "/log/6_filter_vcfs"
bench_dir = out_dir + "/benchmark/6_filter_vcfs"

pairs = [
    (get_donor(sample, matches), sample)
    for sample in samples
    if get_donor(sample, matches) != "gibberish"
]

if pairs:
    donors, filtered_samples = zip(*pairs)
else:
    donors, filtered_samples = [], []

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
    envmodules:
        "bcftools/1.16"
    log:
        log_dir + "/{donor}_decompose.log"
    benchmark:
        bench_dir + "/{donor}_decompose.tsv"
    shell:
        """
        module load bcftools/1.16
        bcftools norm -m - {input.vcf} --threads {threads} -w 10000 -f {input.fasta} -O b -o {output.clean_vcf} > {log} 2>&1
        """

# annotate the vcf with gnomAD
rule gnomad_VCFs:
    input:
        clean_vcf = out_dir + "/vcf/{donor}-clean-var.vcf.gz",
    output:
        annotated_vcf = temp(out_dir + "/vcf/{donor}-annotated-var.vcf.gz"),
    resources:
        mem_mb = mem_large
    threads: 2
    log:
        log_dir + "/{donor}_gnomad.log"
    envmodules:
        "slivar/0.3.1"
    shell:
        """
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
    threads: 4
    log:
        log_dir + "/{donor}_noLCR.log"
    benchmark:
        bench_dir + "/{donor}_noLCR.tsv"
    shell:
        """
        module load bedtools
        bedtools intersect -header -v -a {input.annotated_vcf} -b {input.lcr_bed} | \
        bedtools intersect -header -v -a - -b {input.simplerepeats_bed} | bgzip -c > {output.filtered_vcf}
        tabix -p vcf {output.filtered_vcf}
        """



rule filter_by_depth:
    input:
        filtered_vcf= out_dir + "/vcf/{donor}-annotated-var-noLCR.vcf.gz"
    output:
        depth_vcf= temp(out_dir + "/results/{donor}-depth-filtered.vcf.gz")
    shell:
        """
        ./vcfexpress filter -e 'return all(function (dp) return dp > 5 end, variant:format("DP"))' -o {output.depth_vcf} {input.filtered_vcf}
        """

rule index_depth:
    input:
        depth_vcf= out_dir + "/results/{donor}-depth-filtered.vcf.gz"
    output:
        indexed_vcf= temp(out_dir + "/results/{donor}-depth-filtered.vcf.gz.tbi")
    threads:
        8
    shell:
        """
        tabix -f -p vcf {input.depth_vcf} --threads {threads}
        """

### it shouldn't exist in gnomad (at this point, pre-blood)
rule filter_by_gnomad:
    input:
        depth_vcf= out_dir + "/results/{donor}-depth-filtered.vcf.gz",
        indexed_vcf= out_dir + "/results/{donor}-depth-filtered.vcf.gz.tbi"
    output:
        gnomad_vcf= out_dir + "/results/{donor}-gnomad-filtered.vcf.gz",
        done = out_dir + "/results/{donor}/gnomad.done"
    shell:
        """
        module load bcftools
        bcftools filter -i '(gnomad_popmax_af <= 0 || gnomad_popmax_af == ".")' -Oz -o {output.gnomad_vcf} {input.depth_vcf}
        touch {output.done}
        """

def donor_done_input(wildcards):
    donor = wildcards.donor
    sample = wildcards.sample
    crypt_samples = matches.get(donor, {}).get("crypt_samples", [])
    if sample not in crypt_samples:
        pass
    return os.path.join(out_dir, "results", donor, "gnomad.done")

rule make_lua:
    input:
        gnomad_done = donor_done_input,
    output:
        lua = out_dir + "/results/{donor}/{sample}_soi.lua"
    run:
        with open(output.lua, "w") as f:
            f.write(f'samplesOfInterest = {{"{wildcards.donor}_{wildcards.sample}"}}\n')

### all the other samples must have 0 alt alleles
rule filter_by_sample:
    input:
        lua = rules.make_lua.output.lua,
    output:
        sample_vcf = temp(out_dir + "/results/{donor}/{sample}_filtered_noAD.vcf.gz")
    params:
        gnomad_vcf = out_dir + "/results/{donor}-gnomad-filtered.vcf.gz"
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua \
            -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 0 end, sampleIndexes, variant:format("AD"))' \
            -o {output.sample_vcf} {params.gnomad_vcf}
        """

### this sample must have > # alternate allele
rule filter_by_alt_depth:
    input:
        sample_vcf= out_dir + "/results/{donor}/{sample}_filtered_noAD.vcf.gz",
        lua = rules.make_lua.output.lua,
    output:
        vcf= out_dir + "/results/{donor}/{sample}_filtered.vcf.gz",
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 2 end, sampleIndexes, variant:format("AD"))' -o {output.vcf} {input.sample_vcf}
        """

rule count_snvs:
    input:
        vcf= out_dir + "/results/{donor}/{sample}_filtered.vcf.gz",
    output:
        out_vcf= out_dir + "/results/{donor}/{sample}_filtered_snvs.vcf.gz",
        txt= out_dir + "/results/{donor}/{sample}_snv_count.txt"
    params:
        sample = "{sample}"
    conda:
        "../../envs/cyvcf2.yaml"
    script:
        "../mutations/per-sample-report.py"
