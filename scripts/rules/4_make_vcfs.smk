# Laurel Hiatt 04/14/2025

log_dir = out_dir + "/log/4_make_vcfs"
bench_dir = out_dir + "/benchmark/4_make_vcfs"

rule make_bam_list:
    input:
        bam_sort = expand(out_dir + "/bam/{sample}-sorted.bam",
            sample=samples),
        bai = expand(out_dir + "/bam/{sample}-sorted.bam.bai",
            sample=samples),
        donor = "{donor}"
    params:
        out_dir = out_dir,
        matches = matches
    output:
        bam_file_list = temp(out_dir + "/bam/bam_list_{donor}.txt")
    resources:
        mem_mb = mem_xsmall
    threads: 2
    localrule: True
    log:
        log_dir + "/{donor}_bamlist.log"
    script:
        "../variant_calling/create_bam_list.py"

rule generate_regions:
    input:
        chroms = chroms
    params:
        index = reference_index,
        chunks = chunks,
        out_dir = out_dir
    output:
        regions = out_dir + "/regions/chunk.{chroms}.region.{i}.bed"
    resources:
        mem_mb = mem_small
    threads: 2
    localrule: True
    conda:
         "../../envs/plot_mosdepth.yaml"
    log:
        log_dir + "/{chroms}_{i}_generateregions.log"
    script:
        """
        ../variant_calling/fasta_generate_regions.py
        """

# Current filtering:
### min-alternate-count 2
## qsum 40
rule freebayes_variant_calling:
    input:
        bam_file_list = rules.make_bam_list.output.bam_file_list,
        regions = rules.generate_regions.output.regions,
        ref = reference,
    output:
        full = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf")
    params:
        out_dir = out_dir
    resources:
        mem_mb = mem_xlarge
    threads:
        8
    log:
        log_dir + "/{chroms}_{i}_{donor}_freebayes.log"
    benchmark:
        bench_dir + "/{chroms}_{i}_{donor}_freebayes.tsv"
    envmodules:
        "freebayes/1.3.9"
    shell:
        """
        mkdir -p {params.out_dir}/vcf/{wildcards.chroms}
        freebayes --min-alternate-count 2 --min-alternate-qsum 40 -f {input.ref} -t {input.regions} -L {input.bam_file_list} > {output.full} 2> {log}
        """

rule compress_chunks:
    input:
        chunk_vcf = out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf"
    output:
        chunk_zip = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz")
    resources:
        mem_mb = mem_medium
    threads: 4
    shell:
        """
        bgzip -f {input.chunk_vcf} --threads {threads}
        """

rule index_chunks:
    input:
        chunk_zip = out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz"
    output:
        chunk_zip_index = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi")
    threads: 4
    resources:
        mem_mb = mem_medium
    shell:
        """
        tabix -f -p vcf {input.chunk_zip} --threads {threads}
        """

rule concat_vcfs:
    input:
        chunk_zip = lambda wildcard: expand(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz",
            i = chunks,
            donor = wildcard.donor,
            chroms = chroms),
        chunk_vcf_index = lambda wildcard: expand(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi",
            i = chunks,
            donor = wildcard.donor,
            chroms = chroms),
    output:
        vcf = out_dir + "/vcf/{donor}-var.vcf.gz",
        vcf_index = out_dir + "/vcf/{donor}-var.vcf.gz.tbi"
    threads: 16
    log:
        log_dir + "/{donor}_concat.log"
    benchmark:
        bench_dir + "/{donor}_concat.tsv"
    resources:
        mem_mb = mem_large
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools concat -a {input.chunk_zip} -o {output.vcf} --threads {threads} > {log} 2>&1
        tabix -f -p vcf {output.vcf} --threads {threads}
        """
