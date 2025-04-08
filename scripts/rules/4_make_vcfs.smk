# Laurel Hiatt 04/07/2025

log_dir = out_dir + "/log/4_make_vcfs/"
bench_dir = out_dir + "/benchmark/4_make_vcfs/"

rule make_bam_list:
    input:
        bam_sort = lambda wildcard: expand(out_dir + "/bam/{sample}-sorted.bam",
            sample=samples),
        bai = lambda wildcard: expand(out_dir + "/bam/{sample}-sorted.bam.bai",
            sample=samples),
    params:
        location = out_dir + "/bam/",
        out_dir = out_dir
    output:
        bam_file_list = temp(out_dir + "/bam/bam_list_{donor}.txt")
    resources:
        mem_mb = mem_small
    threads: 4
    log:
        log_dir + "{donor}_bamlist.log"
    benchmark:
        bench_dir + "{donor}_bamlist.tsv"
    shell:
        """
        python ../variant_calling/create_bam_list.py {wildcards.donor} {params.location} {params.out_dir}
        """

rule generate_regions:
    input:
        index = reference_index,
        chroms = chroms
    params:
        chunks = chunks,
        out_dir = out_dir
    output:
        regions = out_dir + "/regions/chunk.{chroms}.region.{i}.bed"
    resources:
        mem_mb = mem_small
    threads: 2
    conda:
         "../../envs/plot_mosdepth.yaml"
    log:
        log_dir + "{chroms}_{i}_generateregions.log"
    benchmark:
        bench_dir + "{chroms}_{i}_generateregions.tsv"
    shell:
        """
        mkdir -p {input.out_dir}/regions
        python ../variant_calling/fasta_generate_regions.py --fai {input.index} --chunks {params.chunks} --bed {params.out_dir}/regions/chunk.{wildcards.chroms}
        """

# Current filtering:
### min-alternate-count 2
## qsum 40
rule freebayes_variant_calling:
    input:
        bam_file_list = rules.make_bam_list.output.bam_file_list,
        ref = reference,
        regions = rules.generate_regions.output.regions
    output:
        full = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf")
    params:
        out_dir = out_dir
    resources:
        mem_mb = mem_xlarge
    threads: 8
    log:
        log_dir + "{chroms}_{i}_{donor}_freebayes.log"
    benchmark:
        bench_dir + "{chroms}_{i}_{donor}_freebayes.tsv"
    shell:
        """
        mkdir -p {params.out_dir}/vcf/{wildcards.chroms}
        module load freebayes
        freebayes --min-alternate-count 2 --min-alternate-qsum 40 -f {input.ref} -t {input.regions} -L {input.bam_file_list} > {output.full}
        """

rule compress_chunks:
    input:
        chunk_vcf = out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf"
    output:
        chunk_zip = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz")
    threads: 8
    log:
        log_dir + "{chroms}_{i}_{donor}_compress.log"
    benchmark:
        bench_dir + "{chroms}_{i}_{donor}_compress.tsv"
    resources:
        mem_mb = mem_medium
    shell:
        """
        bgzip -f {input.chunk_vcf}
        """

rule index_chunks:
    input:
        chunk_zip = out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz"
    output:
        chunk_zip_index = temp(out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi")
    threads: 8
    log:
        log_dir + "{chroms}_{i}_{donor}_index.log"
    benchmark:
        bench_dir + "{chroms}_{i}_{donor}_index.tsv"
    resources:
        mem_mb = mem_medium
    shell:
        """
        tabix -f -p vcf {input.chunk_zip}
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
        log_dir + "{donor}_concat.log"
    benchmark:
        bench_dir + "{donor}_concat.tsv"
    resources:
        mem_mb = mem_large
    shell:
        """
        module load bcftools
        bcftools concat -a {input.chunk_zip} -o {output.vcf}
        tabix -f -p vcf {output.vcf}
        """
