rule make_bam_list:
    input:
        bam_sort = lambda wildcard: expand("{out_dir}/bam/{sample}-sorted.bam", sample=samples, out_dir=out_dir),
        bai = lambda wildcard: expand("{out_dir}/bam/{sample}-sorted.bam.bai", sample=samples, out_dir=out_dir),
    params:
        location = "{out_dir}/bam/"
    output:
        bam_file_list = temp("{out_dir}/bam_list_{donor}.txt")
    shell:
        """
        python variant_calling/create_bam_list.py {wildcards.donor} {params.location} {wildcards.out_dir}
        """
rule generate_regions:
    input:
        index = reference_index,
        out_dir = out_dir,
        chroms = chroms
    params:
        chunks = chunks
    output:
        regions = "{out_dir}/regions/chunk.{chroms}.region.{i}.bed"
    shell:
        """
        mkdir -p {input.out_dir}/regions
        python /variant_calling/fasta_generate_regions.py --fai {input.index} --chunks 9 --bed {input.out_dir}/regions/chunk.{wildcards.chroms}
        """

# Current filtering:
### min-alternate-count 2
## qsum 40
rule freebayes_variant_calling:
    input:
        bam_file_list = "{out_dir}/bam_list_{donor}.txt",
        ref = reference,
        regions = "{out_dir}/regions/chunk.{chroms}.region.{i}.bed"
    output:
        full = temp("{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf"),
    resources:
        mem_mb = 64000
    shell:
        """
        mkdir -p {wildcards.out_dir}/vcf/{wildcards.chroms}
        freebayes --min-alternate-count 2 --min-alternate-qsum 40 -f {input.ref} -t {input.regions} -L {input.bam_file_list} > {output.full}
        """

rule compress_chunks:
    input:
        chunk_vcf = "{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf"
    output:
        chunk_zip = temp("{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf.gz")
    shell:
        """
        bgzip -f {input.chunk_vcf}
        """

rule index_chunks:
    input:
        chunk_zip = "{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf.gz"
    output:
        chunk_zip_index = temp("{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi")
    shell:
        """
        tabix -f -p vcf {input.chunk_zip}
        """

rule concat_vcfs:
    input:
        chunk_zip = lambda wildcard: expand("{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf.gz", i = chunks, out_dir = wildcard.out_dir, donor = wildcard.donor, chroms = chroms),
        chunk_vcf_index = lambda wildcard: expand("{out_dir}/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi", i = chunks, out_dir = wildcard.out_dir, donor = wildcard.donor, chroms = chroms),
    output:
        vcf = "{out_dir}/vcf/{donor}-var.vcf.gz",
        vcf_index = "{out_dir}/vcf/{donor}-var.vcf.gz.tbi"
    shell:
        """
        bcftools concat -a {input.chunk_zip} -o {output.vcf}
        tabix -f -p vcf {output.vcf}
        """

