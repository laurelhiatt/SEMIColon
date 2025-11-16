# Laurel Hiatt 09/10/2025
from itertools import product

log_dir = out_dir + "/log/4_make_vcfs"
bench_dir = out_dir + "/benchmark/4_make_vcfs"

def swap_acct(wc, attempt):
    if attempt == 1:
        return "owner-guest"
    else:
        return "quinlan-rw"

# If first attempt, redwood-guest (matches owner-guest). Else, quinlan.
def swap_part(wc, attempt):
    if attempt == 1:
        return "redwood-shared-guest"
    else:
        return "quinlan-shared-rw"

# If first attempt, 3 hours (short time for guest). Else, two weeks.
def swap_time(wc, attempt):
    if attempt == 1:
        return 600
    else:
        return 20160

# def get_bam_inputs(wildcards):
#     donor_samples = matches.get(wildcards.donor, {}).get("crypt_samples", [])
#     return [
#         os.path.join(out_dir, "bam", f"{s}-sorted.bam") for s in donor_samples
#     ]

def get_bam_inputs(wildcards):
    donor_info = matches.get(wildcards.donor, {})
    bam_files = []

    # Add crypt samples if present
    crypt_samples = donor_info.get("crypt_samples", [])
    bam_files.extend(
        os.path.join(out_dir, "bam", f"{s}-sorted.bam") for s in crypt_samples
    )

    # Add blood sample if present
    blood_sample = donor_info.get("blood_sample")
    if blood_sample:
        bam_files.append(os.path.join(out_dir, "bam", f"{blood_sample}-sorted.bam"))

    return bam_files


# def get_bai_inputs(wildcards):
#     donor_samples = matches.get(wildcards.donor, {}).get("crypt_samples", [])
#     return [
#         os.path.join(out_dir, "bam", f"{s}-sorted.bam.bai") for s in donor_samples
#     ]

def get_bai_inputs(wildcards):
    donor_info = matches.get(wildcards.donor, {})
    bai_files = []

    # Add crypt samples if present
    crypt_samples = donor_info.get("crypt_samples", [])
    bai_files.extend(
        os.path.join(out_dir, "bam", f"{s}-sorted.bam.bai") for s in crypt_samples
    )

    # Add blood sample if present
    blood_sample = donor_info.get("blood_sample")
    if blood_sample:
        bai_files.append(os.path.join(out_dir, "bam", f"{blood_sample}-sorted.bam.bai"))

    return bai_files

def make_chroms_dict(directory_path, chroms):
    chrom_chunks = {}
    try:
        for chrom in chroms:
            for filename in os.listdir(directory_path):
                match = re.search(rf"chunk\.{re.escape(chrom)}\.region\.(\d+)\.bed", filename)
                if match:
                    chunk_number = int(match.group(1))
                    if chrom not in chrom_chunks:
                        chrom_chunks[chrom] = []
                    chrom_chunks[chrom].append(chunk_number)

        # Optional: sort the chunk numbers for each chrom
        for chrom in chrom_chunks:
            chrom_chunks[chrom].sort()

    except FileNotFoundError:
        print(f"Error: Directory not found at {directory_path}")
    except Exception as e:
        print(f"An error occurred: {e}")

    return chrom_chunks

def chunks_n_chroms(*iterables):
    chroms = [i for i in iterables if i[0][0] == "chroms"][0]
    chunks = [i for i in iterables if i[0][0] == "i"][0]
    donor = [i for i in iterables if i[0][0] == "donor"][0] * len(chroms)
    assert len(chroms) == len(chunks)

    combos = list(zip(donor, chroms, chunks))
    ret_list = []
    for combo in combos:
        chunks = list(product([combo[2][0]], combo[2][1]))
        ret_list.extend(product([combo[0]], [combo[1]], chunks))

    return ret_list

chromosome_chunks_dict = make_chroms_dict("/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/regions", chroms = chroms)

rule make_bam_list:
    input:
        bam_sort = get_bam_inputs,
        bai = get_bai_inputs,
    params:
        out_dir = out_dir,
        matches = matches
    output:
        bam_file_list = temp(out_dir + "/bam/bam_list_{donor}.txt")
    resources:
        mem_mb = mem_xsmall
    threads: 2
    localrule: True
    script:
        "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/scripts/variant_calling/create_bam_list.py"

rule generate_regions:
    input:
        chroms = chroms
    params:
        index = reference_index,
        out_dir = out_dir
    output:
        regions = out_dir + "/regions/chunk.{chroms}.region.{i}.bed"
    resources:
        mem_mb = mem_small
    threads: 2
    localrule: True
    conda:
         "../../envs/plot_mosdepth.yaml"
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
        mem_mb = mem_xlarge,
        slurm_account = swap_acct,
        slurm_partition = swap_part,
        runtime = swap_time
    threads:
        1
    log:
        log_dir + "/{chroms}_{i}_{donor}_freebayes.log"
    benchmark:
        bench_dir + "/{chroms}_{i}_{donor}_freebayes.tsv"
    envmodules:
        "freebayes/1.3.4"
    shell:
        """
        module load freebayes/1.3.4
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
        chunk_zip = lambda wildcard: expand(
            out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz",
            chunks_n_chroms,
            donor = wildcard.donor,
            chroms = chromosome_chunks_dict.keys(),
            i = chromosome_chunks_dict.values()
        ),
        chunk_vcf_index = lambda wildcard: expand(
            out_dir + "/vcf/{chroms}/{donor}-variants.{i}.vcf.gz.tbi",
            chunks_n_chroms,
            donor = wildcard.donor,
            chroms = chromosome_chunks_dict.keys(),
            i = chromosome_chunks_dict.values()
        ),
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
        "bcftools/1.16"
    shell:
        """
        module load bcftools/1.16
        bcftools concat -a {input.chunk_zip} -o {output.vcf} --threads {threads} > {log} 2>&1
        tabix -f -p vcf {output.vcf} --threads {threads}
        """
