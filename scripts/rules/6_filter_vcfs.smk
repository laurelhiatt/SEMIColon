# Laurel Hiatt 10/27/2025
log_dir = out_dir + "/log/6_filter_vcfs"
bench_dir = out_dir + "/benchmark/6_filter_vcfs"

from collections import OrderedDict

pairs = [
    (get_donor(sample, matches), sample)
    for sample in samples
    if get_donor(sample, matches) != "gibberish"
]

if pairs:
    donors, filtered_samples = zip(*pairs)
else:
    donors, filtered_samples = [], []

pairs_ds = [
    (get_donor_only_crypt(sample, matches), sample)
    for sample in samples
    if not sample.startswith("GIB") and get_donor_only_crypt(sample, matches) is not None
]

if pairs_ds:
    donors_ds, filtered_samples_ds = zip(*pairs_ds)
else:
    donors_ds, filtered_samples_ds = (), ()



def get_all_samples_for_donor(donor, matches):
    """
    Returns a list of samples for a donor:
    - All crypt_samples
    - Blood sample(s), if present, only if not already in crypt_samples
    """
    rec = matches.get(donor, {})
    samples = list(rec.get("crypt_samples", []))  # copy list

    bs = rec.get("blood_sample", None)
    if bs:
        if isinstance(bs, str) and bs not in samples:
            samples.append(bs)
        elif isinstance(bs, list):
            for b in bs:
                if b not in samples:
                    samples.append(b)

    return samples


def get_all_samples_for_donor_ds(donor, matches):
    """
    Returns a list of samples for a donor:
    - All crypt_samples
    - Blood sample(s), if present, only if not already in crypt_samples
    """
    rec = matches.get(donor, {})
    samples = list(rec.get("crypt_samples", []))  # copy list
    return samples

def get_all_inputs(donors, matches, out_dir):
    """
    Returns all donor filtered VCF paths using the real absolute path.
    Removes duplicates to avoid periodic wildcard detection hang.
    """
    all_inputs = []
    for donor in donors:
        donor_samples = get_all_samples_for_donor(donor, matches)
        donor_paths = [
            f"{out_dir}/results/{donor}/{sample}.filtered.vcf.gz"
            for sample in donor_samples
        ]
        all_inputs.extend(donor_paths)

    from collections import OrderedDict
    all_inputs = list(OrderedDict.fromkeys(all_inputs))

    # Debug print
    #print("DEBUG: find_recurrent input files:", all_inputs)

    return all_inputs


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


# decompose the vcf for downstream
rule decompose_vcfs_ds:
    input:
        vcf = out_dir + "/vcf/merged/{donor}.GRCh38.joint_genotyped.vcf.gz",
        fasta = reference,
        index = rules.index_deepsomatic_donor_vcfs.output.joint_index
    output:
        clean_vcf = temp(out_dir + "/vcf/merged/{donor}-clean-var.vcf.gz")
    resources:
        mem_mb = mem_large
    threads: 8
    envmodules:
        "bcftools/1.16"
    log:
        log_dir + "/{donor}_decompose_ds.log"
    benchmark:
        bench_dir + "/{donor}_decompose_ds.tsv"
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
        annotated_vcf = temp(out_dir + "/vcf/{donor}-annotated-var.vcf.gz")
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

### it shouldn't exist in gnomad
rule filter_by_gnomad:
    input:
        filtered_vcf= out_dir + "/vcf/{donor}-annotated-var-noLCR.vcf.gz"
    output:
        gnomad_vcf= out_dir + "/results/{donor}.gnomad.filtered.vcf.gz"
    shell:
        """
        module load bcftools
        bcftools filter -i '(gnomad_popmax_af <= 0 || gnomad_popmax_af == ".")' -Oz -o {output.gnomad_vcf} {input.filtered_vcf}
        """

rule annotate_gene:
    input:
        gnomad_vcf= out_dir + "/results/{donor}.gnomad.filtered.vcf.gz",
    output:
        vcf = out_dir + "/results/{donor}.annotated.vcf.gz",
        done = out_dir + "/results/{donor}/gnomad.done"
    envmodules:
        "vep/104.2"
    resources:
        mem_mb = mem_large
    params:
        dir_cache = "/scratch/ucgd/lustre/common/data/vep_cache"
    log:
        log_dir + "/{donor}_annotate.log"
    threads:
        8
    envmodules:
        "vep/104.2"
    shell:
        """
        module load vep/104.2
        echo "Running VEP annotation for {wildcards.donor}"
        vep --cache --dir_cache {params.dir_cache} -i {input.gnomad_vcf} --vcf --compress_output bgzip -o {output.vcf} --symbol --force_overwrite --fork 10
        tabix -p vcf {output.vcf}
        echo "Finished VEP annotation for {wildcards.donor}"
        touch {output.done}
        """

def donor_done_input(wildcards):
    donor = wildcards.donor
    sample = wildcards.sample
    donor_info = matches.get(donor, {})
    crypt_samples = donor_info.get("crypt_samples", [])
    blood_sample = donor_info.get("blood_sample")

    # If sample is one of the crypt samples or matches the blood sample
    if sample in crypt_samples or sample == blood_sample:
        return os.path.join(out_dir, "results", donor, "gnomad.done")

    # Otherwise, skip
    return None

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
        sample_vcf = temp(out_dir + "/results/{donor}/{sample}.filtered_noAD.vcf.gz")
    params:
        annotated_vcf = out_dir + "/results/{donor}.annotated.vcf.gz"
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua \
            -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 0 end, sampleIndexes, variant:format("AD"))' \
            -o {output.sample_vcf} {params.annotated_vcf}
        """

rule filter_by_depth:
    input:
        sample_vcf = out_dir + "/results/{donor}/{sample}.filtered_noAD.vcf.gz",
        lua = rules.make_lua.output.lua  # include your lua so sampleIndexes is available
    output:
        depth_vcf= temp(out_dir + "/results/{donor}/{sample}.depth.filtered.vcf.gz")
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua \
          -e 'dps = variant:format("DP"); return any(function(i) local dp = dps[i]; return dp and dp > 7 end, sampleIndexes)' \
          -o {output.depth_vcf} {input.sample_vcf}
        """

#### replace above with this between runs
#
###
rule index_depth:
    input:
        depth_vcf= out_dir + "/results/{donor}/{sample}.depth.filtered.vcf.gz"
    output:
        indexed_vcf= temp(out_dir + "/results/{donor}/{sample}.depth.filtered.vcf.gz.tbi")
    threads:
        8
    shell:
        """
        tabix -f -p vcf {input.depth_vcf} --threads {threads}
        """

### this sample must have > # alternate allele
rule filter_by_alt_depth:
    input:
        sample_vcf = out_dir + "/results/{donor}/{sample}.depth.filtered.vcf.gz",
        index_vcf = out_dir + "/results/{donor}/{sample}.depth.filtered.vcf.gz.tbi",
        lua = rules.make_lua.output.lua,
    output:
        vcf= out_dir + "/results/{donor}/{sample}.filtered.vcf.gz",
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 3 end, sampleIndexes, variant:format("AD"))' -o {output.vcf} {input.sample_vcf}
        """

rule find_recurrent:
    input:
        # Lazy evaluation to ensure DAG sees dependencies
        all_filtered=lambda wc: get_all_inputs(donors, matches, out_dir)
    params:
        min_recurrence = 2
    output:
        recurrent_vcf =  out_dir + "/" + "results/recurrent.vcf",
        report= out_dir + "/" + "results/recurrent.txt"
    conda:
        "../../envs/recurrent.yaml"
    script:
        "../filtering/find_recurrent.py"


rule compress_recurrent:
    input:
        recurrent = out_dir + "/" + "results/recurrent.vcf"
    output:
        recurrent_vcf = out_dir + "/" + "results/recurrent.vcf.gz"
    threads:
        8
    shell:
        """
        module load bcftools
        bcftools sort -Oz -o {output.recurrent_vcf} {input.recurrent}
        tabix -p vcf {output.recurrent_vcf}
        """

rule filter_by_recurrent:
    input:
        sample_vcf = out_dir + "/results/{donor}/{sample}.filtered.vcf.gz",
        recurrent_vcf = out_dir + "/results/recurrent.vcf.gz"
    output:
        vcf= out_dir + "/results/{donor}/{sample}.vcf.gz"
    shell:
        """
        tabix -f -p vcf {input.sample_vcf}
        module load bcftools
        bcftools isec -C -w1 -O z -o {output.vcf} {input.sample_vcf} {input.recurrent_vcf}
        """

rule count_indels:
    input:
        vcf= out_dir + "/results/{donor}/{sample}.vcf.gz"
    output:
        out_vcf= out_dir + "/results/{donor}/{sample}.indels.vcf.gz",
    params:
        sample_name = "{donor}_{sample}",
        ref = reference,
        high_vaf_threshold = 1.1,
        low_vaf_threshold = 0.0
    conda:
        "../../envs/cyvcf2.yaml"
    script:
        "../mutations/per-sample-report-indels.py"

rule count_snvs:
    input:
        vcf= out_dir + "/results/{donor}/{sample}.vcf.gz"
    output:
        out_vcf= out_dir + "/results/{donor}/{sample}.snvs.vcf.gz",
        txt= out_dir + "/results/{donor}/{sample}.snv_count.txt"
    params:
        sample_name = "{donor}_{sample}",
        ref = reference,
        high_vaf_threshold = 1.1,
        low_vaf_threshold = 0.0
    conda:
        "../../envs/cyvcf2.yaml"
    script:
        "../mutations/per-sample-report.py"


# annotate the vcf with gnomAD
rule gnomad_VCFs_ds:
    input:
        clean_vcf = out_dir + "/vcf/merged/{donor}-clean-var.vcf.gz",
    output:
        annotated_vcf = temp(out_dir + "/vcf/merged/{donor}-annotated-var.vcf.gz")
    resources:
        mem_mb = mem_large
    threads: 2
    log:
        log_dir + "/{donor}_gnomad_ds.log"
    envmodules:
        "slivar/0.3.1"
    shell:
        """
        slivar expr -g /scratch/ucgd/lustre/common/data/Slivar/db/gnomad.hg38.genomes.v3.fix.zip -v {input.clean_vcf} -o {output.annotated_vcf}
        """

rule remove_lcr_ds:
    input:
        annotated_vcf = out_dir + "/vcf/merged/{donor}-annotated-var.vcf.gz",
        lcr_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/LCR-hs38.bed.gz",
        simplerepeats_bed = "/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/GRCh38.UCSC.SimpleRepeats.bed.gz"
    output:
        filtered_vcf = out_dir + "/vcf/merged/{donor}-annotated-var-noLCR.vcf.gz"

    resources:
        mem_mb = mem_xlarge
    threads: 4
    log:
        log_dir + "/{donor}_noLCR_ds.log"
    benchmark:
        bench_dir + "/{donor}_noLCR_ds.tsv"
    shell:
        """
        module load bedtools
        bedtools intersect -header -v -a {input.annotated_vcf} -b {input.lcr_bed} | \
        bedtools intersect -header -v -a - -b {input.simplerepeats_bed} | bgzip -c > {output.filtered_vcf}
        tabix -p vcf {output.filtered_vcf}
        """

rule filter_by_gnomad_ds:
    input:
        filtered_vcf= out_dir + "/vcf/merged/{donor}-annotated-var-noLCR.vcf.gz"
    output:
        gnomad_vcf= out_dir + "/results_ds/{donor}.gnomad.filtered.vcf.gz"
    shell:
        """
        module load bcftools
        bcftools filter -i '(gnomad_popmax_af <= 0 || gnomad_popmax_af == ".")' -Oz -o {output.gnomad_vcf} {input.filtered_vcf}
        """

rule annotate_gene_ds:
    input:
        gnomad_vcf= out_dir + "/results_ds/{donor}.gnomad.filtered.vcf.gz",
    output:
        vcf = out_dir + "/results_ds/{donor}.annotated.vcf.gz",
        done = out_dir + "/results_ds/{donor}/gnomad_ds.done"
    envmodules:
        "vep/104.2"
    resources:
        mem_mb = mem_large
    params:
        dir_cache = "/scratch/ucgd/lustre/common/data/vep_cache"
    log:
        log_dir + "/{donor}_annotate_ds.log"
    threads:
        8
    envmodules:
        "vep/104.2"
    shell:
        """
        module load vep/104.2
        echo "Running VEP annotation for {wildcards.donor}"
        vep --cache --dir_cache {params.dir_cache} -i {input.gnomad_vcf} --vcf --compress_output bgzip -o {output.vcf} --symbol --force_overwrite --fork 10
        tabix -p vcf {output.vcf}
        echo "Finished VEP annotation for {wildcards.donor}"
        touch {output.done}
        """

def donor_done_input_ds(wildcards):
    donor = wildcards.donor
    sample = wildcards.sample
    donor_info = matches.get(donor, {})
    crypt_samples = donor_info.get("crypt_samples", [])

    if donor not in matches:
        raise ValueError(f"donor '{donor}' not found in matches")

    if sample not in crypt_samples:
        print(f"Sample '{sample}' not found in crypt_samples for donor '{donor}'")

    # return a single path (string) — valid return type for Snakemake input functions
    return os.path.join(out_dir, "results_ds", donor, "gnomad_ds.done")

rule make_lua_ds:
    input:
        gnomad_done = donor_done_input_ds
    output:
        lua = out_dir + "/results_ds/{donor}/{sample}_soi_ds.lua"
    run:
        if wildcards.sample.startswith("GIB"):
            # mark rule as complete without doing anything
            print("Skipping blood sample:", wildcards.sample)

        with open(output.lua, "w") as f:
            f.write(f'samplesOfInterest = {{"{wildcards.sample}"}}\n')


rule filter_by_sample_ds:
    input:
        lua = rules.make_lua_ds.output.lua,
    output:
        sample_vcf = temp(out_dir + "/results_ds/{donor}/{sample}.filtered_noAD.vcf.gz")
    params:
        annotated_vcf = out_dir + "/results_ds/{donor}.annotated.vcf.gz"
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua \
            -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 0 end, sampleIndexes, variant:format("AD"))' \
            -o {output.sample_vcf} {params.annotated_vcf}
        """

rule filter_by_depth_ds:
    input:
        sample_vcf = out_dir + "/results_ds/{donor}/{sample}.filtered_noAD.vcf.gz",
        lua = rules.make_lua_ds.output.lua  # include your lua so sampleIndexes is available
    output:
        depth_vcf= temp(out_dir + "/results_ds/{donor}/{sample}.depth.filtered.vcf.gz")
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua \
          -e 'dps = variant:format("DP"); return any(function(i) local dp = dps[i]; return dp and dp > 7 end, sampleIndexes)' \
          -o {output.depth_vcf} {input.sample_vcf}
        """

rule index_depth_ds:
    input:
        depth_vcf= out_dir + "/results_ds/{donor}/{sample}.depth.filtered.vcf.gz"
    output:
        indexed_vcf= temp(out_dir + "/results_ds/{donor}/{sample}.depth.filtered.vcf.gz.tbi")
    threads:
        8
    shell:
        """
        tabix -f -p vcf {input.depth_vcf} --threads {threads}
        """

rule filter_by_alt_depth_ds:
    input:
        sample_vcf = out_dir + "/results_ds/{donor}/{sample}.depth.filtered.vcf.gz",
        index_vcf = out_dir + "/results_ds/{donor}/{sample}.depth.filtered.vcf.gz.tbi",
        lua = rules.make_lua_ds.output.lua,
    output:
        vcf= out_dir + "/results_ds/{donor}/{sample}.filtered.vcf.gz",
    shell:
        """
        ./vcfexpress filter -p {input.lua} -p /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/config/sample-groups.lua -e 'return all_none(function(ad) return #ad > 1 and ad[2] > 3 end, sampleIndexes, variant:format("AD"))' -o {output.vcf} {input.sample_vcf}
        """

rule index_by_alt_depth_ds:
    input:
        vcf= out_dir + "/results_ds/{donor}/{sample}.filtered.vcf.gz"
    output:
        index= out_dir + "/results_ds/{donor}/{sample}.filtered.vcf.gz.tbi"
    threads:
        8
    shell:
        """
        tabix -f -p vcf {input.vcf} --threads {threads}
        """

def get_all_inputs_ds(donors, matches, out_dir):
    """
    Returns all donor filtered VCF paths using the real absolute path.
    Removes duplicates to avoid periodic wildcard detection hang.
    """
    all_inputs = []
    for donor in donors:
        donor_samples = get_all_samples_for_donor_ds(donor, matches)
        donor_paths = [
            f"{out_dir}/results_ds/{donor}/{sample}.filtered.vcf.gz"
            for sample in donor_samples
        ]
        all_inputs.extend(donor_paths)

    from collections import OrderedDict
    all_inputs = list(OrderedDict.fromkeys(all_inputs))

    return all_inputs


rule find_recurrent_ds:
    input:
        all_filtered=lambda wc: get_all_inputs_ds(donors, matches, out_dir),
    params:
        min_recurrence = 2
    output:
        recurrent_vcf =  out_dir + "/" + "results_ds/recurrent.vcf",
        report= out_dir + "/" + "results_ds/recurrent.txt"
    conda:
        "../../envs/recurrent.yaml"
    script:
        "../filtering/find_recurrent.py"


rule compress_recurrent_ds:
    input:
        recurrent = out_dir + "/" + "results_ds/recurrent.vcf"
    output:
        recurrent_vcf = out_dir + "/" + "results_ds/recurrent.vcf.gz"
    threads:
        8
    shell:
        """
        module load bcftools
        bcftools sort -Oz -o {output.recurrent_vcf} {input.recurrent}
        tabix -p vcf {output.recurrent_vcf}
        """

rule filter_by_recurrent_ds:
    input:
        sample_vcf = out_dir + "/results_ds/{donor}/{sample}.filtered.vcf.gz",
        index = out_dir + "/results_ds/{donor}/{sample}.filtered.vcf.gz.tbi",
        recurrent_vcf = out_dir + "/results_ds/recurrent.vcf.gz"
    output:
        vcf= out_dir + "/results_ds/{donor}/{sample}.vcf.gz"
    shell:
        """
        tabix -p vcf {input.sample_vcf}
        module load bcftools
        bcftools isec -C -w1 -O z -o {output.vcf} {input.sample_vcf} {input.recurrent_vcf}
        """

rule count_indels_ds:
    input:
        vcf= out_dir + "/results_ds/{donor}/{sample}.vcf.gz"
    output:
        out_vcf= out_dir + "/results_ds/{donor}/{sample}.indels.vcf.gz",
    params:
        sample_name = "{sample}",
        ref = reference,
        high_vaf_threshold = 1.1,
        low_vaf_threshold = 0.0
    conda:
        "../../envs/cyvcf2.yaml"
    script:
        "../mutations/per-sample-report-indels.py"

rule count_snvs_ds:
    input:
        vcf= out_dir + "/results_ds/{donor}/{sample}.vcf.gz"
    output:
        out_vcf= out_dir + "/results_ds/{donor}/{sample}.snvs.vcf.gz",
        txt= out_dir + "/results_ds/{donor}/{sample}.snv_count.txt"
    params:
        sample_name = "{sample}",
        ref = reference,
        high_vaf_threshold = 1.1,
        low_vaf_threshold = 0.0
    conda:
        "../../envs/cyvcf2.yaml"
    script:
        "../mutations/per-sample-report.py"