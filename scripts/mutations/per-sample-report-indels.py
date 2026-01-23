# import sys
# from cyvcf2 import VCF, Writer
# from pyfaidx import Fasta

# vcf_path = snakemake.input["vcf"]
# fname = snakemake.output["out_vcf"]
# high_vaf_threshold = snakemake.params["high_vaf_threshold"]
# low_vaf_threshold = snakemake.params["low_vaf_threshold"]
# sample_name = snakemake.params["sample_name"]
# reference = snakemake.params["ref"]

# def count_unique_indels(vcf_path, fname, sample_name):
#     vcf = VCF(vcf_path)
#     # sample_names = vcf.samples
#     # unique_snvs = {sample: 0 for sample in sample_names}
#     w = Writer(fname, vcf)

#     try:
#         sample_idx = vcf.samples.index(sample_name)
#     except ValueError:
#         raise RuntimeError(f"Sample {sample_name} not found in {vcf_path}")

#     passing_indels = []

#     for variant in vcf:
#         if not variant.is_indel:
#             continue  # Skip non-indels

#         genotypes = variant.gt_types  # 0=hom-ref, 1=het, 2=unknown, 3=hom-alt

#         vafs = variant.gt_alt_freqs

#         # Ensure no sample has an unknown genotype (2)
#         #if 2 in genotypes:
#         #    continue
#         if (genotypes == 2).sum() > 1:
#             continue

#         # Make sure only sample of interest in non-ref
#         other_non_ref = [
#             i for i, gt in enumerate(genotypes) if gt in {1, 3} and i != sample_idx
#         ]
#         if other_non_ref:
#             continue  # not unique


#         # Count as unique only if specific sample is non-ref
#         if genotypes[sample_idx] in {1, 3} and low_vaf_threshold < vafs[sample_idx] < high_vaf_threshold:
#             passing_indels.append(
#                 (variant.CHROM, variant.POS, variant.REF, ",".join(variant.ALT), vafs[sample_idx])
#             )
#             w.write_record(variant)

#     w.close()

# count_unique_indels(vcf_path, fname, sample_name)



import sys
from cyvcf2 import VCF, Writer
from pyfaidx import Fasta

# Path/params from Snakemake
vcf_path = snakemake.input["vcf"]
sample_name = snakemake.params["sample_name"]
reference = snakemake.params["ref"]
high_vaf_threshold = snakemake.params["high_vaf_threshold"]
low_vaf_threshold = snakemake.params["low_vaf_threshold"]

def _get_out(key, ext):
    try:
        return snakemake.output[key]
    except Exception:
        return f"{sample_name}.{key}.{ext}"

out_vcfs = {
    "nuclear_all": _get_out("nuclear_all_vcf", "vcf.gz"),
    "nuclear_clonal": _get_out("nuclear_clonal_vcf", "vcf.gz"),
    "nuclear_subclonal": _get_out("nuclear_subclonal_vcf", "vcf.gz"),
    "mito_all": _get_out("mito_all_vcf", "vcf.gz"),
    "mito_clonal": _get_out("mito_clonal_vcf", "vcf.gz"),
    "mito_subclonal": _get_out("mito_subclonal_vcf", "vcf.gz"),
}

out_txts = {
    "nuclear_all": _get_out("nuclear_all_txt", "txt"),
    "nuclear_clonal": _get_out("nuclear_clonal_txt", "txt"),
    "nuclear_subclonal": _get_out("nuclear_subclonal_txt", "txt"),
    "mito_all": _get_out("mito_all_txt", "txt"),
    "mito_clonal": _get_out("mito_clonal_txt", "txt"),
    "mito_subclonal": _get_out("mito_subclonal_txt", "txt"),
}

MITO_NAMES = {"MT", "chrM", "M", "mitochondrion"}

def count_and_split(vcf_path, sample_name):
    vcf = VCF(vcf_path)

    try:
        sample_idx = vcf.samples.index(sample_name)
    except ValueError:
        raise RuntimeError(f"Sample {sample_name} not found in {vcf_path}")

    writers = {k: Writer(path, vcf) for k, path in out_vcfs.items()}

    genome = Fasta(reference, as_raw=True)

    passing = {k: [] for k in out_txts.keys()}

    for variant in vcf:
        if not variant.is_indel:
            continue

        chrom = variant.CHROM
        pos = variant.POS

        genotypes = variant.gt_types
        if (genotypes == 2).sum() > 1:
            continue

        other_non_ref = [i for i, gt in enumerate(genotypes) if gt in {1, 3} and i != sample_idx]
        if other_non_ref:
            continue

        if genotypes[sample_idx] not in {1, 3}:
            continue

        vafs = variant.gt_alt_freqs
        vaf = None
        if vafs is not None and len(vafs) > sample_idx:
            try:
                vaf = float(vafs[sample_idx])
            except Exception:
                vaf = None

        alt_str = ",".join(variant.ALT)
        is_mito = chrom in MITO_NAMES

        key_all = "mito_all" if is_mito else "nuclear_all"
        writers[key_all].write_record(variant)
        passing[key_all].append((chrom, pos, variant.REF, alt_str, vaf))

        # Clonal: low_vaf_threshold < vaf <= high_vaf_threshold (inclusive upper bound)
        if vaf is not None and (low_vaf_threshold < vaf <= high_vaf_threshold):
            key_clonal = "mito_clonal" if is_mito else "nuclear_clonal"
            writers[key_clonal].write_record(variant)
            passing[key_clonal].append((chrom, pos, variant.REF, alt_str, vaf))

        # Subclonal: inclusive lower (vaf <= low_vaf_threshold)
        if vaf is not None and vaf <= low_vaf_threshold:
            key_sub = "mito_subclonal" if is_mito else "nuclear_subclonal"
            writers[key_sub].write_record(variant)
            passing[key_sub].append((chrom, pos, variant.REF, alt_str, vaf))

    for w in writers.values():
        w.close()

    def _write_txt(path, records):
        with open(path, "w") as f:
            f.write(f"# {sample_name} ({path})\n")
            f.write(f"# VAF filter: clonal when {low_vaf_threshold} < VAF <= {high_vaf_threshold}\n")
            f.write("CHROM\tPOS\tREF\tALT\tVAF\n")
            for chrom, pos, ref, alt, vaf in records:
                vaf_str = f"{vaf:.4f}" if (vaf is not None) else "NA"
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vaf_str}\n")
            f.write(f"\n# Total: {len(records)}\n")

    for k, path in out_txts.items():
        _write_txt(path, passing[k])

# Run
count_and_split(vcf_path, sample_name)
