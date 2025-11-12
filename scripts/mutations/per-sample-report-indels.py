import sys
from cyvcf2 import VCF, Writer
from pyfaidx import Fasta

vcf_path = snakemake.input["vcf"]
fname = snakemake.output["out_vcf"]
high_vaf_threshold = snakemake.params["high_vaf_threshold"]
low_vaf_threshold = snakemake.params["low_vaf_threshold"]
sample_name = snakemake.params["sample_name"]
reference = snakemake.params["ref"]

def count_unique_indels(vcf_path, fname, sample_name):
    vcf = VCF(vcf_path)
    # sample_names = vcf.samples
    # unique_snvs = {sample: 0 for sample in sample_names}
    w = Writer(fname, vcf)

    try:
        sample_idx = vcf.samples.index(sample_name)
    except ValueError:
        raise RuntimeError(f"Sample {sample_name} not found in {vcf_path}")

    passing_indels = []

    for variant in vcf:
        if not variant.is_indel:
            continue  # Skip non-indels

        genotypes = variant.gt_types  # 0=hom-ref, 1=het, 2=unknown, 3=hom-alt

        vafs = variant.gt_alt_freqs

        # Ensure no sample has an unknown genotype (2)
        if 2 in genotypes:
            continue

        # Make sure only sample of interest in non-ref
        other_non_ref = [
            i for i, gt in enumerate(genotypes) if gt in {1, 3} and i != sample_idx
        ]
        if other_non_ref:
            continue  # not unique


        # Count as unique only if specific sample is non-ref
        if genotypes[sample_idx] in {1, 3} and low_vaf_threshold < vafs[sample_idx] < high_vaf_threshold:
            passing_indels.append(
                (variant.CHROM, variant.POS, variant.REF, ",".join(variant.ALT), vafs[sample_idx])
            )
            w.write_record(variant)

    w.close()

    count_unique_indels(vcf_path, fname, sample_name)