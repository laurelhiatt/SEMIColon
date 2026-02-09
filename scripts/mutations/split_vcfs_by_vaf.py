#!/usr/bin/env python3

import os
import sys
import re
from cyvcf2 import VCF, Writer

VAF_THRESHOLD = 0.2


def infer_sample_name(filename):
    """
    Sample name = everything before '.nuclear*'
    """
    match = re.match(r"(.+?)\.nuclear.*", filename)
    if not match:
        raise RuntimeError(
            f"Could not infer sample name from filename: {filename}"
        )
    return match.group(1)


def get_vaf_for_sample(variant, sample_idx):
    """
    Safely extract VAF for a specific sample index.
    """
    vafs = variant.gt_alt_freqs
    if vafs is None or len(vafs) <= sample_idx:
        return None
    try:
        return float(vafs[sample_idx])
    except Exception:
        return None


def split_vcf(vcf_path, out_dir):
    base = os.path.basename(vcf_path)
    sample_name = infer_sample_name(base)

    vcf = VCF(vcf_path)

    try:
        sample_idx = vcf.samples.index(sample_name)
    except ValueError:
        raise RuntimeError(
            f"Sample '{sample_name}' not found in VCF samples: {vcf.samples}"
        )

    subclonal_path = os.path.join(
        out_dir, f"{sample_name}.subclonal.vcf"
    )
    other_path = os.path.join(
        out_dir, f"{sample_name}.other.vcf"
    )

    sub_writer = Writer(subclonal_path, vcf)
    other_writer = Writer(other_path, vcf)

    n_sub = 0
    n_other = 0

    for variant in vcf:
        vaf = get_vaf_for_sample(variant, sample_idx)

        if vaf is not None and vaf < VAF_THRESHOLD:
            sub_writer.write_record(variant)
            n_sub += 1
        else:
            other_writer.write_record(variant)
            n_other += 1

    sub_writer.close()
    other_writer.close()
    vcf.close()

    print(
        f"{base} → {sample_name} | "
        f"subclonal: {n_sub} | other: {n_other}"
    )


def main(vcf_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    vcfs = [
        os.path.join(vcf_dir, f)
        for f in os.listdir(vcf_dir)
        if f.endswith(".vcf")
    ]

    if not vcfs:
        sys.exit("No VCF files found.")

    for vcf_path in sorted(vcfs):
        split_vcf(vcf_path, out_dir)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(
            f"Usage: {sys.argv[0]} <vcf_directory> <output_directory>"
        )

    main(sys.argv[1], sys.argv[2])
