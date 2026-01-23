#!/usr/bin/env python3
"""
donor_shared_variants.py

Given a directory of bgzipped VCF files (*.vcf.gz), where each VCF is a joint-called donor file
(containing multiple samples), produce for each donor a TSV file containing variant alleles
and how many samples within that donor carry each allele.

Output columns (TSV):
  chrom, pos, end, ref, alt, allele_idx, samples_with_allele, total_samples, sample_list, qual, filter, info

Defaults:
  - Only record alleles present in at least 2 samples (shared). Change with --min-samples.

Notes:
  - Works with multi-allelic sites and handles genotype values like 0/1, 1|1, 0/2, ./., etc.
  - End is computed as pos + max(len(ref), len(alt)) - 1 (useful for indels).
  - No external VCF libraries required; uses gzip for .vcf.gz reading.
"""

import argparse
import gzip
import os
import sys
from pathlib import Path
from typing import List, Tuple, Dict

def parse_args():
    p = argparse.ArgumentParser(description="Produce per-donor shared-variant TSVs from a directory of .vcf.gz files.")
    p.add_argument("vcf_dir", nargs='?',default="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results_ds", help="Directory containing .vcf.gz files (one per donor).")
    p.add_argument("-o", "--outdir", nargs='?', default="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results_ds/within-donor/", help="Output directory (will be created).")
    p.add_argument("--min-samples", nargs='?', type=int, default=2, help="Minimum number of samples that must have the allele to be included (default=2).")
    p.add_argument("--pattern", nargs='?', default="GB*.annotated.vcf.gz", help="Filename pattern to search for in the input directory (default: *.vcf.gz).")
    p.add_argument("--include-filtered", action="store_true", help="Include variants even if FILTER != PASS (default: only include all records regardless of FILTER).")
    p.add_argument("--compress-output", action="store_true", help="Compress output TSV as .tsv.gz")
    return p.parse_args()

def open_vcf(path: Path):
    # Use gzip for .vcf.gz
    return gzip.open(path, "rt") if path.suffixes[-2:] == ['.vcf', '.gz'] or path.suffix == '.gz' else open(path, "r")

def compute_end(pos: int, ref: str, alt: str) -> int:
    # crude end calculation: covers the full length of the longest allele
    return pos + max(len(ref), len(alt)) - 1

def parse_genotype_field(format_field: str, sample_field: str) -> Dict[str,str]:
    """
    Split FORMAT and sample string into a dict of key->value.
    """
    keys = format_field.split(":")
    vals = sample_field.split(":")
    d = {}
    for i,k in enumerate(keys):
        d[k] = vals[i] if i < len(vals) else ""
    return d

def sample_has_allele(gt_str: str, allele_idx_str: str) -> bool:
    """
    Decide whether a sample genotype string (GT) indicates presence of a particular allele.
    gt_str examples: "0/1", "1|1", "0/2", "./.", "0", "1"
    allele_idx_str: e.g. "1", "2"
    Returns True if any allele index in GT equals allele_idx_str.
    Missing genotypes ('.' or './.') return False.
    """
    if not gt_str:
        return False
    # handle phased/unphased separators
    sep = "/" if "/" in gt_str else ("|" if "|" in gt_str else None)
    if sep is None:
        alleles = [gt_str]
    else:
        alleles = gt_str.split(sep)
    for a in alleles:
        if a == "." or a == "":
            continue
        # Some callers use numeric allele indices; others might encode ploidy differently.
        if a == allele_idx_str:
            return True
    return False

def process_vcf(path: Path, outdir: Path, min_samples: int = 2, include_filtered: bool = False, compress_output: bool = False):
    print(f"Processing {path} ...")
    with open_vcf(path) as fh:
        sample_names = []
        header_lines = []
        for line in fh:
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                header_lines.append(line)
                header_cols = line.strip().split("\t")
                # sample columns start at index 9
                if len(header_cols) > 9:
                    sample_names = header_cols[9:]
                else:
                    sample_names = []
                break

        total_samples = len(sample_names)
        if total_samples == 0:
            print(f"Warning: no sample columns found in {path}. Skipping.")
            return

        outname = path.stem  # removes .gz => gives e.g. sample.vcf
        # if file was sample.vcf.gz, path.stem is sample.vcf, so remove .vcf if present
        if outname.endswith(".vcf"):
            outname = outname[:-4]
        outfile = outdir / f"{outname}.shared.tsv"
        if compress_output:
            outfile = Path(str(outfile) + ".gz")

        if compress_output:
            out_fh = gzip.open(outfile, "wt")
        else:
            out_fh = open(outfile, "w")

        header = "\t".join([
            "chrom","pos","end","ref","alt","allele_index","samples_with_allele","total_samples","sample_list","qual","filter","info"
        ]) + "\n"
        out_fh.write(header)

        # iterate remaining variant lines
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0]
            pos = int(cols[1])
            _id = cols[2]
            ref = cols[3]
            alt_field = cols[4]
            qual = cols[5]
            filt = cols[6]
            info = cols[7]
            fmt = cols[8] if len(cols) > 8 else ""
            samples = cols[9:] if len(cols) > 9 else []

            if (not include_filtered) and (filt != "PASS" and filt != "."):
                # if user does not want filtered records, skip
                # NOTE: user can set include_filtered to True to include them
                pass  # we'll still process; remove this if you want to skip filtered
            # Process multi-allelic ALTs separately
            alts = alt_field.split(",") if alt_field != "." else []
            # If there are no ALT alleles, skip
            if not alts:
                continue

            # find index of GT in FORMAT
            fmt_keys = fmt.split(":") if fmt else []
            try:
                gt_index = fmt_keys.index("GT")
            except ValueError:
                gt_index = None

            # For each alt allele, count samples that have it
            for ai, alt in enumerate(alts, start=1):
                allele_idx_str = str(ai)  # 1-based allele index in GT fields
                present_samples = []
                for sname, sfield in zip(sample_names, samples):
                    # extract GT value
                    if gt_index is not None:
                        # careful splitting - sample fields may be "." or fewer subfields
                        s_parts = sfield.split(":")
                        gt = s_parts[gt_index] if gt_index < len(s_parts) else ""
                    else:
                        # no GT available; attempt to infer from SAMPLE column content? skip
                        gt = ""
                    if sample_has_allele(gt, allele_idx_str):
                        present_samples.append(sname)

                n_present = len(present_samples)
                if n_present >= min_samples:
                    end = compute_end(pos, ref, alt)
                    sample_list_str = ",".join(present_samples) if present_samples else ""
                    out_fields = [
                        chrom,
                        str(pos),
                        str(end),
                        ref,
                        alt,
                        allele_idx_str,
                        str(n_present),
                        str(total_samples),
                        sample_list_str,
                        qual,
                        filt,
                        info
                    ]
                    out_fh.write("\t".join(out_fields) + "\n")

        out_fh.close()
        print(f"Wrote {outfile}.")

def main():
    args = parse_args()
    vcf_dir = Path(args.vcf_dir)
    if not vcf_dir.exists() or not vcf_dir.is_dir():
        print(f"Error: {vcf_dir} not found or not a directory.", file=sys.stderr)
        sys.exit(1)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pattern = args.pattern
    vcf_paths = sorted(vcf_dir.glob(pattern))
    if not vcf_paths:
        print(f"No files matching {pattern} in {vcf_dir}.", file=sys.stderr)
        sys.exit(1)

    for p in vcf_paths:
        process_vcf(p, outdir, min_samples=args.min_samples, include_filtered=args.include_filtered, compress_output=args.compress_output)

if __name__ == "__main__":
    main()
