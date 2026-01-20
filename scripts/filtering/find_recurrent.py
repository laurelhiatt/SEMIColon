import pysam
import os
from collections import defaultdict
from datetime import datetime

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
input_files = snakemake.input["all_filtered"]
min_recurrence= snakemake.params.get("min_recurrence", 2)
report = snakemake.output["report"]
output_file = snakemake.output["recurrent_vcf"]

def create_recurrent_vcf(input_files, output_file, min_recurrence, report):
    """
    Creates a VCF with only recurrent variants.
    Removes all sample/genotype columns.
    Adds RECURRENCE + FILES fields to INFO.
    Produces a text report listing variants and files they were found in.
    """

    if not input_files:
        print("No input VCF files provided.")
        return

    # ---------------------------
    # Build a clean header
    # ---------------------------
    first_vcf = pysam.VariantFile(input_files[0])
    header = pysam.VariantHeader()

    # Copy contigs
    for ctg in first_vcf.header.contigs:
        header.contigs.add(ctg, length=first_vcf.header.contigs[ctg].length)

    # Copy INFO definitions
    for info_id in first_vcf.header.info:
        header.info.add(
            info_id,
            first_vcf.header.info[info_id].number,
            first_vcf.header.info[info_id].type,
            first_vcf.header.info[info_id].description
        )

    # ---- Copy FILTER definitions so we can re-use record filters ----
    # This fixes the KeyError: 'Invalid filter: ...' when adding filters to new_rec
    try:
        for fid in first_vcf.header.filters:
            fobj = first_vcf.header.filters[fid]
            # add(filter_id, number, type, description) - number/type may be None in some versions
            # Use safe defaults if attributes missing
            number = getattr(fobj, "number", ".")
            ftype = getattr(fobj, "type", "String")
            desc = getattr(fobj, "description", "")
            # pysam header.filters.add expects (id, number, type, description)
            # number should be string-like per pysam, but passing as-is usually works
            header.filters.add(fid, number, ftype, desc)
    except Exception:
        # If header.filters not available or copying fails, continue without raising;
        # we'll still handle unknown filters below when adding them to records.
        pass

    # Add custom INFO fields
    header.add_line('##INFO=<ID=RECURRENCE,Number=1,Type=Integer,Description="Number of samples this variant was found in.">')
    header.add_line('##INFO=<ID=FILES,Number=.,Type=String,Description="List of input files containing this variant.">')

    # Add FORMAT column (required, even if no samples)
    header.formats.add("GT", 1, "String", "Dummy genotype format (no samples)")

    # Define columns (no samples!)
    header.add_meta("fileformat", value="VCFv4.2")
    header.add_sample("DUMMY")  # We need one dummy sample to keep FORMAT column; we’ll drop its genotype later.

    # -----------------------------------
    # Gather variant counts + file list
    # -----------------------------------
    variant_data = defaultdict(lambda: {"count": 0, "record": None, "files": set()})

    for vcf_file in input_files:
        vcf_in = pysam.VariantFile(vcf_file)
        for record in vcf_in:
            # record.alts may be None for some records — guard against that
            if record.alts is None:
                continue
            for alt in record.alts:
                key = (record.chrom, record.pos, record.ref, alt)

                # Store first record for INFO copying
                if variant_data[key]["record"] is None:
                    variant_data[key]["record"] = record

                variant_data[key]["count"] += 1
                variant_data[key]["files"].add(vcf_file)

        vcf_in.close()

    # -----------------------------------
    # Write output VCF (no samples)
    # -----------------------------------
    with pysam.VariantFile(output_file, "w", header=header) as vcf_out:

        for (chrom, pos, ref, alt), info in variant_data.items():

            if info["count"] < min_recurrence:
                continue

            old_record = info["record"]

            # Create record
            new_rec = vcf_out.new_record(
                contig=chrom,
                start=pos - 1,
                alleles=(ref, alt)
            )

            # Copy INFO fields
            for key in old_record.info:
                try:
                    new_rec.info[key] = old_record.info[key]
                except Exception:
                    pass  # ignore fields incompatible with copied header

            # Add recurrence count + file list
            new_rec.info["RECURRENCE"] = info["count"]
            new_rec.info["FILES"] = list(info["files"])

            # Copy QUAL and FILTER
            new_rec.qual = old_record.qual

            # Add filters safely: if filter not present in output header, try to add a minimal header entry;
            # otherwise skip it (so we won't raise KeyError).
            for f in old_record.filter.keys():
                try:
                    new_rec.filter.add(f)
                except KeyError:
                    # Try to add a minimal FILTER definition to header, then add to record.
                    # Minimal description - adjust if you want more accurate descriptions.
                    try:
                        header.filters.add(f, ".", "String", "Copied filter (added dynamically)")
                        new_rec.filter.add(f)
                    except Exception:
                        # If adding to header fails, just skip that filter
                        continue

            # Remove sample data — set dummy genotype to missing
            new_rec.samples["DUMMY"]["GT"] = (None, None)

            vcf_out.write(new_rec)

    # -----------------------------------
    # Write report
    # -----------------------------------
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(report, "w") as f:
        f.write(f"Recurrent Variants Report: {timestamp}\n")
        f.write(f"Output VCF: {output_file}\n")
        f.write(f"Minimum recurrence: {min_recurrence}\n\n")

        for (chrom, pos, ref, alt), info in variant_data.items():
            if info["count"] >= min_recurrence:
                f.write(f"{chrom}:{pos} {ref}>{alt}\t{info['count']} samples\tFiles: {', '.join(info['files'])}\n")


create_recurrent_vcf(input_files, output_file, min_recurrence, report)
