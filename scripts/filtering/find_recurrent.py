import pysam
import os
from collections import defaultdict
from datetime import datetime

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


input_files = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/GB130/24_130DC_filtered.vcf.gz,/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/GB74/D27_074_SI_filtered.vcf.gz".split(",")
output_file = "recurrent.vcf"
min_recurrence= 2
report = "report.txt"

# input_files = snakemake.input
# output_file = snakemake.output["recurrent_vcf"]
# min_recurrence= snakemake.params["recurrence_threshold"]
# report = snakemake.output["report"]

def create_recurrent_vcf(input_files, output_file, min_recurrence, report):
    """
    Creates a VCF file containing only recurrent mutations from a cohort.
    """
    if not input_files:
        print("No input VCF files provided.")
        return

    # Use the header from the first VCF file as a template
    first_vcf = pysam.VariantFile(input_files[0])
    header = first_vcf.header.copy()

    # Add a custom INFO field to the header for "Recurrence Count"
    header.add_line('##INFO=<ID=RECURRENCE,Number=1,Type=Integer,Description="Number of samples this variant was found in.">')

    variant_counts = defaultdict(lambda: [0, None])


    for vcf_file in input_files:
        vcf_in = pysam.VariantFile(vcf_file)
        for record in vcf_in:
            # Create a unique key for the variant (handling multiple alts)
            for alt in record.alts:
                variant_key = (record.chrom, record.pos, record.ref, alt)
                if variant_counts[variant_key][0] == 0:
                    # Store the first occurrence's record object
                    variant_counts[variant_key][1] = record
                variant_counts[variant_key][0] += 1
        vcf_in.close()

    # Write recurrent variants to the output VCF
    with pysam.VariantFile(output_file, 'w', header=header) as vcf_out:
        for (chrom, pos, ref, alt), (count, old_record) in variant_counts.items():

            if count < min_recurrence:
                continue

        # Create new record tied to output header
            new_rec = vcf_out.new_record(
                contig=chrom,
                start=pos - 1,
                stop=pos - 1 + len(ref),
                alleles=(ref, alt))

        # Copy INFO fields
            for key in old_record.info:
                if key in vcf_out.header.info:       # safety check
                    new_rec.info[key] = old_record.info[key]

        # Add RECURRENCE
            new_rec.info["RECURRENCE"] = count

        # Copy QUAL
            new_rec.qual = old_record.qual

        # Copy FILTER
            for f in old_record.filter.keys():
                new_rec.filter.add(f)

        # Write record
            vcf_out.write(new_rec)



    with open(report, "w") as f:
        f.write(f"Recurrent Variants Report for {timestamp} \n")
        f.write(f"Recurrent variants VCF created at: {output_file}\n")
        f.write(f"Minimum recurrence threshold: {min_recurrence}\n")
        f.write(f"Total unique variants found: {len(variant_counts)}\n")
        f.write(f"Variants meeting min_recurrence ({min_recurrence}): {sum(1 for c, r in variant_counts.values() if c >= min_recurrence)}\n")

# Create the recurrent VCF
create_recurrent_vcf(input_files, output_file, min_recurrence, report)