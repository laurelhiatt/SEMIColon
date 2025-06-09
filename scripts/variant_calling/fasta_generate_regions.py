#!/usr/bin/env python3
import os
from math import ceil

chroms = snakemake.input["chroms"]

index = snakemake.params["index"]

out_dir = snakemake.params["out_dir"]


directory_path = os.path.join(out_dir, "regions")
os.makedirs(directory_path, exist_ok=True)

bed_prefix = os.path.join(directory_path, f"chunk.{chroms}")

def generate_chunked_regions(fai_path, chromosome, bed_prefix):
    if not fai_path.endswith(".fai"):
        fai_path += ".fai"

    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            if chrom_name != chromosome:
                continue
            chrom_length = int(fields[1])
            print(f"Processing {chrom_name} with length {chrom_length}")


            chunk_num = ceil(chrom_length / 5000000)

            for i in range(chunk_num):
                start = i * 5000000
                end = min(start + 5000000, chrom_length)
                bed_file = f"{bed_prefix}.region.{i+1}.bed"
                with open(bed_file, "w") as bed:
                    bed.write(f"{chrom_name}\t{start}\t{end}\n")

# Run
generate_chunked_regions(index, chroms, bed_prefix)