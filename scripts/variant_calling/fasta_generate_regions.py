#!/usr/bin/env python3
import os
from math import ceil

chroms = snakemake.input["chroms"]

chunks = int(snakemake.params["chunks"])
index = snakemake.params["index"]
out_dir = snakemake.params["out_dir"]

directory_path = os.path.join(out_dir, "regions")
os.makedirs(directory_path, exist_ok=True)

bed_prefix = os.path.join(directory_path, f"chunk.{chroms}")

def generate_chunked_regions(fai_path, num_chunks, chromosome, bed_prefix):
    if not fai_path.endswith(".fai"):
        fai_path += ".fai"

    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])

            if chrom_name != chromosome:
                continue

            chunk_size = ceil(chrom_length / num_chunks)
            for i in range(num_chunks):
                start = i * chunk_size
                end = min(start + chunk_size, chrom_length)
                bed_file = f"{bed_prefix}.region.{i+1}.bed"
                with open(bed_file, "w") as bed:
                    bed.write(f"{chrom_name}\t{start}\t{end}\n")

# Run
generate_chunked_regions(index, chunks, chroms, bed_prefix)