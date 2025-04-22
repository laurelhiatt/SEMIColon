# Laurel Hiatt
# Last validated: 04/14/2025

sys.stdout = sys.stderr = open(snakemake.log.stdio, "w")

# create_bam_list.py
import json
import sys
import yaml
import os

# Access matches
matches = snakemake.params.matches

out_dir = snakemake.params.out_dir

donor= snakemake.input.donor

location = out_dir + "/bam/"

crypt_samples = matches[donor]["crypt_samples"]

with open(location + "bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
