# Laurel Hiatt
# Last validated: 04/14/2025

# create_bam_list.py
import json
import sys
import yaml
import os

# Access matches
matches = snakemake.params.matches

out_dir = snakemake.params.out_dir

donor= snakemake.wildcards.donor

location = out_dir + "/bam/"

crypt_samples = matches[donor]["crypt_samples"]

with open(location + f"bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
