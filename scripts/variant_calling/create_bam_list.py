# create_bam_list.py
import json
import sys
import yaml
import os

# Get the absolute path to config.yaml
config_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../data/config/snakemake_config", "CCconfig.yaml"))

# Load the YAML file
with open(config_path, "r") as file:
    config = yaml.safe_load(file)

# Access matches
matches = config["matches"]

print(matches)  # Debugging: See the loaded dictionary


donor = sys.argv[1]
location = sys.argv[2]
out_dir = sys.argv[3]

crypt_samples = matches[donor]["crypt_samples"]

with open(f"{out_dir}/bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
