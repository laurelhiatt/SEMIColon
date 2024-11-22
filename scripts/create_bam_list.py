# create_bam_list.py
import json
import sys

donor = sys.argv[1]
location = sys.argv[2]
out_dir = sys.argv[3]
matches = {
    'GB115': {'crypt_samples': ['A1', 'A2', "A3", "A4", "A5", "A6", "A7", "B1", "B2", "B3", "B4"]},
}

crypt_samples = matches[donor]["crypt_samples"]

with open(f"{out_dir}/bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
