# create_bam_list.py
import json
import sys

donor = sys.argv[1]
location = sys.argv[2]
out_dir = sys.argv[3]
matches = {
    'HLS': {'crypt_samples': ['HLS_1C_30_B5']},
    'PD34199': {'crypt_samples': ['PD34199a_t28']},
    }

crypt_samples = matches[donor]["crypt_samples"]

with open(f"{out_dir}/bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
