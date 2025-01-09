# create_bam_list.py
import json
import sys

donor = sys.argv[1]
location = sys.argv[2]
out_dir = sys.argv[3]
matches = {
    'GB115': {'crypt_samples': ['Laurel-1', 'Laurel-2', "Laurel-3", "Laurel-5", "Laurel-6", "Laurel-7",
    "Laurel-8", "Laurel-10", "Laurel-11", "Laurel-12", "Laurel-13", "Laurel-14", "Laurel-15", "Laurel-16", "Laurel-17"]},
}
crypt_samples = matches[donor]["crypt_samples"]

with open(f"{out_dir}/bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
