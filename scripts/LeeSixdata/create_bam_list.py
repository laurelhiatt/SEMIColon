# create_bam_list.py
import json
import sys

donor = sys.argv[1]
location = sys.argv[2]
out_dir = sys.argv[3]
matches = {
    'HLS': {'crypt_samples': ['HLS_1C_30_B5', 'HLS_1C_30_D5', 'HLS_1C_30_G5', 'HLS_1C_30_H5', 'HLS_2C_30_D6', 'HLS_2C_30_E6']},
    'PD28690': {'crypt_samples': ['PD28690bx_2_a5', 'PD28690bx_2_d5', 'PD28690cb_2_g3', 'PD28690cb_2_h3']},
    'PD34199': {'crypt_samples': ['PD34199a_c17', 'PD34199a_c22', 'PD34199a_c27', 'PD34199a_c28', 'PD34199a_c31', 'PD34199a_c3', 'PD34199a_c5', 'PD34199a_c7', 'PD34199a_t31', 'PD34199a_t28']},
    'PD34202': {'crypt_samples': ['PD34202a_s9']}
    }

crypt_samples = matches[donor]["crypt_samples"]

with open(f"{out_dir}/bam_list_{donor}.txt", "w") as f:
    for sample in crypt_samples:
        f.write(f"{location}{sample}-sorted.bam\n")
