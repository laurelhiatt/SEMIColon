import os
import gzip

def fix_fastq_format(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        line_num = 0
        for line in infile:
            line_num += 1
            # Check for the sequence identifier line
            if line_num % 4 == 1:
                if not line.startswith('@'):
                    print(f"Error: Expected '@' at line {line_num} in {input_file}, got {line.strip()}")
                    # Fix the error if necessary
                    line = '@' + line[1:]
            # Check for the plus line
            elif line_num % 4 == 3:
                if not line.startswith('+'):
                    print(f"Error: Expected '+' at line {line_num} in {input_file}, got {line.strip()}")
                    # Fix the error if necessary
                    line = '+\n'
            outfile.write(line)

fix_fastq_format("../data/test_R1.fastq.gz", "../data/fixed_test_R1.fastq.gz")