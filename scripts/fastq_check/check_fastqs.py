import os
import gzip

# Directory containing the FASTQ files
input_directory='/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/spermseq/wgs/19841R/Fastq/merged_trimmed_fastq'
output_directory='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/'

print("Starting the script...", flush=True)
print(f"Input directory: {input_directory}", flush=True)
print(f"Output directory: {output_directory}", flush=True)



def fix_fastq_format(input_file, output_file):
    with gzip.open(input_file, 'rt') as infile, gzip.open(output_file, 'wt') as outfile:
        line_num = 0
        for line in infile:
            line_num += 1
            # Check for the sequence identifier line
            if line_num % 4 == 1:
                if not line.startswith('@'):
                    print(f"Error: Expected '@' at line {line_num} in {input_file}, got {line.strip()}", flush=True)
                    # Fix the error if necessary
                    line = '@' + line[1:]
            # Check for the plus line
            elif line_num % 4 == 3:
                if not line.startswith('+'):
                    print(f"Error: Expected '+' at line {line_num} in {input_file}, got {line.strip()}", flush=True)
                    # Fix the error if necessary
                    line = '+\n'
            outfile.write(line)

def process_directory(directory, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(directory):
        if filename.endswith(".fastq.gz"):
            input_file = os.path.join(directory, filename)
            output_file = os.path.join(output_directory, f"fixed_{filename}")
            fix_fastq_format(input_file, output_file)
            print(f"Processed {input_file}, fixed file saved as {output_file}", flush=True)


# Process all FASTQ files in the directory
process_directory(input_directory, output_directory)
