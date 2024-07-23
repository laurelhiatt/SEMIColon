import os
import subprocess
import glob

def run_kraken2(input_fastq1, input_fastq2, db_path, output_report, output_classification):
    command = [
        "kraken2",
        "--db", db_path,
        "--report", output_report,
        "--output", output_classification,
        "--paired", input_fastq1, input_fastq2
    ]

    # Use Popen to capture stdout and stderr
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Decode output and error to strings
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')

    # Print stdout and stderr
    print("STDOUT:\n", stdout)
    print("STDERR:\n", stderr)

    # Check if the command was successful
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, command)

def find_paired_fastq_files(directory):
    fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz"))
    paired_files = {}

    for file in fastq_files:
        base_name = os.path.basename(file)
        print(base_name)
        if "R1" in base_name:
            pair = base_name.replace("R1", "R2")
            if os.path.join(directory, pair) in fastq_files:
                paired_files[file] = os.path.join(directory, pair)

    return paired_files

def main():
    #directory = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/cryptdata/fastq/"
    directory = "/scratch/ucgd/lustre-labs/quinlan/data-shared/SEMIColon/24286R/24286R/Fastq/shortname_fastqs/"
    print(directory)
    db_path = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/standard_db/"
    print(db_path)

    paired_files = find_paired_fastq_files(directory)
    print(paired_files)

    for input_fastq1, input_fastq2 in paired_files.items():
        base_name = os.path.basename(input_fastq1).replace("R1", "")
        output_report = "{}_kraken2_report.txt".format(base_name)
        output_classification = "{}_kraken2_classification.txt".format(base_name)

        run_kraken2(input_fastq1, input_fastq2, db_path, output_report, output_classification)

if __name__ == "__main__":
    main()
