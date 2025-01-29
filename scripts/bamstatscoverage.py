# for stats_file in /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam/*-sorted.stats; do
#     grep "^COV" "$stats_file" | sed 's/\[[^]]*\]//g' | cut -f 2- > "${stats_file%.stats}_coverage_distribution.txt"
# done


import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Directory containing the coverage distribution files
directory = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam"

# Process each coverage distribution file
for file in glob.glob(os.path.join(directory, "*_coverage_distribution.txt")):
    print(f"Processing file: {file}")
    coverage_data = []

    # Load coverage values
    with open(file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) == 2:  # Ensure the line has expected columns
                coverage = int(fields[0])  # Extract the coverage value
                count = int(fields[1])     # Extract the count
                coverage_data.extend([coverage] * count)  # Repeat `coverage` value `count` times

    # Check if coverage data is populated
    if not coverage_data:
        print(f"No data found in {file}, skipping...")
        continue

    # Compute statistics
    mean_coverage = np.mean(coverage_data)
    median_coverage = np.median(coverage_data)
    max_coverage = np.max(coverage_data)

    # Print statistics
    print(f"{os.path.basename(file)}: Mean={mean_coverage}, Median={median_coverage}, Max={max_coverage}")

    # Generate histogram
    plt.hist(coverage_data, bins=50, color='blue', alpha=0.7, log=True)
    plt.title(f"Coverage Distribution for {os.path.basename(file)}")
    plt.xlabel("Coverage")
    plt.ylabel("Frequency (log scale)")

    # Save histogram as PNG
    plt.savefig(file.replace("_coverage_distribution.txt", "_coverage_histogram.png"))
    plt.clf()  # Clear the figure for the next plot
