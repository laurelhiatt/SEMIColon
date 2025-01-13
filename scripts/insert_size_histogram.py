import os
import glob
import numpy as np
import csv
import matplotlib.pyplot as plt
from collections import Counter


# Directory containing the insert size files
directory = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam"

# # Process each insert size file
# for file in glob.glob(os.path.join(directory, "*_insert_sizes.txt")):
#     # Load insert sizes
#     with open(file, 'r') as f:
#         sizes = [int(line.strip()) for line in f]

#     # Compute statistics
#     min_size = min(sizes)
#     max_size = max(sizes)
#     median_size = np.median(sizes)

#     # Print statistics
#     print(f"{os.path.basename(file)}: Min={min_size}, Max={max_size}, Median={median_size}")

    # # Create histogram
    # plt.hist(sizes, bins=50, color='blue', alpha=0.7)
    # plt.title(f"Insert Size Distribution for {os.path.basename(file)}")
    # plt.xlabel("Insert Size")
    # plt.ylabel("Frequency")
    # plt.xlim(0, 10000)  # Limit the x-axis to 0-10,000

    # Save histogram as PNG
    # plt.savefig(file.replace("_insert_sizes.txt", "_histogram.png"))
    # plt.clf()  # Clear the figure for the next plot

# Process each insert size file
for file in glob.glob(os.path.join(directory, "*_insert_sizes.txt")):
    # Load insert sizes
    with open(file, 'r') as f:
        sizes = [int(line.strip()) for line in f]

    # Count frequencies of each insert size
    size_counts = Counter(sizes)

    # Output file for the table
    output_file = file.replace("_insert_sizes.txt", "_size_counts.csv")

    # Write to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Insert Size", "Count"])
        for size, count in sorted(size_counts.items()):
            writer.writerow([size, count])

    print(f"Count table saved to {output_file}")
