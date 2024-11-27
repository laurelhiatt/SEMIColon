#!/bin/bash

# Define input and output directories
input_dir="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/SlideScrape/bam"
output_dir="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/SlideScrape/mosdepth"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all sorted BAM files in the input directory
for bam_file in "$input_dir"/*-sorted.bam; do
  # Extract the base filename without extension
  base_name=$(basename "$bam_file" -sorted.bam)

  # Run mosdepth with full output
  mosdepth -n --fast-mode --by 500 \
           "$output_dir/$base_name" \
           "$bam_file"
done

echo "Mosdepth processing completed. Outputs are in $output_dir"
