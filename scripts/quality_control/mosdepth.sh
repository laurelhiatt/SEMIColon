#!/bin/bash

# Define input and output directories
input_dir="../../data/output/CellCut/bam"
output_dir="../../data/output/CellCut/mosdepth"
log=${snakemake_log[0]}

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo "Processing BAM files in $input_dir" > $log
echo "Output will be saved to $output_dir" >> $log
# Loop through all sorted BAM files in the input directory
for bam_file in "$input_dir"/*-sorted.bam; do
  # Extract the base filename without extension
  base_name=$(basename "$bam_file" -sorted.bam)

  # Run mosdepth with full output
  mosdepth -n --fast-mode --by 500 \
           "$output_dir/$base_name" \
           "$bam_file" 2>> "$log"
done

echo "Mosdepth processing completed. Outputs are in $output_dir" >> $log