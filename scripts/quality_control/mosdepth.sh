#!/bin/bash

set -euo pipefail

out_dir="${snakemake_params[out_dir]}"
BAM_FILE="${snakemake_input[bam_sort]}"
log="${snakemake_log[0]}"

# Define input and output directories
input_dir="$out_dir/bam"
output_dir="$out_dir/mosdepth"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo "Processing $BAM_FILE in $input_dir" > "$log"
echo "Output will be saved to $output_dir" >> "$log"

base_name=$(basename "$BAM_FILE" -sorted.bam)

# Run mosdepth with full output
mosdepth -n --fast-mode --by 500 \
           "$output_dir/$base_name" \
           "$BAM_FILE" 2>> "$log"

echo "Mosdepth processing completed. Outputs are in $output_dir" >> "$log"