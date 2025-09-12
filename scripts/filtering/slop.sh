#!/usr/bin/env bash
set -euo pipefail

# ==== USER VARIABLES ====
LCRBED="/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/LCR-hs38.bed.gz"
SRBED="/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/GRCh38.UCSC.SimpleRepeats.simple.bed"
GENOME="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/genome.chrom.sizes"                     # chrom sizes file
VCFDIR="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/GB115"                           # directory of .vcf.gz files
DISTANCES=(100 300 500)
BEDS=("LCR" "SR")                          # flank sizes to check

# ==== PREP OUTPUT DIRS ====
for d in "${DISTANCES[@]}"; do
    mkdir -p results_${d}bp
done

# ==== MAKE LCR FLANKS ====
for d in "${DISTANCES[@]}"; do
    bedtools flank -i "$REFBED" -g "$GENOME" -l "$d" -r "$d" > flanks_${d}bp.bed
done

# ==== MAKE LCR FLANKS ====
for d in "${DISTANCES[@]}"; do
    bedtools flank -i "$REFBED" -g "$GENOME" -l "$d" -r "$d" > flanks_${d}bp.bed
done

# ==== LOOP THROUGH VCFs ====
for vcf in "${VCFDIR}"/*snvs.vcf.gz; do
    sample=$(basename "$vcf" .vcf.gz)

    # Convert VCF to BED (0-based start, end = POS if END missing)
    tmpbed=$(mktemp)
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$vcf" \
        | awk -v OFS="\t" '{ if ($3=="") $3=$2+1; print }' > "$tmpbed"

    for d in "${DISTANCES[@]}"; do
        # Find closest variant to each flank interval
        bedtools closest -a {bed}_flanks_${d}bp.bed -b "$tmpbed" -d \
        | awk -v maxd="$d" '$NF <= maxd {print}' \
        > results_${d}bp/${sample}_closest_${d}bp.tsv
    done

    rm "$tmpbed"
done