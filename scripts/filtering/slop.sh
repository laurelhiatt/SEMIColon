#!/usr/bin/env bash
set -euo pipefail

LCRBED="/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/LCR-hs38.bed.gz"
SRBED="/scratch/ucgd/lustre-labs/quinlan/data-shared/annotations/GRCh38.UCSC.SimpleRepeats.simple.bed"
GENOME="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/genome.chrom.sizes"                     # chrom sizes file
VCFDIR="/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/vaf_spectra_under885"                           # directory of .vcf.gz files
DISTANCES=(100 300 500)                       # flank sizes to check

for d in "${DISTANCES[@]}"; do
    mkdir -p results_${d}bp
done

# LCR slops
for d in "${DISTANCES[@]}"; do
    bedtools slop -i "$LCRBED" -g "$GENOME" -b "$d" > LCR_slops_${d}bp.bed
done

# simple repeat slops
for d in "${DISTANCES[@]}"; do
    bedtools slop -i "$SRBED" -g "$GENOME" -b "$d" > SR_slops_${d}bp.bed
done

# loo pthrough vcf
for vcf in "${VCFDIR}"/*snvs.vcf; do
    sample=$(basename "$vcf" .vcf)

    # Convert VCF to BED (0-based start, end = POS if END missing)
    tmpbed=$(mktemp)
    bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' "$vcf" \
        | awk -v OFS="\t" '{ if ($3=="") $3=$2+1; print }' > "$tmpbed"

    for d in "${DISTANCES[@]}"; do
        # Find closest variant to each flank interval
        bedtools intersect -a LCR_slops_${d}bp.bed -b "$tmpbed" \
        | awk -v maxd="$d" '$NF <= maxd {print}' \
        > results_${d}bp/${sample}_LCRoverlap_${d}bp.tsv

        bedtools intersect -a SR_slops_${d}bp.bed -b "$tmpbed" \
        | awk -v maxd="$d" '$NF <= maxd {print}' \
        > results_${d}bp/${sample}_SRoverlap_${d}bp.tsv
    done

    rm "$tmpbed"
done