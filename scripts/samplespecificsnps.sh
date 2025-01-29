#!/bin/bash

### Sample 6
# Input VCF file
VCF_FILE="data/output/CellCut/vcf/toedit_GB115-annotated-var.vcf"

# Output file for SNPs unique to GB115_Laurel-6
OUTPUT_FILE="data/output/CellCut/vcf/GB115_Laurel-6_unique_snps.vcf"

# Extract SNPs unique to GB115_Laurel-6
awk '
    BEGIN {OFS="\t"}
    /^#/ {print $0; next}  # Print header lines as-is
    {
        # Ensure TYPE=snp in the INFO field
        if ($8 ~ /TYPE=snp/) {
            # Extract the genotype for GB115_Laurel-6 (column 16 in this case)
            laurel6 = $16

            # Skip if Laurel-6 has no evidence of the alternate allele
            if (laurel6 ~ /^0\/0/ || laurel6 == "./.") next

            # Check the genotypes of all other samples
            unique = 1
            for (i = 10; i <= NF; i++) {
                if (i != 16 && $i == laurel6) {
                    unique = 0
                    break
                }
            }

            # Print the line if Laurel-6 has a unique genotype
            if (unique) print $0
        }
    }
' "$VCF_FILE" > "$OUTPUT_FILE"

echo "Unique SNPs for GB115_Laurel-6 have been written to $OUTPUT_FILE"


#!/bin/bash

###Sample 11
# Input VCF file
VCF_FILE="data/output/CellCut/vcf/toedit_GB115-annotated-var.vcf"

# Output file for SNPs unique to GB115_Laurel-6
OUTPUT_FILE="data/output/CellCut/vcf/GB115_Laurel-11_unique_snps.vcf"

# Extract SNPs unique to GB115_Laurel-6
awk '
    BEGIN {OFS="\t"}
    /^#/ {print $0; next}  # Print header lines as-is
    {
        # Ensure TYPE=snp in the INFO field
        if ($8 ~ /TYPE=snp/) {
            # Extract the genotype for GB115_Laurel-11 (column 14 in this case)
            laurel11 = $14

            # Skip if Laurel-11 has no evidence of the alternate allele
            if (laurel11 ~ /^0\/0/ || laurel11 == "./.") next

            # Check the genotypes of all other samples
            unique = 1
            for (i = 10; i <= NF; i++) {
                if (i != 14 && $i == laurel11) {
                    unique = 0
                    break
                }
            }

            # Print the line if Laurel-11 has a unique genotype
            if (unique) print $0
        }
    }
' "$VCF_FILE" > "$OUTPUT_FILE"

echo "Unique SNPs for GB115_Laurel-11 have been written to $OUTPUT_FILE"
