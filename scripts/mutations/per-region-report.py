import sys
from cyvcf2 import VCF

# Define regions and their corresponding samples
REGIONS = {
    "ascending_colon": ["GB115_Laurel-1", "GB115_Laurel-2", "GB115_Laurel-3"],
    "cecum": ["GB115_Laurel-5", "GB115_Laurel-6", "GB115_Laurel-7", "GB115_Laurel-8", "GB115_Laurel-10", "GB115_Laurel-11", "GB115_Laurel-12"]
}

def count_unique_snvs_per_region(vcf_path):
    vcf = VCF(vcf_path, strict_gt=True)
    sample_names = vcf.samples

    # Convert region definitions to sample indices
    region_samples = {region: [sample_names.index(s) for s in samples if s in sample_names] for region, samples in REGIONS.items()}

    unique_snvs_per_region = {region: 0 for region in REGIONS}

    for variant in vcf:
        if not variant.is_snp:
            continue  # Skip non-SNVs

        genotypes = variant.gt_types  # 0=hom-ref, 1=het, 2=unknown, 3=hom-alt

        # Ensure no sample has an unknown genotype (2)
        if 2 in genotypes:
            continue

        region_variant_presence = {region: [i for i in indices if genotypes[i] in {1, 3}] for region, indices in region_samples.items()}

        # Count as unique only if exactly one region has non-ref alleles
        non_ref_regions = [region for region, samples in region_variant_presence.items() if len(samples) > 0]

        if len(non_ref_regions) == 1:
            unique_snvs_per_region[non_ref_regions[0]] += 1

    for region, count in unique_snvs_per_region.items():
        print(f"{region}: {count} unique SNVs")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python per-region-report.py <input.vcf.gz>")
        sys.exit(1)

    count_unique_snvs_per_region(sys.argv[1])
