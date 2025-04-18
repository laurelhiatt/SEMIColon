import sys
from cyvcf2 import VCF

def count_unique_snvs(vcf_path):
    vcf = VCF(vcf_path)
    sample_names = vcf.samples
    unique_snvs = {sample: 0 for sample in sample_names}

    for variant in vcf:
        if not variant.is_snp:
            continue  # Skip non-SNVs

        genotypes = variant.gt_types  # 0=hom-ref, 1=het, 2=unknown, 3=hom-alt

        # Ensure no sample has an unknown genotype (2)
        if 2 in genotypes:
            continue

        # Identify samples with non-ref alleles (het or hom-alt)
        non_ref_samples = [i for i, gt in enumerate(genotypes) if gt in {1, 3}]

        # Count as unique only if exactly one sample is non-ref
        if len(non_ref_samples) == 1:
            unique_snvs[sample_names[non_ref_samples[0]]] += 1

    for sample, count in unique_snvs.items():
        print(f"{sample}: {count} unique SNVs")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input.vcf.gz>")
        sys.exit(1)
    count_unique_snvs(sys.argv[1])
