import sys
from cyvcf2 import VCF

def count_unique_snvs(vcf_path):
    vcf = VCF(vcf_path)
    sample_names = vcf.samples
    unique_snvs = {sample: 0 for sample in sample_names}

    for variant in vcf:
        if not variant.is_snp:
            continue  # Skip non-SNVs

        genotypes = variant.gt_types  # List of genotypes per sample (0=hom-ref, 1=het, 2=hom-alt, 3=unknown)

        # Find samples that have the variant uniquely (only one sample has non-ref alleles)
        non_ref_samples = [i for i, gt in enumerate(genotypes) if gt in {1, 2}]
        if len(non_ref_samples) == 1:
            unique_snvs[sample_names[non_ref_samples[0]]] += 1

    for sample, count in unique_snvs.items():
        print(f"{sample}: {count} unique SNVs")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input.vcf.gz>")
        sys.exit(1)
    count_unique_snvs(sys.argv[1])