#!/usr/bin/env python
#usage: python return_homopol_denovos.py -i <input_vcf> -k <kid_id> -m <mom_id> -d <dad_id> > <output_file>

import argparse
from cyvcf2 import VCF
import pysam


def get_variant_type(v):
    """
    Returns true if the variant is a hompolymer,
    which HipSTR annotates as PERIOD = 1.
    """
    if variant.INFO.get('PERIOD') == 1:
        return "homopolymer"
    else:
        return "non_homopolymer"


def get_sample_gt_column(samples, desired_sample):
    """
    Return the 0-based column of the desired sample name
    in the list of sample GT information.
    For example, the first sample in the VCF genotype column
    would have position 0 in the genotype arrays returned by cyvcf2
    """
    return samples.index(sample)


def is_gt_unknown(variant, mom, dad, kid, samples):
    """
    Return true if genotype of mom, dad, or kid is unknown
    Return false if genotypes are all known (0, 1, or 3)
    gt_types = 2 = unknown
    param samples is needed for cyvcf2
    """
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)
    gts = variant.gt_types
    # 0 = HOMREF, 1 = HET, 2 = UNKNOWN, 3 = HOMALT
    if gts[mom_col] == 2 or gts[dad_col] == 2 or gts[kid_col] == 2:
        return True
    else:
        return False


def is_this_denovo(variant, mom, dad, kid, samples):
    """
    Return true if denovo and false if not denovo
    param samples is needed for cyvcf2
    """
    #for HCH: if passed in as a parameter, can use again
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)
    gts = variant.gt_types
    # 0 = HOMREF, 1 = HET, 2 = UNKNOWN, 3 = HOMALT
    #if mom and dad are HOMREF and kid is HET, this is denovo
    if gts[mom_col] == 0 and gts[dad_col] == 0 and gts[kid_col] == 1:
        return True
    else:
        return False


def get_allele_lens(variant, mom, dad, kid, samples):
    """
    Returns allele lengths for mom, dad, and kid
    param samples is needed for cyvcf2
    """

    #get index positions for mom, dad, kid
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)

    #access gt bases for mom, get length of 2 alleles
    mom_gt_bases = variant.gt_bases[mom_col].split('|')
    mom_allele1_len = len(mom_gt_bases[0])
    mom_allele2_len = len(mom_gt_bases[1])
    mom_allele_lens = str(mom_allele1_len) + "," + str(mom_allele2_len)

    #access gt bases for dad, get length of 2 alleles
    dad_gt_bases = variant.gt_bases[dad_col].split('|')
    dad_allele1_len = len(dad_gt_bases[0])
    dad_allele2_len = len(dad_gt_bases[1])
    dad_allele_lens = str(dad_allele1_len) + "," + str(dad_allele2_len)

    #access gt bases for kid, get length of 2 alleles
    kid_gt_bases = variant.gt_bases[kid_col].split('|')
    kid_allele1_len = len(kid_gt_bases[0])
    kid_allele2_len = len(kid_gt_bases[1])
    kid_allele_lens = str(kid_allele1_len) + "," + str(kid_allele2_len)

    return mom_allele_lens, dad_allele_lens, kid_allele_lens


def get_genotypes(variant, mom, dad, kid, samples):
    """
    Returns genotypes for mom, dad, and kid
    :param samples is needed for cyvcf2
    """

    #get index positions for mom, dad, kid
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)

    #access genotypes for mom, dad, kid
    mom_gt = variant.gt_types[mom_col]
    dad_gt = variant.gt_types[dad_col]
    kid_gt = variant.gt_types[kid_col]

    return mom_gt, dad_gt, kid_gt


def get_candidate_dnms(variant, mom, dad, kid, samples):
    """
    Accesses gt bases for mom, dad, kid
    For each of kid's alleles, asks if it matches parental gt
    if yes, returns "inherited"
    if no, returns "candidate dnm"
    param samples is needed for cyvcf2
    """

    #get index positions for mom, dad, kid
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)

    #access gt bases for mom, dad, kid
    mom_gt_bases = variant.gt_bases[mom_col].split('|')
    dad_gt_bases = variant.gt_bases[dad_col].split('|')
    kid_gt_bases = variant.gt_bases[kid_col].split('|')

    #make a list to append if criteria are met
    kid_dnms = []
    for allele in kid_gt_bases:
        #if allele differs from mom and differs from dad, append to dnm list
        #if one allele differs from mom/dad, dnm list will have a length of 1
        #if both alleles differ from mom and dad, dnm list will have a length of 2
        if allele not in mom_gt_bases and allele not in dad_gt_bases:
            kid_dnms.append(allele)

    if len(kid_dnms) > 0:
        return "candidate_dnm"
    else:
        return "inherited"


def is_perfect_repeat(variant):
    """
    returns True if repeat allele is perfect
    returns False is repeat allele is interrupted
    """
    unique_characters = set(variant.REF)
    if len(unique_characters) == 1:
        return True
    else:
        return False


# def get_ref_coords(variant):
#     """
#     returns genomic coordinates for a variant
#     """
#     #get reference coordinates
#     chrom = variant.CHROM
#     ref_start = variant.start
#     ref_end = ref_start + len(variant.REF)
#
#     return chrom, ref_start, ref_end


# def calculate_read_depth(bam_file, chrom, start, end):
#     bam = pysam.AlignmentFile(bam_file, "rb")
#     read_depth = 0
#     for pileupcolumn in bam.pileup(chrom, start, end):
#         read_depth += pileupcolumn.nsegments
#         bam.close()
#     return read_depth


def get_coverage_info(variant, mom, dad, kid, samples):
    """
    For each variant, gets DP for mom, dad, kid
    HipSTR FORMAT DP = "number of valid reads used for sample's genotype"
    if number of valid reads used to genotype < 10 for mom, dad, or kid,
    return "not_enough_covg"
    """

    #get index positions for mom, dad, kid
    mom_col = samples.index(mom)
    dad_col = samples.index(dad)
    kid_col = samples.index(kid)

    #get DP and PDP for mom, dad, kid
    mom_DP = variant.format('DP')[mom_col]
    dad_DP = variant.format('DP')[dad_col]
    kid_DP = variant.format('DP')[kid_col]

    if mom_DP < 12 or dad_DP < 12 or kid_DP < 12:
        return "insufficient_covg"
    else:
        return "sufficient_covg"


def get_VAF(variant, kid, samples):
    """
    For each variant, gets PDP for kid
    HipSTR FORMAT PDP = "fractional reads supporting each haploid genotype"
    PDP divided by DP gives allele fraction
    If VAF is < 0.3, returns "low_VAF"
    """
    #get index positions for kid
    kid_col = samples.index(kid)

    #get PDP and DP for kid
    kid_PDP = variant.format('PDP')[kid_col].split('|')
    kid_DP = variant.format('DP')[kid_col]

    #PDP divided by DP gives allele fraction
    #if variant allele fraction is < 0.3, return "low_VAF"
    kid_PDP_1 = kid_PDP[0]
    kid_VAF_1 = float(kid_PDP_1) / float(kid_DP)
    kid_PDP_2 = kid_PDP[1]
    kid_VAF_2 = float(kid_PDP_2) / float(kid_DP)

    kid_VAF_list = [kid_VAF_1, kid_VAF_2]

    return kid_VAF_list
#      if kid_VAFs[1] < 0.3:
#          return "low_VAF"


#at least 5 reads supporting DNM



if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input_vcf", required=True, help="HipSTR VCF file")
    ap.add_argument("-k", "--kid_id", required=True, help="sample id of proband/kid in which to call denovos")
    ap.add_argument("-m", "--mom_id", required=True, help="sample id of mom")
    ap.add_argument("-d", "--dad_id", required=True, help="sample id of dad")

    args = ap.parse_args()

    vcf = VCF(args.input_vcf)
    kid_id = args.kid_id
    mom_id = args.mom_id
    dad_id = args.dad_id

    #testing calculate_read_depth
    kid_bam = "/scratch/ucgd/lustre-work/quinlan/data-shared/datasets/elementbio/CEPH/element_sequenced/BAMs/merged/2187-E_merged_sort.bam"


    # Print header
    print("chrom\tpos\tfamily\tchild\tchild_gt\tmat_gt\tpat_gt\trepeat_unit\tlen_repeat_unit\tref_allele_len\tis_perfect\tkid_VAFs")

    for variant in vcf:
        variant_type = get_variant_type(variant)
        gt_unknown = is_gt_unknown(variant, mom_id, dad_id, kid_id, vcf.samples)
        valid_reads = get_coverage_info(variant, mom_id, dad_id, kid_id, vcf.samples)

        #impose exclusion rules
        #first, skip nonhomopolymers (will exclude the majority)
        if variant_type != "homopolymer": continue
        #skip if gt is unknown for anyone in trio
        if gt_unknown == True: continue
        #skip if number of reads used for hipSTR genotyping in mom, dad, or kid < 8
        if valid_reads != "sufficient_covg": continue
#         #if kid VAF is < 0.3, skip
#         if kid_VAF == "low_VAF": continue


        #if one of proband's alleles is not present in mom or dad, this is a candidate dnm
        #this is very liberal and candidate dnms will need additional filtering
        if get_candidate_dnms(variant, mom_id, dad_id, kid_id, vcf.samples) == "candidate_dnm":

            #get kid_VAF
            kid_VAFs = get_VAF(variant, kid_id, vcf.samples)
#             print(kid_VAF)

            #get allele lens as input for Michelle's validation script
            mat_gt = get_allele_lens(variant, mom_id, dad_id, kid_id, vcf.samples)[0]
            pat_gt = get_allele_lens(variant, mom_id, dad_id, kid_id, vcf.samples)[1]
            child_gt = get_allele_lens(variant, mom_id, dad_id, kid_id, vcf.samples)[2]

            #get the rest of the info needed as input for Michelle's validation script
            chrom = variant.CHROM
            pos = variant.POS
            family = kid_id.split("-")[0]
            child = kid_id
            repeat_unit = variant.REF[0]
            len_repeat_unit = variant.INFO.get('PERIOD')
            ref_allele_len = len(variant.REF)
            is_perfect = is_perfect_repeat(variant)


            #print in this format as input for Michelle's validation script
            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}'.format(chrom, pos, family, child, child_gt, mat_gt, pat_gt, repeat_unit, len_repeat_unit, ref_allele_len, is_perfect, kid_VAFs))


