#!/bin/bash
sites='/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/sites.hg38.vcf.gz'
fasta='/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta'
# file_path='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/bam'
out_dir='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut'

file_path='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/vcf'

# samples=(
#     "D091_IX8Y012"
#     "D100_VG6U010"
#     "D127_V5ET001"
#     "D130_WN11-001"
#     "D212_92N8014"
# )

donors=(
    "GB115"
)

# for sample in "${samples[@]}"; do
#     file="${file_path}/${sample}-sorted.bam"
#     echo $file
#     /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier extract $file --sites $sites --fasta $fasta -d $out_dir
# done


for donor in "${donors[@]}"; do
    file="${file_path}/${donor}-var.vcf.gz"
    echo $file
    /uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier extract $file --sites $sites --fasta $fasta -d $out_dir
done