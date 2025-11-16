ref="/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta"

for bam in *bam; do
    echo "Processing $bam"
    # Extract the sample name by removing the -sorted.bam suffix
    sample=$(basename "$bam" -sorted.bam)
    # Run inSTRbility
    python /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/inSTRbility/inSTRbility/core.py -ref $ref -bam $bam -bed /uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/inSTRbility/inSTRbility/MSI-catalog.bed --min-reads 8 --threads 2 --reads-out -o $sample.tsv
done

