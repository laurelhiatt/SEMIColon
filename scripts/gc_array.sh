#!/bin/bash
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-rw
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --array=1-417%40
#SBATCH --output=logs/gc_%A_%a.out
#SBATCH --error=logs/gc_%A_%a.err


set -euo pipefail

#source /uufs/chpc.utah.edu/common/HIPAA/u1264408/software/pkg/miniforge3/etc/profile.d/conda.sh
#conda activate gc_callable
PYTHON=/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/tools/software/pkg/miniforge3/envs/gc_callable/bin/python
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" depth_files.txt)

BASENAME=$(basename $FILE .per-base.bed.gz)
OUTFILE=/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/scripts/results/${BASENAME}.gc.txt
REFERENCE=/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta


if [[ -f "$OUTFILE" ]]; then
    echo "Output exists for $BASENAME, skipping."
    exit 0
fi

$PYTHON gc_callable.py \
    --depth_file $FILE --reference $REFERENCE --output $OUTFILE