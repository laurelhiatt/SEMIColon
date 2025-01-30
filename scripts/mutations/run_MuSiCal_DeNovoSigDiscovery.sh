#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -c 10
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=laurel.hiatt@hsc.utah.edu

python /uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/MuSiCal_DeNovoSigDiscovery.py 
