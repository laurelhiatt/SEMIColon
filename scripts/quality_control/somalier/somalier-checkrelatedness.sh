#!/bin/bash
# SAMPLES=~/somalier_files.txt
input_dir='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/somalier'
groups='/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/somalier-groups.txt'

/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier relate ${input_dir}/*.somalier -g $groups -o groups
#/uufs/chpc.utah.edu/common/HIPAA/u1264408/tools/somalier/somalier relate ${input_dir}/*.somalier  -o nogroups