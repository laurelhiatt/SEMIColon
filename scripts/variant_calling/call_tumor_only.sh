#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/GEMS/SEMIColon/deep_variant_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deepsomatic \
                --model_type WGS_TUMOR_ONLY \
                --num_shards ${snakemake[threads]} \
                --reads_tumor ${snakemake_input[tumor]} \
                --output_vcf ${snakemake_output[vcf]} \
                --output_gvcf ${snakemake_output[gvcf]} \
                --sample_name_tumor ${snakemake_wildcards[CRYPT]} \
                --ref ${snakemake_input[ref]} \
                --regions ${snakemake_wildcards[CHROM]} \
