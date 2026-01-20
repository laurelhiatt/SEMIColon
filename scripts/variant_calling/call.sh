#!/bin/bash
set -e

module load singularity

export SINGULARITYENV_TMPDIR=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/GEMS/SEMIColon/deep_variant_tmp/

singularity exec --cleanenv \
        -H $SINGULARITYENV_TMPDIR \
        -B /usr/lib/locale/:/usr/lib/locale/ \
            ${snakemake_input[sif]} \
            run_deepsomatic \
                --model_type ${snakemake_params[model_type]} \
                --num_shards ${snakemake[threads]} \
                ${snakemake_params[normal_cmd]} \
                --reads_tumor ${snakemake_input[tumor]} \
                --output_vcf ${snakemake_output[vcf]} \
                --output_gvcf ${snakemake_output[gvcf]} \
                ${snakemake_params[normal_name]} \
                --sample_name_tumor ${snakemake_wildcards[crypt]} \
                --ref ${snakemake_input[ref]} \
                --regions ${snakemake_wildcards[chrom]} \
