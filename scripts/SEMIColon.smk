# /scratch/ucgd/lustre-work/quinlan/data-shared/datasets/spermseq/wgs

# https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html unrecommended mamba install

# config file at some point

# Pipeline plans
# - Align fastqs to reference genome
#     - Lets do whatever Lee-Six did, to hg38
# - Variant calling (Call mutations)
#     - FreeBayes?
#     - DeepVariant?
#     - Likely SNVs, maybe indels?
# - Pipeline to process called variants/Vcf
#     - PASS variants, filter
# - Mutational signatures
#     - https://github.com/HLee-Six/colon_microbiopsies what Lee-Six did (in R)
# - Mutational spectra
# - Make plots and save to rational places


rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"