# SEMIColon
Scripts and data for SEMIColon: Somatic Exploration of Mosaicism in Colon

# Contents

We have a scripts folder and a data folder.

Within the scripts folder, there are (planned) subfolders specific to different processes

# Cleaning Branch To-Do's
- [X] Containerize software
  - [X] Create conda environments in container
  - [X] Create environments in container for respective envmodules
- [ ] Hardcode paths
  - [ ] Remove hardcoded paths
  - [ ] Replace hardcoded paths with relative paths, input parameters, etc.
- [ ] Replace redundant copy+paste rules with reused rules using `use rule ... as ... with:`
- [ ] Determine need of symlinks for inSTRbility, MuSiCal.
- [ ] Determine how `figures`, `filtering`, `mutations`, and `pks` relate to the pipeline.
  - [ ] Include `filtering` and `pks` Snakefiles as modules in main Snakefile?
- [ ] Get Laurel's feedback on which scripts to remove and which scripts to keep.


## Directions
Please review the following files and directories. To mark those that you want removed (or those that should not show up in the final manuscript repository), please fill in the brackets with an X (i.e. [X]). If you mark a directory with this, I will assume that the directory and all contents within the directory are fine to remove. If you do this, don't worry about marking each independent file in that directory. 

| Folder / File | Remove |
| ---           | ---           |
| ├── data  | [ ] |
| │   ├── clb_genes.gtf  | [X] |
| │   ├── config  | [ ] |
| │   │   ├── sample-groups.lua  | [ ] |
| │   │   ├── slurm_config  | [ ] |
| │   │   │   └── config.yaml  | [ ] |
| │   │   ├── snakemake_config  | [ ] |
| │   │   │   ├── CCconfig.yaml  | [ ] |
| │   │   │   └── toyconfig.yaml  | [X] |
| │   │   ├── spectra_config  | [X] |
| │   │   │   ├── config.yaml  | [X] |
| │   │   │   └── spectraconfig.yaml  | [X] |
| │   │   └── test_config  | [X] |
| │   │       ├── config.yaml  | [X] |
| │   │       └── testconfig.yaml  | [X] |
| │   ├── COSMIC_v3.4_ID_GRCh37.txt  | [ ] |
| │   ├── COSMIC_v3.4_SBS_GRCh38.txt  | [X] |
| │   ├── COSMIC_v3.5_SBS_GRCh38.txt  | [ ] |
| │   ├── GCF_000025745.1_ASM2574v1_genomic.fna  | [X] |
| │   ├── genomic.gtf  | [X] |
| │   ├── refcds_GRCh38_hg38.rda  | [X] |
| │   └── somalier-groups.txt  | [X] |
| ├── envs  | [X] |
| │   ├── alfred2.yaml  | [X] |
| │   ├── alfred.yaml  | [X] |
| │   ├── cyvcf2.yaml  | [X] |
| │   ├── fastp.yaml  | [X] |
| │   ├── make_bams.yaml  | [X] |
| │   ├── plot_mosdepth.yaml  | [X] |
| │   ├── recurrent.yaml  | [X] |
| │   └── vcfstats.yaml  | [X] |
| ├── genome.chrom.sizes  | [X] |
| ├── inSTRbility  | [X] |
| ├── MuSiCal  | [X] |
| ├── README.  | [ ] |
| └── scripts  | [ ] |
|     ├── figures  | [ ] |
|     │   ├── archive  | [X] |
|     │   │   ├── Figure1.ipynb  | [ ] |
|     │   │   ├── Figure2.ipynb  | [ ] |
|     │   │   ├── Figure3.ipynb  | [ ] |
|     │   │   ├── Figure4.ipynb  | [ ] |
|     │   │   ├── other_sigs.ipynb  | [ ] |
|     │   │   ├── Statistics.ipynb  | [ ] |
|     │   │   ├── subclonal_sigs.ipynb  | [ ] |
|     │   │   └── Supp.ipynb  | [ ] |
|     │   ├── Figure1  | [ ] |
|     │   │   ├── base_quality.svg  | [ ] |
|     │   │   ├── coverage.svg  | [ ] |
|     │   │   ├── Hiatt_SBS96_profile.svg  | [ ] |
|     │   │   ├── indels_accumulation.svg  | [ ] |
|     │   │   ├── SNVs_accumulation.svg  | [ ] |
|     │   │   └── VAF_Hiatt_median_publication_stair.svg  | [ ] |
|     │   ├── Figure1.ipynb  | [ ] |
|     │   ├── Figure2  | [ ] |
|     │   │   ├── exp_v_con_mean.svg  | [ ] |
|     │   │   ├── exp_v_con_sum.svg  | [ ] |
|     │   │   ├── Indels_accumulation_side.svg  | [ ] |
|     │   │   ├── msi_by_region.svg  | [ ] |
|     │   │   ├── mut_prop_side.svg  | [ ] |
|     │   │   └── SNVs_accumulation_side.svg  | [ ] |
|     │   ├── Figure2.ipynb  | [ ] |
|     │   ├── Figure3  | [ ] |
|     │   │   ├── cohortsigs_SBS.svg  | [ ] |
|     │   │   ├── cornish_SBS18.svg  | [ ] |
|     │   │   ├── cornish_SBS1_subclonal.svg  | [ ] |
|     │   │   ├── cornish_SBS1.svg  | [ ] |
|     │   │   ├── cornish_SBS5_subclonal.svg  | [ ] |
|     │   │   ├── cornish_SBS5.svg  | [ ] |
|     │   │   ├── cornish_SNVs_subclonal.svg  | [ ] |
|     │   │   ├── cornish_SNVs.svg  | [ ] |
|     │   │   ├── SBS18_clonal.svg  | [ ] |
|     │   │   ├── SBS1_clonal.svg  | [ ] |
|     │   │   ├── SBS5_clonal.svg  | [ ] |
|     │   │   └── unique_SNVs_clonal.svg  | [ ] |
|     │   ├── Figure3_4.ipynb  | [ ] |
|     │   ├── Figure4  | [ ] |
|     │   │   ├── cohortsigs_ID.svg  | [ ] |
|     │   │   ├── cornish_ID18.svg  | [ ] |
|     │   │   ├── cornish_ID1.svg  | [ ] |
|     │   │   ├── cornish_ID2_subclonal.svg  | [ ] |
|     │   │   ├── cornish_ID2.svg  | [ ] |
|     │   │   ├── cornish_INDELs_subclonal.svg  | [ ] |
|     │   │   ├── cornish_INDELs.svg  | [ ] |
|     │   │   ├── ID18_clonal.svg  | [ ] |
|     │   │   ├── ID1_clonal.svg  | [ ] |
|     │   │   ├── ID2_clonal.svg  | [ ] |
|     │   │   └── total_indels_clonal.svg  | [ ] |
|     │   └── ProcessData.ipynb  | [ ] |
|     ├── filtering  | [ ] |
|     │   ├── filter.ipynb  | [X] |
|     │   ├── find_recurrent.py  | [ ] |
|     │   ├── recurrent.ipynb  | [X] |
|     │   ├── slop.sh  | [X] |
|     │   └── Snakefile  | [X] |
|     ├── mutations  | [ ] |
|     │   ├── dataframe.ipynb  | [ ] |
|     │   ├── gene-specific.ipynb  | [X] |
|     │   ├── instability  | [ ] |
|     │   │   ├── ADtodataframe.py  | [ ] |
|     │   │   ├── instability.tsv  | [X] |
|     │   │   └── inSTRbility.sh  | [ ] |
|     │   ├── matrix.ipynb  | [ ] |
|     │   ├── per-sample-report-indels.py  | [ ] |
|     │   ├── per-sample-report.py  | [ ] |
|     │   ├── selection.ipynb  | [X] |
|     │   ├── signatures  | [ ] |
|     │   │   ├── LookingatMuSiCalModelPlots.ipynb  | [X] |
|     │   │   ├── MuSiCal_DeNovoSigAnalysis.py  | [X] |
|     │   │   ├── MuSiCal_DeNovoSigDiscovery.py  | [X] |
|     │   │   ├── MuSiCal_DeNovoSigMatching.py  | [X] |
|     │   │   ├── MuSiCal_SigMatchingOnly.py  | [X] |
|     │   │   ├── mutation_signatures.Rmd  | [X] |
|     │   │   ├── SBS88.tsv  | [X] |
|     │   │   ├── sigassignment.py  | [X] |
|     │   │   ├── sigpro_dbs.py  | [ ] |
|     │   │   ├── SigProfiler.py  | [X] |
|     │   │   ├── sigpro_IDs.py  | [ ] |
|     │   │   └── sigpro.py  | [ ] |
|     │   ├── split_vcfs_by_vaf.py  | [X] |
|     │   └── within-donor  | [ ] |
|     │       └── donor_shared.py  | [X] |
|     ├── overlap.ipynb  | [X] |
|     ├── pks  | [ ] |
|     │   ├── 2026_01_15_indel_loads.csv  | [X] |
|     │   ├── 2026_01_26_ID18_refit_abs.csv  | [X] |
|     │   ├── 2026_01_26_ID18_refit.csv  | [X] |
|     │   ├── 2026_01_26_indel_loads.csv  | [X] |
|     │   ├── 2026_01_26_SBS88_refit_abs.csv  | [X] |
|     │   ├── 2026_01_26_SBS88_refit.csv  | [X] |
|     │   ├── 2026_01_26_snv_loads.csv  | [ X |
|     │   ├── colibactin_all.ipynb  | [ ] |
|     │   ├── colibactin_coding_genic.csv  | [X] |
|     │   ├── colibactin_coding.ipynb  | [ ] |
|     │   ├── colibactin_driver_coding_gene_snvs.csv  | [X] |
|     │   ├── colibactin_driver_gene_IDs.csv  | [X] |
|     │   ├── colibactin_driver_gene_snvs.csv  | [X] |
|     │   ├── colibactinFig2.ipynb  | [ ] |
|     │   ├── colibactin_id_burden_coding.csv  | [X] |
|     │   ├── colibactin_id_burden_genic.csv  | [X] |
|     │   ├── colibactin.ipynb  | [ ] |
|     │   ├── colibactin_snv_burden_coding.csv  | [X] |
|     │   ├── colibactin_snv_burden_genic.csv  | [X] |
|     │   ├── ID18_refit_coding.csv  | [X] |
|     │   ├── ID18_refit_genic.csv  | [X] |
|     │   ├── Output  | [X] |
|     │   │   └── Supplementary_tables  | [ ] |
|     │   │       └── Supplementary_Table_1.doc  | [ ] |
|     │   ├── Processed_data  | [X] |
|     │   │   └── Contexts  | [ ] |
|     │   │       ├── adenoma_unique_contexts_all.txt  | [ ] |
|     │   │       ├── adenoma_unique_contexts.txt  | [ ] |
|     │   │       ├── carcinoma_unique_contexts_all.txt  | [ ] |
|     │   │       ├── normal_unique_contexts_all.txt  | [ ] |
|     │   │       └── normal_unique_contexts.txt  | [ ] |
|     │   ├── SBS88_refit_coding.csv  | [X] |
|     │   ├── SBS88_refit_genic.csv  | [X] |
|     │   └── Snakefile  | [X] |
|     ├── pks.sh  | [X] |
|     ├── quality_control   | [ ] |
|     │   ├── plot-dist.py  | [ ] |
|     │   └── stats.R  | [ ] |
|     ├── rules  | [ ] |
|     │   ├── 0_merge_fastq.smk  | [ ] |
|     │   ├── 1_fastp.smk  | [ ] |
|     │   ├── 2_make_bams.smk  | [ ] |
|     │   ├── 3_check_bams.smk  | [ ] |
|     │   ├── 4_make_vcfs.smk  | [ ] |
|     │   ├── 5_check_vcfs.smk  | [ ] |
|     │   ├── 6_filter_vcfs.smk  | [ ] |
|     │   └── 7_mutations.smk  | [ ] |
|     ├── Snakefile  | [ ] |
|     └── variant_calling  | [ ] |
|         ├── call.sh  | [ ] |
|         ├── call_tumor_only.sh  | [ ] |
|         ├── create_bam_list.py  | [ ] |
|         ├── fasta_generate_regions.py  | [ ] |
|         ├── GenerateFreebayesRegions.R  | [X] |
|         └── joint.sh  | [ ] |
