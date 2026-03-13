# SEMIColon
Scripts and data for SEMIColon: Somatic Exploration of Mosaicism in Colon

# Contents

We have a scripts folder and a data folder.

Within the scripts folder, there are (planned) subfolders specific to different processes

# Cleaning Branch To-Do's
- [ ] Containerize software
  - [ ] Create conda environments in container
  - [ ] Create environments in container for respective envmodules
- [ ] Hardcode paths
  - [ ] Remove hardcoded paths
  - [ ] Replace hardcoded paths with relative paths, input parameters, etc.
- [ ] Replace redundant copy+paste rules with reused rules using `use rule ... as ... with:`
- [ ] Determine need of symlinks for inSTRbility, MuSiCal.
- [ ] Determine how `figures`, `filtering`, `mutations`, and `pks` relate to the pipeline.
  - [ ] Include `filtering` and `pks` Snakefiles as modules in main Snakefile?
- [ ] Get Laurel's feedback on which scripts to remove and which scripts to keep.


| Folder / File | Remove |
| ---           | ---           |
| ├── data  | [ ] |
| │   ├── clb_genes.gtf  | [ ] |
| │   ├── config  | [ ] |
| │   │   ├── sample-groups.lua  | [ ] |
| │   │   ├── slurm_config  | [ ] |
| │   │   │   └── config.yaml  | [ ] |
| │   │   ├── snakemake_config  | [ ] |
| │   │   │   ├── CCconfig.yaml  | [ ] |
| │   │   │   └── toyconfig.yaml  | [ ] |
| │   │   ├── spectra_config  | [ ] |
| │   │   │   ├── config.yaml  | [ ] |
| │   │   │   └── spectraconfig.yaml  | [ ] |
| │   │   └── test_config  | [ ] |
| │   │       ├── config.yaml  | [ ] |
| │   │       └── testconfig.yaml  | [ ] |
| │   ├── COSMIC_v3.4_ID_GRCh37.txt  | [ ] |
| │   ├── COSMIC_v3.4_SBS_GRCh38.txt  | [ ] |
| │   ├── COSMIC_v3.5_SBS_GRCh38.txt  | [ ] |
| │   ├── GCF_000025745.1_ASM2574v1_genomic.fna  | [ ] |
| │   ├── genomic.gtf  | [ ] |
| │   ├── refcds_GRCh38_hg38.rda  | [ ] |
| │   └── somalier-groups.txt  | [ ] |
| ├── envs  | [X] |
| │   ├── alfred2.yaml  | [X] |
| │   ├── alfred.yaml  | [X] |
| │   ├── cyvcf2.yaml  | [X] |
| │   ├── fastp.yaml  | [X] |
| │   ├── make_bams.yaml  | [X] |
| │   ├── plot_mosdepth.yaml  | [X] |
| │   ├── recurrent.yaml  | [X] |
| │   └── vcfstats.yaml  | [X] |
| ├── genome.chrom.sizes  | [ ] |
| ├── inSTRbility  | [X] |
| ├── MuSiCal  | [X] |
| ├── README.  | [ ] |
| └── scripts  | [ ] |
|     ├── figures  | [ ] |
|     │   ├── archive  | [ ] |
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
|     │   ├── filter.ipynb  | [ ] |
|     │   ├── find_recurrent.py  | [ ] |
|     │   ├── recurrent.ipynb  | [ ] |
|     │   ├── slop.sh  | [ ] |
|     │   └── Snakefile  | [ ] |
|     ├── mutations  | [ ] |
|     │   ├── dataframe.ipynb  | [ ] |
|     │   ├── gene-specific.ipynb  | [ ] |
|     │   ├── instability  | [ ] |
|     │   │   ├── ADtodataframe.py  | [ ] |
|     │   │   ├── instability.tsv  | [ ] |
|     │   │   └── inSTRbility.sh  | [ ] |
|     │   ├── matrix.ipynb  | [ ] |
|     │   ├── per-sample-report-indels.py  | [ ] |
|     │   ├── per-sample-report.py  | [ ] |
|     │   ├── selection.ipynb  | [ ] |
|     │   ├── signatures  | [ ] |
|     │   │   ├── LookingatMuSiCalModelPlots.ipynb  | [ ] |
|     │   │   ├── MuSiCal_DeNovoSigAnalysis.py  | [ ] |
|     │   │   ├── MuSiCal_DeNovoSigDiscovery.py  | [ ] |
|     │   │   ├── MuSiCal_DeNovoSigMatching.py  | [ ] |
|     │   │   ├── MuSiCal_SigMatchingOnly.py  | [ ] |
|     │   │   ├── mutation_signatures.Rmd  | [ ] |
|     │   │   ├── SBS88.tsv  | [ ] |
|     │   │   ├── sigassignment.py  | [ ] |
|     │   │   ├── sigpro_dbs.py  | [ ] |
|     │   │   ├── SigProfiler.py  | [ ] |
|     │   │   ├── sigpro_IDs.py  | [ ] |
|     │   │   └── sigpro.py  | [ ] |
|     │   ├── split_vcfs_by_vaf.py  | [ ] |
|     │   └── within-donor  | [ ] |
|     │       └── donor_shared.py  | [ ] |
|     ├── overlap.ipynb  | [ ] |
|     ├── pks  | [ ] |
|     │   ├── 2026_01_15_indel_loads.csv  | [ ] |
|     │   ├── 2026_01_26_ID18_refit_abs.csv  | [ ] |
|     │   ├── 2026_01_26_ID18_refit.csv  | [ ] |
|     │   ├── 2026_01_26_indel_loads.csv  | [ ] |
|     │   ├── 2026_01_26_SBS88_refit_abs.csv  | [ ] |
|     │   ├── 2026_01_26_SBS88_refit.csv  | [ ] |
|     │   ├── 2026_01_26_snv_loads.csv  | [ ] |
|     │   ├── colibactin_all.ipynb  | [ ] |
|     │   ├── colibactin_coding_genic.csv  | [ ] |
|     │   ├── colibactin_coding.ipynb  | [ ] |
|     │   ├── colibactin_driver_coding_gene_snvs.csv  | [ ] |
|     │   ├── colibactin_driver_gene_IDs.csv  | [ ] |
|     │   ├── colibactin_driver_gene_snvs.csv  | [ ] |
|     │   ├── colibactinFig2.ipynb  | [ ] |
|     │   ├── colibactin_id_burden_coding.csv  | [ ] |
|     │   ├── colibactin_id_burden_genic.csv  | [ ] |
|     │   ├── colibactin.ipynb  | [ ] |
|     │   ├── colibactin_snv_burden_coding.csv  | [ ] |
|     │   ├── colibactin_snv_burden_genic.csv  | [ ] |
|     │   ├── ID18_refit_coding.csv  | [ ] |
|     │   ├── ID18_refit_genic.csv  | [ ] |
|     │   ├── Output  | [ ] |
|     │   │   └── Supplementary_tables  | [ ] |
|     │   │       └── Supplementary_Table_1.doc  | [ ] |
|     │   ├── Processed_data  | [ ] |
|     │   │   └── Contexts  | [ ] |
|     │   │       ├── adenoma_unique_contexts_all.txt  | [ ] |
|     │   │       ├── adenoma_unique_contexts.txt  | [ ] |
|     │   │       ├── carcinoma_unique_contexts_all.txt  | [ ] |
|     │   │       ├── normal_unique_contexts_all.txt  | [ ] |
|     │   │       └── normal_unique_contexts.txt  | [ ] |
|     │   ├── SBS88_refit_coding.csv  | [ ] |
|     │   ├── SBS88_refit_genic.csv  | [ ] |
|     │   └── Snakefile  | [ ] |
|     ├── pks.sh  | [ ] |
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
|         ├── GenerateFreebayesRegions.R  | [ ] |
|         └── joint.sh  | [ ] |