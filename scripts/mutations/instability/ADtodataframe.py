import pandas as pd
import glob
import os

indir = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/bam"
outfile = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/scripts/mutations/instability.tsv"
catalog_file = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/inSTRbility/inSTRbility/MSI-catalog.bed"

# Read catalog file (no header)
catalog_cols = ["CHROM", "START", "END", "MOTIF", "NAME"]
catalog = pd.read_csv(catalog_file, sep="\t", names=catalog_cols)

# Dictionary to store all data per sample
sample_dict = {}

# Iterate over catalog
for idx, row in catalog.iterrows():
    chrom = str(row["CHROM"]).replace("chr", "")
    start = int(row["START"])
    end = int(row["END"])
    name = row["NAME"]

    for tsv_file in glob.glob(os.path.join(indir, "*.tsv")):
        tsv_df = pd.read_csv(tsv_file, sep="\t")
        if tsv_df.shape[0] <= 1:
            continue
        tsv_df.columns = tsv_df.columns.str.strip()

        locus_split = tsv_df["#locus_id"].str.strip().str.split("[:-]", expand=True)
        tsv_df["CHROM"] = locus_split[0].str.replace("chr", "")
        tsv_df["LOC_START"] = locus_split[1].astype(int)
        tsv_df["LOC_END"] = locus_split[2].astype(int)

        tsv_chr = tsv_df[tsv_df["CHROM"] == chrom]
        overlap = tsv_chr[(tsv_chr["LOC_START"] <= end+20) & (tsv_chr["LOC_END"] >= start-20)]

        if overlap.empty:
            continue

        sample_name = os.path.basename(tsv_file).replace(".tsv", "")

        mean_ad_list = [f"{o_row['mean_ad']}({o_row['haplogroup']})" for _, o_row in overlap.iterrows()]
        median_ad_list = [f"{o_row['median_ad']}({o_row['haplogroup']})" for _, o_row in overlap.iterrows()]

        mean_ad_str = ", ".join(mean_ad_list)
        median_ad_str = ", ".join(median_ad_list)

        # Initialize sample entry if it doesn't exist
        if sample_name not in sample_dict:
            sample_dict[sample_name] = {}

        # Store columns by locus name
        sample_dict[sample_name][f"{name}_mean_ad"] = mean_ad_str
        sample_dict[sample_name][f"{name}_median_ad"] = median_ad_str

# Convert to DataFrame (one row per sample)
output_df = pd.DataFrame.from_dict(sample_dict, orient="index")
output_df.index.name = "Sample"
output_df.reset_index(inplace=True)

# Save to file
output_df.to_csv(outfile, sep="\t", index=False)