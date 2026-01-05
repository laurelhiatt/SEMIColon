from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze
genInstall.install('GRCh38')
import pandas as pd
import numpy as np
from scipy.optimize import nnls
from pathlib import Path


#input_data = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/vaf_spectra/snvs/output/SBS/snvs.SBS96.all"
input_matrix_csv = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/Hiatt_continueif2unknownDec82025/snvs/output/SBS/snvs.SBS96.all"
# input_matrix_csv = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/CellCut/results/Hiatt_continueif2unknownDec82025/indels/output/ID/indels.ID83.all"
input_data = input_matrix_csv
#output_dir = "snvs"
output_dir = "snvs_clonal"
output = output_dir
colibactin_sig_csv = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/scripts/mutations/signatures/SBS88.tsv"
colibactin_signature = "/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/scripts/mutations/signatures/snvs_clonal/Colibactin_plus_UNDEFINED.tsv"

# # load data
# counts = pd.read_csv(input_matrix_csv, index_col=0, sep = '\t')   # index = 96 contexts, cols = samples
# sig_colib = pd.read_csv(colibactin_sig_csv, index_col=0, sep = '\t')  # index = same 96 contexts

# # ensure contexts align
# if not counts.index.equals(sig_colib.index):
#     raise SystemExit("Row (context) labels do not match between input matrix and signature file. Align/transpose/collapse first.")

# # signature matrix (96 x 1)
# S = sig_colib.values.astype(float)  # shape (96,1)

# # compute residuals per sample
# residuals_accum = np.zeros(S.shape[0], dtype=float)

# reconstructions = {}
# exposures = {}

# for sample in counts.columns:
#     y = counts[sample].values.astype(float)   # shape (96,)

#     # fit non-negative least squares to single signature
#     x, rnorm = nnls(S, y)   # x is exposure(s), here length 1
#     recon = S.dot(x)        # reconstructed counts (96,)
#     resid = y - recon
#     resid_clipped = np.clip(resid, 0, None)   # negative residuals -> 0

#     residuals_accum += resid_clipped
#     reconstructions[sample] = recon
#     exposures[sample] = x

# # create UNDEFINED signature from accumulated residuals
# undefined_sig = residuals_accum.copy()
# if undefined_sig.sum() == 0:
#     print("WARNING: residuals sum to zero. The provided signature explains all counts (within NNLS). No meaningful UNDEFINED signature.")
# else:
#     undefined_sig /= undefined_sig.sum()   # normalize to sum = 1 (probability profile)

# # save a two-column signature file: Colibactin + UNDEFINED
# sig_df = pd.DataFrame({
#     sig_colib.columns[0]: sig_colib.iloc[:,0],
#     "UNDEFINED": undefined_sig
# }, index=counts.index)

# sig_out = output_dir + "/Colibactin_plus_UNDEFINED.csv"
# sig_df.to_csv(sig_out)
# print("Saved new signature file to:", sig_out)

# # Also optionally save per-sample exposures & reconstructions for inspection
# exp_df = pd.DataFrame({s: exposures[s][0] for s in exposures}, index=[sig_colib.columns[0]])
# exp_df.loc["UNDEFINED"] = np.nan  # UNDEFINED exposures are unknown until you re-run fit with two signatures
# exp_df.to_csv(output_dir + "/initial_single_sig_exposures.csv")


Analyze.cosmic_fit(input_data, output, input_type="matrix", context_type="96", cosmic_version=3.4, exome=False,
                   genome_build="GRCh38", signature_database=colibactin_signature,collapse_to_SBS96 = True,
                   exclude_signature_subgroups=None, export_probabilities=True,
                   export_probabilities_per_mutation=False, make_plots=False,
                   sample_reconstruction_plots=False, verbose=True)
