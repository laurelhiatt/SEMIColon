import musical
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-pcklein', '-p', help='pickle input, the output from de novo signature discovery', required=True)
parser.add_argument('-H_Sout', '-H', help='H_S output, signature assignments', required=True)
parser.add_argument('-W_Sout', '-W', help='W_S output, signatures', required=True)
# parser.add_argument('--clean_Ws', '-C', choices=(True, False),  help= "Whether or not we are using the clean_W_s option to avoid overfitting", type = bool)
parser.add_argument('--clean', action='store_true')
parser.add_argument('--no-clean', dest='clean', action='store_false')
parser.set_defaults(feature=True)

args = parser.parse_args()
pcklein = args.pcklein
H_Sout = args.H_Sout
W_Sout = args.W_Sout
clean_Ws = args.clean

# print(clean_Ws)
# optional, used when I was checking to make sure that the Boolean argument was functioning

with open(pcklein, 'rb') as f:
    model = pickle.load(f)


thresh_grid = np.array([
    0.0001, 0.0002, 0.0005,
    0.001, 0.002, 0.005,
    0.01, 0.02, 0.05,
    0.1, 0.2, 0.5,
    1., 2., 5.
])

catalog = musical.load_catalog('/uufs/chpc.utah.edu/common/HIPAA/u1264408/u1264408/Git/SEMIColon/data/output/MuSiCal/musical/data/COSMIC_v3p2_SBS_WGS.csv')
W_catalog = catalog.W
print(W_catalog.shape[1])
# when/if we decide to work with restricted signatures, we want to make sure we have the right number being input as the catalog

model.assign_grid(W_catalog,
                  method_assign='likelihood_bidirectional', # Method for performing matching and refitting
                  thresh_match_grid=thresh_grid, # Grid of threshold for matchinng
                  thresh_refit_grid=thresh_grid, # Grid of threshold for refitting
                  thresh_new_sig=0.0, # De novo signatures with reconstructed cosine similarity below this threshold will be considered novel
                  connected_sigs=False, # Whether or not to force connected signatures to co-occur
                  clean_W_s= clean_Ws # clean_W_s=False # An optional intermediate step to avoid overfitting to small backgrounds in de novo signatures for 96-channel SBS signatures
                 )



model.validate_grid(validate_n_replicates=1, # Number of simulation replicates to perform for each grid point
                    grid_selection_method='pvalue', # Method for selecting the best grid point
                    grid_selection_pvalue_thresh=0.05 # Threshold used for selecting the best grid point
                   )

print(model.best_grid_point)
print(model.thresh_match)
print(model.thresh_refit)

W_s = model.W_s
H_s = model.H_s

H_s.to_csv(H_Sout)
W_s.to_csv(W_Sout)

print(W_s.columns.tolist())