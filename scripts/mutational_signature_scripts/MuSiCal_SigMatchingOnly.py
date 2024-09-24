import argparse 
parser = argparse.ArgumentParser()
parser.add_argument('-matrix', '-m', help='The Mutation Count Matrix for the data set', required=True)
args = parser.parse_args()
matrix = args.matrix

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

X = pd.read_csv(matrix, index_col=0)


catalog = musical.load_catalog('/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/MuSiCal/MuSiCal/musical/data/COSMIC-MuSiCal_v3p2_SBS_WGS.csv')
W_catalog = catalog.W
W = catalog.W
print(W_catalog.shape[1])


# We can restrict the catalog to specific tumor types if we want
# print(catalog.show_tumor_type_options().tolist())
##this will let us see the tumor type options, like ColoRec.AdenoCA and so forth

## this is how we actually restrict tumor type 
#catalog.restrict_catalog(tumor_type='Skin.Melanoma')
#print(catalog.W.shape[1])

Hnaive, model = musical.refit.refit(X, W, method='thresh_naive', thresh=0)

print(Hnaive.head())

Hnaive.to_csv('/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/MuSiCal/MuSiCal/20220629Hnaivetable.csv')

H, model = musical.refit.refit(X, W, method='likelihood_bidirectional', thresh=0.02)

print(H.head())

H.to_csv('/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/MuSiCal/MuSiCal/20220629H_0.001_table.csv')


# Associated signatures (e.g., APOBEC signatures SBS2 and SBS13) can be forced to co-occur using the option connected_sigs=True (by default it is set to False).