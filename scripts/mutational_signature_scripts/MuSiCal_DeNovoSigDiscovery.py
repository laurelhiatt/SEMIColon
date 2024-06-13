#import sys
#print(sys.path)
#'/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/tools/ENTER/envs/python37_musical/lib/python3.7
# 

import argparse 
parser = argparse.ArgumentParser()
parser.add_argument('-matrix', '-m', help='The Mutation Count Matrix for the data set', required=True)
parser.add_argument('-pckle', '-p', help='pickle output, please use data name and date', required=True)
args = parser.parse_args()
matrix = args.matrix
pckle = args.pckle 

#conda activate python37_musical
import musical
print('hello')

#%matplotlib inline
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle
#%load_ext autoreload
#%autoreload 2

x = pd.read_csv(matrix, index_col=0)
# print(x.head())

model = musical.DenovoSig(x, 
                          min_n_components=1, # Minimum number of signatures to test
                          max_n_components=20, # Maximum number of signatures to test
                          init='random', # Initialization method
                          method='mvnmf', # mvnmf or nmf
                          n_replicates=20, # Number of mvnmf/nmf replicates to run per n_components
                          ncpu=10, # Number of CPUs to use
                          max_iter=100000, # Maximum number of iterations for each mvnmf/nmf run
                          bootstrap=True, # Whether or not to bootstrap X for each run
                          tol=1e-8, # Tolerance for claiming convergence of mvnmf/nmf
                          verbose=1, # Verbosity of output
                          normalize_X=False # Whether or not to L1 normalize each sample in X before mvnmf/nmf
                         )
model.fit()

with open(pckle, 'wb') as f:
    pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

print('goodbye')