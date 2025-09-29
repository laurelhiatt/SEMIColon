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

with open('/uufs/chpc.utah.edu/common/HIPAA/u1264408/lustre/u1264408/MuSiCal/19610DeNovoSig_screen.pkl', 'rb') as f:
    model = pickle.load(f)

print(model.n_components)

model.plot_selection()

fig = musical.sigplot_bar(model.W)