import numpy as np
import pandas as pd
import random
import scipy.stats as stats
from pyopenms import *

# load th ms file
def get_ms1_data(mzml_file):

    exp = MSExperiment()
    MzMLFile().load(mzml_file, exp)

    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() == 1:
             spec.append(s)

    mz_list = []
    i_list = []
    for ind_spec in spec:
         mz, i = ind_spec.get_peaks()
         mz_list.append(mz)
         i_list.append(i)

    int_list = np.concatenate(i_list).ravel()
    mz_list_c = np.concatenate(mz_list).ravel()

    mz_series = pd.Series(mz_list_c)
    intensity_series = pd.Series(int_list)
    finale_dataframe = pd.concat([mz_series, intensity_series], axis = 1)
    return finale_dataframe

run1_kleiner_full = get_ms1_data(mzml_file = "../data/kleiner_data/Run1_C4_832ng_PP.mzML")
s03_f30_full = get_ms1_data(mzml_file = "../../ross-sea-meta-omics/data/mzML-converted/180412_0749_097_S03_Rep_1_PP.mzML")

run1_kleiner = run1_kleiner_full.sample(500000)
s03_f30 = s03_f30_full.sample(500000)

run1_kleiner.to_csv("../data/ms1_dist_run1_kleiner.csv")
s03_f30.to_csv("../data/ms1_dist_s03_rep1.csv")

