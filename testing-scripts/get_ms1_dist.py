# python script that gest the total number of ms1 peaks, the total number of ms1 spectra, and then prints thes
# with the file name

import numpy as np
import pandas as pd
import random
import scipy.stats as stats
from pyopenms import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('mzML_file')
args = parser.parse_args()

# load th ms file
def get_number_ms1_peaks(mzml_file):
# load the ms file
    exp = MSExperiment()
    MzMLFile().load(mzml_file, exp)

# for every spectra, if it's MSLevel is 1, then append to the vector spec
    spec = []
    for s in exp.getSpectra():
        if s.getMSLevel() == 1:
             spec.append(s)

# declare three variables
    total_ms1_peaks = 0
    mz_list = []
    i_list = []

# for each spectra in spec, colelct the mz and intensity values, the number of ms1_peaks (add those to 
# total_ms1_peaks), and then append the mz list onto the mz_list object
    for ind_spec in spec:
         mz, i = ind_spec.get_peaks()
         ms1_peaks = len(mz)
         total_ms1_peaks += ms1_peaks
         mz_list.append(mz)

    print(mzml_file, total_ms1_peaks, len(mz_list))

get_number_ms1_peaks(mzml_file = args.mzML_file)
