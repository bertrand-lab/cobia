# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 10:39:17 2017

Input files: a csv file with mass, peptide sequence, retention times.

# The script calculated number of cofragmenting peptides. There is an option for global or targeted.

Output files: (1) cofragmentation scores associated with peptides (2) parameter settings file.

@author: J. Scott P. McCain

"""
import os
import pandas as pd
import numbers
import numpy as np
import time
import multiprocessing as mp
import csv
from sys import argv
import math

from cobia_parallel import parallel_split
from cobia_functions import cofrag_ion_counter_sparse_para
from cobia_functions import cofrag_ion_counter_targeted

def cobia_cofrag(lcfilename, ddaparams, globalcofrag, output_name, target_df):

	print('Loading data set...')
	dda_params = pd.read_csv(ddaparams)
	peptide_master = pd.read_csv(lcfilename)
	custom_name = output_name

	if 'mass' not in peptide_master.columns:
		raise NameError('Error in peptide data file, mass column not found.')
        all_required_params = ['max_injection_time',
                               'precursor_selection_window',
                               'ion_peak_width',
                               'number_of_parallel',
                               'every_nth']
        # checking that all required parameters are present:
        if sorted(all_required_params) != sorted(list(dda_params.keys())):
            raise NameError('Error in parameter input file, check for typos or missing parameter.')
	pep_list = peptide_master['peptide_sequence'].str.split('-').tolist()
        charge_state = []
        three_indicators = ['H', 'RP', 'KP']
        print(pep_list[1:100])
        # removing the -OH and just getting the peptide
	pep_list_real = []
	for element_i in pep_list:
            if isinstance(element_i, list):
                pep_list_real.append(element_i[0])
            elif math.isnan(element_i):
                pep_list_real.append('')

	# getting charge state estimates
        for element_i in pep_list_real:
            if isinstance(element_i, basestring):
                if any(x in three_indicators for x in element_i):
                    charge_state.append(3)
                else:
                    charge_state.append(2)
            elif math.isnan(element_i):
                charge_state.append('')
	print(charge_state[1:100])
	print(type(charge_state[1]))
	# converting charge state integers to floats
	charge_state_float = [float(i) for i in charge_state]
	#adjusting mz for different charge states by peptide
        peptide_master['mz'] = peptide_master['mass']/charge_state_float

	print(peptide_master['mz'].tolist()[1:100])
	print(peptide_master['mass'].tolist()[1:100])

	# Filtering out ions between 50 and 2000 mz
	peptide_mz_filter = peptide_master.query('mz > 50 & mz < 2000')
	# Removing duplicate peptides
	peptide_unique = peptide_mz_filter.drop_duplicates(['peptide_sequence'])
	peptide_sequences = peptide_unique['peptide_sequence']
	#%% MS-sim settings:
	print('Cofragmentation-prediction with the following parameters:')
	print(dda_params)
	# Binning parameter for the injection time in minutes.
	max_injection_time = dda_params['max_injection_time'][0]
	if isinstance(max_injection_time, numbers.Number) != True:
	    raise NameError('Error in parameter input file, max_injection_time takes only Numeric.')

	# Width of precursor selection window.
	precursor_selection_window = dda_params['precursor_selection_window'][0]
	if isinstance(precursor_selection_window, numbers.Number) != True:
	    raise NameError('Error in parameter input file, precursor_selection_window takes only Numeric.')

	# Average peak width (eluting from LC column)
	ion_peak_width = dda_params['ion_peak_width'][0] # e.g. 0.11 minutes
	if isinstance(ion_peak_width, numbers.Number) != True:
	    raise NameError('Error in parameter input file, ion_peak_width takes only Numeric.')

	number_of_parallel = dda_params['number_of_parallel'][0]
        if isinstance(number_of_parallel, (int, np.integer)) != True:
	    raise NameError('Must specify an integer value for number_of_parallel nodes.')

	every_nth = dda_params['every_nth'][0]
	if isinstance(every_nth, (int, np.integer)) != True:
	    raise NameError('Must specifiy an integer value for every_nth sampling.')

	#%% Modifying the DataFrame to include a range for ion peak width and the injection bins.
	print('Computing injection bins...')

	# making bins for the injection times.
	# the stop is slightly more to make the last time inclusive for the pd.cut function below.
	injection_bins_ranges = np.arange(start = peptide_unique['rts'].min(),
		                   stop = peptide_unique['rts'].max(),# + peptide_unique['rts'].max()*0.05,
	                           step = max_injection_time)

	# sometimes the injection bins miss the last little bit, so this is creating another injection bins
	# that is the maximum rention time
	if injection_bins_ranges.max() < peptide_unique['rts'].max():
	    injection_bins_ranges = np.append(injection_bins_ranges, peptide_unique['rts'].max())

	injection_bins = pd.cut(peptide_unique['rts'].tolist(),
                        injection_bins_ranges,
                        right = True,
                        include_lowest = True)

	# Adding injection bins to the peptide_unique dataframe.
	peptide_unique = peptide_unique.assign(injection_bins = injection_bins)

	#%% Range for ion peak width.
	peptide_unique = peptide_unique.assign(rt_upper = pd.Series(peptide_unique['rts'] + ion_peak_width/2),
	                                       rt_lower = pd.Series(peptide_unique['rts'] - ion_peak_width/2))

	# sorting peptide df first by retention time, before it's split up for parallel
	peptide_unique_sorted = peptide_unique.sort_values(['rts'])

	# if doing a global cofragmentation assessment

	if globalcofrag == 'global':
	    # Splitting the dataset for parallelization:
	    para_test_df = parallel_split(p_comp = number_of_parallel,
	                              peptide_df = peptide_unique_sorted)
	    # Define an output queue
	    parallel_list = list(range(1, number_of_parallel + 1))
	    print('initializing parallel...')
            # Setting up dictionary object to return results
	    manager = mp.Manager()
            return_dict = manager.dict()
	    procs = []

            print('instantiating parallel...')
	    for core_val in parallel_list:
	        print(core_val)
	        proc = mp.Process(target = cofrag_ion_counter_sparse_para,
                      args = (core_val, max_injection_time,
                              para_test_df,
                              precursor_selection_window,
                              return_dict,
                              every_nth))
	        procs.append(proc)
	        proc.start()

	    print('joining parallel')
	    for proc in procs:
	        proc.join()

	    return_dict2 = dict(return_dict)

	    cofrag_ion_series = pd.Series(name = 'cofrag_ions')

	    for dict_series in range(0, len(return_dict2)):
	        cofrag_ion_series = cofrag_ion_series.append(return_dict2[dict_series])

        # if doing a targeted cofragmentation assessment

	if globalcofrag == 'targeted':

	    target_df = pd.read_csv(target_df)

	    cofrag_ion_series = cofrag_ion_counter_targeted(max_injection_time, precursor_selection_window, list_of_targets = target_df, peptide_df = peptide_unique_sorted)

	    # Removing zeros, as this is an artifact of the targeted subsampling. 
	    cofrag_ion_series = cofrag_ion_series[cofrag_ion_series!=0]

	# converting key observations to pandas.series
	sim_series_mean = cofrag_ion_series.groupby(cofrag_ion_series.index).mean()
	sim_series_median = cofrag_ion_series.groupby(cofrag_ion_series.index).median()
	sim_series_sd = cofrag_ion_series.groupby(cofrag_ion_series.index).std()

	sim_series_mean = sim_series_mean.rename('mean_cofrag_score')
	sim_series_median = sim_series_median.rename('median_cofrag_score')
	sim_series_sd = sim_series_sd.rename('sd_cofrag_score')

	if globalcofrag == 'targeted':
	    cofrag_df = pd.concat([sim_series_mean, sim_series_median, sim_series_sd, target_df], axis = 1)
	if globalcofrag == 'global':
	    cofrag_df = pd.concat([sim_series_mean, sim_series_median, sim_series_sd, peptide_master], axis = 1)

	# Writing file names with key parameters in the actual name.
	file_name = str(custom_name) + '_mi-' + str(max_injection_time) + '_ipw-' + str(ion_peak_width) + '_para-' + str(number_of_parallel) + '_co-sim.csv'
	param_name = str(custom_name) + '_mi-' + str(max_injection_time) + '_ipw-' + str(ion_peak_width) + '_para-' + str(number_of_parallel) + '_params.csv'

	# Saving the file.
	cofrag_df.to_csv(file_name)
	with open(param_name, 'wb') as f:  # Just use 'w' mode in 3.x
	    w = csv.DictWriter(f, dda_params.keys())
	    w.writeheader()
	    w.writerow(dda_params)
