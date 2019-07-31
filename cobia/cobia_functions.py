# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 10:04:01 2017

These are functions for calculating cofragmentation scores (as in cobia_cofrag.py).

_filter_ion_parcel() is used to select an ion parcel from a given retention time bin.
cofrag_ion_counter_targeted() is used for targeted calculations of cofragmentation scores.
cofrag_ion_counter_sparse_para() is used for global calculations of cofragmention scores, and therefore
uses sparse sampling and parallel computing option. This function requires cobia_parallel.py

@author: J. Scott P. McCain
"""
import pandas as pd
import numbers
import numpy as np

# helper function to filter a peptide dataframe by an injection bin
def _filter_ion_parcel(peptide_df, filtered_injection_bin_i):
    if type(filtered_injection_bin_i) != pd._libs.interval.Interval:
        raise NameError('Error with retention time binning. Check pandas version for interval class support (>0.2.0).')
    ion_parcel = peptide_df[(peptide_df.rt_upper > filtered_injection_bin_i.left) &
                        (peptide_df.rt_upper < filtered_injection_bin_i.right) |
                        (peptide_df.rt_lower > filtered_injection_bin_i.left) &
                        (peptide_df.rt_lower < filtered_injection_bin_i.right) |
                        (peptide_df.rt_lower <= filtered_injection_bin_i.left) &
                        (peptide_df.rt_upper >= filtered_injection_bin_i.right)]
    return(ion_parcel)

# targeted function just for a list of peptides:
def cofrag_ion_counter_targeted(max_injection_time, precursor_selection_window, list_of_targets, peptide_df):
    # Checking variable types:
    if isinstance(precursor_selection_window, numbers.Number) != True:
        raise NameError('precursor_selection_window must be a number.')
    # Initializing returned object
    cofrag_ion_series = pd.Series(name = 'cofrag_ions')
    ## go through each target peptide
    for peptide in list_of_targets.itertuples(index = True, name = "Pandas"):
        # subset the target peptide
        temp_pep = getattr(peptide, "pep_seq")
        print(temp_pep)
        # subset the entire potential metaproteome (input fasta file w/ lc predictions).
        temp_pep_df = peptide_df[peptide_df['peptide_sequence'].str.contains(temp_pep + '-')]
        if len(temp_pep_df) > 1:
            print(temp_pep_df)
            temp_pep_df = temp_pep_df.loc[temp_pep_df['mz'] == temp_pep_df['mz'].min()]
            print('its working')
#            raise NameError('Eek! this is something Ive been meaning to fix. The peptide matches more than one potential peptide, likely because its a subsequence. Please annoy me to figure it out')
        if not len(temp_pep_df) == 1:
            raise NameError('Eek! One of the target peptides you supplied is not in the database of potential peptides, check to see you would actually expect it:' + temp_pep)
        list(temp_pep_df)
        # find the first time the peptide enters the MS, and the last.
        pep_rt_lower = float(temp_pep_df['rt_lower'])
        pep_rt_upper = float(temp_pep_df['rt_upper'])
        pep_mz = float(temp_pep_df['mz'])
        # subset all co-eluting peptides into a new df.
        sub_pep_df = peptide_df[(peptide_df.rt_upper > pep_rt_lower) &
                                (peptide_df.rt_lower < pep_rt_lower) |
                                (peptide_df.rt_upper > pep_rt_upper) &
                                (peptide_df.rt_lower < pep_rt_upper) |
                                (peptide_df.rt_lower <= pep_rt_lower) &
                                (peptide_df.rt_upper >= pep_rt_upper)]
        # subset all injection bins found in this subset.
        sub_pep_df_unique_injection_bins = sub_pep_df['injection_bins'].unique()
        # low hanging fruit: this loops through injection bins that don't contain the peptide (cofrag score = 0).
        # for each injection bin, count the number of cofragmenting ions the target peptide has.
        for injection_bin_i in range(0, len(sub_pep_df_unique_injection_bins)):
            filtered_injection_bin = sub_pep_df_unique_injection_bins[injection_bin_i]
            # Ion parcel includes any ion which has an overlapping retention time bin (referred to with
            # rt_upper/rt_lower columns) within the injection_bin.
            ion_parcel = _filter_ion_parcel(peptide_df = sub_pep_df,
                                        filtered_injection_bin_i = filtered_injection_bin)
            cofrag_ion_count = ion_parcel.mz[(ion_parcel.mz < pep_mz + precursor_selection_window/2) &
                                             (ion_parcel.mz > pep_mz - precursor_selection_window/2)].shape[0]
            ion_cofrag = pd.Series([cofrag_ion_count], index = [peptide.Index], name = 'cofrag_count')
            cofrag_ion_series = cofrag_ion_series.append(ion_cofrag)
    return(cofrag_ion_series)

def cofrag_ion_counter_sparse_para(para_input_no, max_injection_time,
#                             sorted_unique_injection_bins,
                              peptide_unique_dict,
                              precursor_selection_window,
                              return_dict,
                              every_nth):

    #Go through and do a manual de-bug
    #para_input_no = 2
    #peptide_unique_dict = para_test_df
    #precursor_selection_window = 3
    #every_nth = 4
    #
    #print('starting off')
    # Checking variable types:
    if isinstance(precursor_selection_window, numbers.Number) != True:
        raise NameError('precursor_selection_window must be a number.')
    if isinstance(every_nth, (int, np.integer)) != True:
        raise NameError('every_nth must be an integer.')
#    if type(sorted_unique_injection_bins) != list:
#        raise NameError('input must be a sorted list')
    if type(peptide_unique_dict) != dict:
        raise NameError('function designed for a pd.DataFrame')

    parallel_zerod = para_input_no - 1

    cofrag_ion_series = pd.Series(name = 'cofrag_ions')

    try:
        peptide_unique = peptide_unique_dict[parallel_zerod][1]
    except:
        print('Incorrect number of parallel chunks')

    #print(peptide_unique_dict_sub)
    #
    #print(peptide_unique)

    subset_injection_bins_sorted = sorted(peptide_unique_dict[parallel_zerod][0])
    #print(subset_injection_bins_sorted)
    #subset_injection_bins_sorted = sorted(peptide_unique['injection_bins'].unique())
    #print(type(subset_injection_bins_sorted))
    #print(uniform_injection_bins[0:5])

    # If you just take the injection bins from the observations, you can get gaps in injection bins. This recreates a uniform injection bin
    # across the entire retention time window. It's a bit slower, but you don't miss any peptides.
    uniform_injection_bin_ranges = np.arange(start = subset_injection_bins_sorted[0].left, stop = subset_injection_bins_sorted[-1].right, step = max_injection_time)
    temp_series = np.arange(start = subset_injection_bins_sorted[0].left, stop = subset_injection_bins_sorted[-1].right, step = max_injection_time)
    uniform_injection_bins = pd.cut(x = temp_series.tolist(), bins = uniform_injection_bin_ranges, right = True, include_lowest = True)

    #print(uniform_injection_bins[0:5])
    # It's redundant to go into every bin, given the ratio of ion_peak_width and
    # max_injection_bin. So we remove every_nth bin to subsample to approximate
    # the actual number of cofragmenting ions.
    #if every_nth <= 1:
    #    subset_injection_bins_sorted_sparse = subset_injection_bins_sorted
    #else:
    #    subset_injection_bins_sorted_sparse = subset_injection_bins_sorted[0::every_nth]

    if every_nth <= 1:
        subset_injection_bins_sorted_sparse = uniform_injection_bins
    else:
        subset_injection_bins_sorted_sparse = uniform_injection_bins[0::every_nth]

    # Add the end of the list to get another observation, in case of a non-divisible
    # number of bins.
    #if not subset_injection_bins_sorted_sparse[-1] == uniform_injection_bins[-1]:
    #    print('im here')
    #    print(uniform_injection_bins[-1])
    #    print(type(subset_injection_bins_sorted_sparse))
    #    subset_injection_bins_sorted_sparse.add_categories(uniform_injection_bins[-1])
        #
        #subset_injection_bins_sorted_sparse.append(uniform_injection_bins[-1])
    #if parallel_zerod == 0:
    #	peptide_unique.to_csv('tester.csv')
    #print(subset_injection_bins_sorted[0:29])
    #print(subset_injection_bins_sorted_sparse)
    #counter_i = 1
    for injection_bin_i in range(0, len(set(subset_injection_bins_sorted_sparse))):

        filtered_injection_bin = subset_injection_bins_sorted_sparse[injection_bin_i]

        # Ion parcel includes any ion which has an overlapping retention time bin (referred to with
        # rt_upper/rt_lower columns) within the injection_bin.
        ion_parcel = _filter_ion_parcel(peptide_df = peptide_unique,
                                        filtered_injection_bin_i = filtered_injection_bin)
	#if parallel_zerod == 4:
	    #print(ion_parcel.shape)
            #if counter_i == 1:
                #print(ion_parcel)
                #print(ion_parcel[(ion_parcel.mz > 330) & (ion_parcel.mz < 350)])
                #counter_i += 1
            #print(filtered_injection_bin)
        print('injection bin i:', injection_bin_i)
        for ion in ion_parcel.itertuples():
            cofrag_ion_count = ion_parcel.mz[(ion_parcel.mz < ion.mz + precursor_selection_window/2) &
                                             (ion_parcel.mz > ion.mz - precursor_selection_window/2)].shape[0]
            ion_cofrag = pd.Series([cofrag_ion_count], index = [ion.Index], name = 'cofrag_count')
            cofrag_ion_series = cofrag_ion_series.append(ion_cofrag)

    return_dict[parallel_zerod] = cofrag_ion_series

    return(return_dict)
#%%
if __name__ == "__main__":
    _filter_ion_parcel()
    cofrag_ion_counter_sparse_para()
    cofrag_ion_counter_targeted()
