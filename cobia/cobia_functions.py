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
        # subset the entire potential metaproteome (input fasta file w/ lc predictions).
        temp_pep_df = peptide_df[peptide_df['peptide_sequence'].str.contains(temp_pep + '-')]
        if len(temp_pep_df) > 1:
            print(temp_pep_df)
            temp_pep_df = temp_pep_df.loc[temp_pep_df['mz'] == temp_pep_df['mz'].min()]
        if not len(temp_pep_df) == 1:
            print('Eek! One of the target peptides you supplied is not in the database of potential peptides, check to see you would actually expect it:' + temp_pep)
            # if this condition is met, move to the next peptide
            continue
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

    # Checking variable types:
    if isinstance(precursor_selection_window, numbers.Number) != True:
        raise NameError('precursor_selection_window must be a number.')
    if isinstance(every_nth, (int, np.integer)) != True:
        raise NameError('every_nth must be an integer.')
    if type(peptide_unique_dict) != dict:
        raise NameError('function designed for a pd.DataFrame')

    parallel_zerod = para_input_no - 1

    cofrag_ion_series = pd.Series(name = 'cofrag_ions')

    try:
        peptide_unique = peptide_unique_dict[parallel_zerod][1]
    except:
        print('Incorrect number of parallel chunks')

    subset_injection_bins_sorted = sorted(peptide_unique_dict[parallel_zerod][0])

    # If you just take the injection bins from the observations, you can get gaps in injection bins. This recreates a uniform injection bin
    # across the entire retention time window. It's a bit slower, but you don't miss any peptides.
    uniform_injection_bin_ranges = np.arange(start = subset_injection_bins_sorted[0].left, stop = subset_injection_bins_sorted[-1].right, step = max_injection_time)
    temp_series = np.arange(start = subset_injection_bins_sorted[0].left, stop = subset_injection_bins_sorted[-1].right, step = max_injection_time)
    uniform_injection_bins = pd.cut(x = temp_series.tolist(), bins = uniform_injection_bin_ranges, right = True, include_lowest = True)

    # It's redundant to go into every bin, given the ratio of ion_peak_width and
    # max_injection_bin. So we remove every_nth bin to subsample to approximate
    # the actual number of cofragmenting ions.
    if every_nth <= 1:
        subset_injection_bins_sorted_sparse = uniform_injection_bins
    else:
        subset_injection_bins_sorted_sparse = uniform_injection_bins[0::every_nth]

    for injection_bin_i in range(0, len(set(subset_injection_bins_sorted_sparse))):

        filtered_injection_bin = subset_injection_bins_sorted_sparse[injection_bin_i]

        # Ion parcel includes any ion which has an overlapping retention time bin (referred to with
        # rt_upper/rt_lower columns) within the injection_bin.
        ion_parcel = _filter_ion_parcel(peptide_df = peptide_unique,
                                        filtered_injection_bin_i = filtered_injection_bin)
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
