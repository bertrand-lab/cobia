# -*- coding: utf-8 -*-
"""
Created on Tue Sep 05 12:14:22 2017

Functions for parallel computing option.

@author: J. Scott P. McCain
"""

import pandas as pd
import numpy as np
import numbers
import multiprocessing as mp
import warnings

#tester2 = peptide_unique
#tester3 = tester2.sort_values(by = 'rts')

#%%

def parallel_split(p_comp,
                   peptide_df):

    '''
    This function splits a dataset into subsets for parallel computing.

    Input is a dataset (peptide_df), and the number of subsets (p_comp) to
    split the data into. overlap_rt should be the ion_peak_width value.

    Returns a dictionary, with keys of the dictionary corresponding
    to numbers of p_comp. The items in the dictionary are lists, the first
    element of the list are the coordinates of the unbiased dataframe
    in terms of rt (retention time). The second item of the list
    is a pandas.DataFrame, containing the raw data.
    '''

    #p_comp = 3
    #peptide_df = tester3

    # Checking variable types:
    if isinstance(p_comp, (int, np.integer)) != True:
        raise NameError('p_comp must be an integer >= 1.')
    if type(peptide_df) != pd.core.frame.DataFrame:
        raise NameError('function designed for a pd.DataFrame')
    if p_comp == 0:
        raise NameError('p_comp value must be > 0')
    print('Number of parallel CPUS running:', p_comp)

    # Safety statement for number of cores available
    available_cpus = mp.cpu_count()
    if p_comp > available_cpus:
        warnings.warn("You don't have enough cores available, cores adjusted to (maximmum allowable - 1).")
        p_comp = available_cpus - 1

    # assigning quantiles based on center of ion_peak_width
    #df_column_added = peptide_df.assign(p_comp = pd.qcut(peptide_df['rts'],
    #                                                     p_comp))
    #print(peptide_df['injection_bins'].unique())
    parallel_groupings = np.array_split(peptide_df['injection_bins'].unique(),
                                        p_comp)
    #print(parallel_groupings[1:10])
    #print('parallel 2')
   # print(parallel_groupings[2])
    #unique_quantiles = df_column_added['p_comp'].unique()

    # dictionary of chunked datasets
    chunked_datasets = {}

    #print('line 69')
    #for i in range(0, len(unique_quantiles)):
    if p_comp > 1:
        dict_number = 0
        for i in parallel_groupings:
            #quantile_x = unique_quantiles[i]
            #print(i)
            print('in loop')
            quantile_upper_bound = max(i).right
            print(quantile_upper_bound)
            quantile_lower_bound = min(i).left
            print(quantile_lower_bound)
            # Quantiles need to be adjusted if they are the first or last one.
            # This is because the quantiles are computed from 'rts', which is the
            # center of the ion_peak_width.
            #if i == 0:
            #    quantile_lower_bound = df_column_added.rt_lower.min()
            #else:
            #    quantile_lower_bound = quantile_x.left

            # If this is the last quantile
            #if i == len(unique_quantiles) - 1:
            #    quantile_upper_bound = df_column_added.rt_upper.max()
            #else:
            #    quantile_upper_bound = quantile_x.right

            # Subset the dataframe if the rt_lower lies within the quantiles, or
            # the rt_upper lies within the quantiles.
            sub_df = peptide_df[(peptide_df.rt_lower > quantile_lower_bound) &
                                (peptide_df.rt_lower < quantile_upper_bound) |
                                (peptide_df.rt_upper > quantile_lower_bound) &
                                (peptide_df.rt_upper < quantile_upper_bound) |
                                (peptide_df.rt_lower < quantile_lower_bound) &
                                (peptide_df.rt_upper > quantile_upper_bound)]

            # Converting lower bounds, converting upper bounds
            #sub_df.loc[(sub_df.rt_lower < quantile_lower_bound),'rt_lower'] = quantile_lower_bound
            #sub_df.loc[(sub_df.rt_upper > quantile_upper_bound),'rt_upper'] = quantile_upper_bound
	    #print(sub_df)
            chunked_datasets[dict_number] = [i, sub_df]
            dict_number += 1
    elif p_comp == 1:
        chunked_datasets[0] = [peptide_df['injection_bins'].unique(), peptide_df]

    return(chunked_datasets)

if __name__ == "__main__":
    parallel_split()
