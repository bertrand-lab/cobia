# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 13:24:41 2018

@author: Scott
"""

import pandas as pd
import math as math
import string as str


#%%

peptide_master = pd.read_csv('Run4_C4_2000ng_oligo_lc-retention-times.csv')

#%%

pep_list = peptide_master['peptide_sequence'].tolist()

#%%

pep_list_sub = pep_list[1:100]

#%%

pep_list = peptide_master['peptide_sequence'].str.split('-').tolist()

pep_list_real = []
for element_i in pep_list:
    if isinstance(element_i, list):
        pep_list_real.append(element_i[1])
        continue
    elif math.isnan(element_i):
        pep_list_real.append('')
        continue
    else:
        print(element_i)

#%%
pep_list_real_sub = pep_list_real[1:100]

#%%
charge_state = []
three_indicators = ['H', 'RP', 'KP']
for element_i in pep_list_real_sub:
    if isinstance(element_i, basestring):
        if any(x in three_indicators for x in element_i):
            charge_state.append(3)
            print(3)
        else:
            charge_state.append(2)
            print(2)
    elif math.isnan(element_i):
        charge_state.append('')
    else:
        print(element_i)
#%%
peptide_master['mz'] = peptide_master['mass']/charge_state