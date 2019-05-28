# -*- coding: utf-8 -*-
"""
This script is for:

Taking the input from OpenMS modelled retention times, and outputting a csv file with peptides, rts,
and masses which would then be used for cofragmentation prediction.

"""

import pyteomics
from pyteomics import electrochem
from pyteomics import biolccc
from pyteomics import mass
from pyopenms import *
import re
import pandas as pd
import time
from sys import argv
import numbers
import os

#from argparse import ArgumentParser

#inputting command line arguments

#parser = ArgumentParser()
#parser.add_argument("-f", "--file", dest = "rtfilename", help = "Input of csv file that contains tryptic peptides and retention times from database_tryptic.py")
#parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of output file name")
#args = parser.parse_args()

# read in txt file with RTs and Peptide sequences as seq_vec (list) and rt_vec (list)

def openms_modelled_rt(rtfilename, output_name):

	seq_rt_df = pd.read_csv(rtfilename, names = ['seq_rt'])
	df = pd.DataFrame(seq_rt_df.seq_rt.str.split(' ', 1).tolist(), columns = ['pep_seq', 'rts'])
	seq_vec = df['pep_seq'].tolist()
	peptide_rts = df['rts'].tolist()

	print('Removing xs and *s from seqs...')
	#contig_vec_pd = pd.Series(contig_vec, name = 'contig')
	seq_vec_terms = [central_pep + '-OH' for central_pep in seq_vec]
	# removing contigs with unknown amino acid (X) or selenocysteine (U)
	stars_removed_peps = []
	for starred_peptide in seq_vec_terms:
            line_new = starred_peptide
            # some peptides have unknown amino acids denoted as *, remove them. 
            if '*' in line_new:
                continue
	    #some peptides have unknown amino acids, remove them.
	    if 'X' in line_new:
	        continue
	    if 'U' in line_new:
	        continue
	    stars_removed_peps.append(line_new)

	#changing B to asparagine
	b_removed_peps = []
	for b_peptide in stars_removed_peps:
	    line_new = re.sub('B', 'N', b_peptide)
	    b_removed_peps.append(line_new)

	#changing Z to glutamine
	z_removed_peps = []
	for z_peptide in b_removed_peps:
	    line_new = re.sub('Z', 'Q', z_peptide)
	    z_removed_peps.append(line_new)

	# #modifying peptides: oxidation of methionine, carbamidomethylation of cysteine, acetylation of N terminal (this one was done upstream)

	print('Modifying peptides...')
	mod_pep = []
	for tryp_pep in z_removed_peps:
	    test_iso = pyteomics.parser.isoforms(tryp_pep,
	                                         fixed_mods = {'ox':['M'], 'cam':['C']},
	                                         show_unmodified_termini = True)
	    for blah in test_iso:
	        mod_pep.append(blah)

	#%%

	# # modified amino acid dictionary for mass calculation

	aa_comp = dict(mass.std_aa_comp)
	aa_comp['Ac-'] = mass.Composition({'C': 2, 'H': 3, 'N': 0, 'O': 1, 'P': 0})
	aa_comp['cam'] = mass.Composition({'C': 2, 'H': 3, 'N': 1, 'O': 1, 'P': 0})
	aa_comp['ox'] = mass.Composition({'O':1})

	#%%
	# calculate peptide isoelectric points, masses, and charge at pH = 7

	print('Calculating peptide physicochemical properties...')
	iso_electric_points = []
	pep_charges = []
	pep_mass = []
	i = 0

	# print(mod_pep)

	for peptide in mod_pep:
	    peptide_isoelectric_point = electrochem.pI(peptide)
	    peptide_charge = electrochem.charge(peptide, 7)
	    peptide_mass = mass.calculate_mass(sequence = peptide, aa_comp = aa_comp)
	    pep_charges.append(peptide_charge)
	    iso_electric_points.append(peptide_isoelectric_point)
	    pep_mass.append(peptide_mass)
	    # print(i)
	    i += 1

	# Combining the sequences, times, and physicochemical characteristics.
	peptides_pd = pd.Series(z_removed_peps, name = 'peptide_sequence')
	peptide_rts = pd.Series(peptide_rts, name = 'rts')
	iso_electric_points_pd = pd.Series(iso_electric_points, name = 'iso_point')
	pep_charges_pd = pd.Series(pep_charges, name = 'charge')
	pep_mass_pd = pd.Series(pep_mass, name = 'mass')

	peptide_dataframe = pd.concat([peptides_pd,
	                               peptide_rts,
	                               iso_electric_points_pd,
	                               pep_charges_pd,
	                               pep_mass_pd], axis = 1)

	current_date = time.strftime("%Y-%m-%d")

	custom_name = output_name
	file_name = custom_name + '_lc-retention-times.csv'
	print(file_name)
	peptide_dataframe.to_csv(file_name)
