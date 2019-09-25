# -*- coding: utf-8 -*-
"""
This script is for:

Inputting a fasta file that is the database of potential protein products. AA sequences are expected.

Output a txt file of tryptic peptides.

"""

import pyteomics
from pyteomics import electrochem
#from pyteomics import biolccc
from pyteomics import mass
from pyopenms import *
import re
import pandas as pd
import time
from sys import argv
import numbers
import os
import csv

def database_trypsin(filename, output_name, contig_ids):

	seq_vec = []
	contig_vec = []

	# initalized variable that will contain the name of each sequence:
	last_seq = None

	#reading in fasta file
	fasta_in = open(filename, 'r')

	#%%
	for line in fasta_in:
	    line = line.strip()
	    if len(line) == 0: # blank line
	        continue # skip to next iteration of loop
	    elif line[0] == ">": # header-line
	        last_seq = line[1:]
	    else: # sequence line
	        # separate if statements for if the fasta file was input as amino acids or as genes or as mrna
	        cleaved_line = pyteomics.parser.cleave(str(line), pyteomics.parser.expasy_rules['trypsin'])
	        cleaved_line = list(cleaved_line)
	        for tryp_pep in cleaved_line:
	            if len(tryp_pep) < 5:
	                continue
	            seq_vec.append(tryp_pep)
	            contig_vec.append(last_seq)
	fasta_in.close()

	with open(output_name, 'w') as output:
	    writer = csv.writer(output, lineterminator = '\n')
            i = 0
	    for val in seq_vec:
                if contig_ids == 'write':
	            writer.writerow([val, contig_vec[i]])
                elif contig_ids == 'no-write':
                    writer.writerow([val])
                i += 1
