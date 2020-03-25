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
from Bio import SeqIO

def database_trypsin(filename, output_name, contig_ids):

    seq_vec = []
    contig_vec = []

    # initalized variable that will contain the name of each sequence:
    last_seq = None

    #reading in fasta file
    for seq_record in SeqIO.parse(filename, "fasta"):
        cleaved_line = pyteomics.parser.cleave(str(seq_record.seq), pyteomics.parser.expasy_rules['trypsin'])
        cleaved_line = list(cleaved_line)
        for tryp_pep in cleaved_line:
            if len(tryp_pep) < 5:
                continue
            seq_vec.append(tryp_pep)
            contig_vec.append(seq_record.id)

    with open(output_name, 'w') as output:
        writer = csv.writer(output, lineterminator = '\n')
        i = 0
        for val in seq_vec:
            if contig_ids == 'write':
                writer.writerow([val, contig_vec[i]])
            elif contig_ids == 'no-write':
                writer.writerow([val])
            i += 1