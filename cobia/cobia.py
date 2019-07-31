# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 10:39:17 2017

This script directs the command line tool, cobia, to other sub-programs.

@author: Scott McCainS
"""
import pandas as pd
import numpy as np
import numbers
import multiprocessing as mp
import warnings
import os
import time
import csv
from sys import argv
import math

import pyteomics
from pyteomics import electrochem
#from pyteomics import biolccc
from pyteomics import mass
import re
from argparse import ArgumentParser
from pyopenms import *
import argparse

from .cobia_cofrag import cobia_cofrag
from .cobia_parallel import parallel_split
from .cobia_functions import cofrag_ion_counter_sparse_para
from .cobia_functions import cofrag_ion_counter_targeted
from ._peptide_mod_biolccc_rt_prediction import peptide_mod_biolccc_rt_prediction
from ._openms_modelled_rt import openms_modelled_rt
from ._database_trypsin import database_trypsin
from ._open_ms_rt import subsample_idxml

def main():

	parser = argparse.ArgumentParser(prog = 'main')

	# Subparser for cofragmentation prediction
	subparsers = parser.add_subparsers(title = 'subcommands', description = "The following subcommands are possible: cofrag_prediction, peptide_mod_biolccc_rt_prediction, database_trypsin, openms_modelled_rt",
						dest = 'subparser_name')

	cobia_cofrag_parser = subparsers.add_parser("cofrag_prediction")
	cobia_cofrag_parser.add_argument("-f", "--file", dest = "lcfilename", help = "Input is a .csv file of peptides with predicted retention times ('rts'), peptide mass ('mass').")
	cobia_cofrag_parser.add_argument("-l", "--ddafile",
                    dest = "ddaparams",
                    help = "Input of csv file that has the MS sampling parameters: max_injection_time, precursor_selection_window, ion_peak_width, number_of_parallel, and every_nth")
	# add argument to do global or targeted
	cobia_cofrag_parser.add_argument("--global", dest = "globalcofrag", help = "A string, either 'global' or 'targeted' (no quotations), denoting which approach to use.")
	cobia_cofrag_parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of output file name")
	cobia_cofrag_parser.add_argument("-t", "--targets", dest = "target_df", help = "File with target peptides. .csv required, with only one column, with the header pep_seq")

	# Subparser for peptide modification and biolccc retention time prediction
	peptide_mod_biolccc_rt_prediction_parser = subparsers.add_parser("peptide_mod_biolccc_rt_prediction")
	peptide_mod_biolccc_rt_prediction_parser.add_argument("-f", "--file", dest = "fastafilename", help = "Input is a .fasta file of protein coding sequences (Amino acids accepted only).")
	peptide_mod_biolccc_rt_prediction_parser.add_argument("-l", "--lcfile",
	                    dest = "lcparams",
                            help = "Input of csv file that has the lc parameters: column_length, column_diameter, column_pore_size, second_solvent_concentration_a, second_solvent_concentration_b, gradient_0, gradient_1, gradient_2, flow_rate, code_format, linear, model")
	peptide_mod_biolccc_rt_prediction_parser.add_argument('-g', '--custom_gradient', help = 'Input of csv file which defines the custom gradient (non linear LC gradient)', dest = "customgradient")
	peptide_mod_biolccc_rt_prediction_parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of output file name")

	# Subparser for tryptic peptides from fasta file
	database_trypsin_parser = subparsers.add_parser("database_trypsin")
	database_trypsin_parser.add_argument("-f", "--file", dest = "filename", help = "Input fasta file of database (AA sequences)")
	database_trypsin_parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of .txt file output with all tryptic peptides from database")
        database_trypsin_parser.add_argument("-c", "--add_contigs", dest = "add_contigs", help = '"write" or "no-write" regarding contig IDs to be added to the output file. For RTPredict, choose FALSE. This setting is required!')

	# Subparser for peptide modification script after RTModel from OpenMS
	openms_modelled_rt_parser = subparsers.add_parser("openms_modelled_rt")
	openms_modelled_rt_parser.add_argument("-f", "--file", dest = "rtfilename", help = "Input of csv file that contains tryptic peptides and retention times from database_tryptic.py")
	openms_modelled_rt_parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of output file name")

	subsample_idxml_parser = subparsers.add_parser("subsample_idxml")
	subsample_idxml_parser.add_argument("-f", "--file", dest = "filename", help = "Input idXML file for subsetting data")
	subsample_idxml_parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of .txt file output")
	subsample_idxml_parser.add_argument("-s", "--subset", dest = "subsetnth", help = "The degree of subsetting, every s-th peptide from the idXML is selected", type = int)

	args = parser.parse_args()

	if args.subparser_name == "cofrag_prediction":
		cobia_cofrag(lcfilename = args.lcfilename, ddaparams = args.ddaparams, globalcofrag = args.globalcofrag, output_name = args.output_name, target_df = args.target_df)
	elif args.subparser_name == "peptide_mod_biolccc_rt_prediction":
		peptide_mod_biolccc_rt_prediction(lc_params_file = args.lcparams, fasta_file_name = args.fastafilename, custom_gradient = args.customgradient, output_name = args.output_name)
	elif args.subparser_name == "database_trypsin":
		database_trypsin(filename = args.filename, output_name = args.output_name, contig_ids = args.add_contigs)
	elif args.subparser_name == "openms_modelled_rt":
		openms_modelled_rt(rtfilename = args.rtfilename, output_name = args.output_name)
	elif args.subparser_name == "subsample_idxml":
		subsample_idxml(filename = args.filename, output_name = args.output_name, subset = args.subsetnth)

if __name__ == "__main__":
	main()
