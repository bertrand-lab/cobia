# -*- coding: utf-8 -*-
"""
This script is for:

Retention time prediction using BioLCCC model, getting mass estimates, and modifying peptides.
The input is a fasta file, and output is a .csv dataframe with the AA sequence, charge at pH = 7, isoelectric point,
mass, and RT.

"""

import pyteomics
from pyteomics import electrochem
from pyteomics import biolccc
from pyteomics import mass
import re
from argparse import ArgumentParser

#Dictionary of the Genetic Code#
### This function (Codon_to_Aminoacid) is written by Noor Youssef. Perhaps later on we would want to use codon fasta files with
### nucleotides as an input. For now, it's just fasta files with AA inputs.

Codon_AA =  {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
             "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
             "TAT":"Y", "TAC":"Y", "TAA" :"", "TAG" : "",
             "TGT":"C", "TGC":"C", "TGG":"W", "TGA":"",
             "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
             "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
             "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
             "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
             "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
             "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
             "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
             "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
             "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
             "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
             "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
             "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

def Codon_to_Aminoacid(CodonSeq):
    '''
        Given a codon sequence, convert it to an amino acid sequence.
    '''
    AASeq = "".join([Codon_AA[CodonSeq[i:i+3]] for i in range(0,len(CodonSeq),3)])
    return AASeq


# read in lc paramter file
def peptide_mod_biolccc_rt_prediction(lc_params_file, fasta_file_name, custom_gradient, output_name):
	lc_params = pd.read_csv(lc_params_file)

	all_required_params = ['column_length',
	                       'column_diameter',
	                       'column_pore_size',
	                       'second_solvent_concentration_a',
	                       'second_solvent_concentration_b',
	                       'gradient_0',
	                       'gradient_1',
	                       'gradient_2',
	                       'flow_rate',
	                       'code_format',
	                       'linear',
			       'model']

	# Check if all parameters are in parameter input file.
	if sorted(all_required_params) != sorted(list(lc_params.keys())):
	    raise NameError('Error in parameter LC input file, check for typos or missing parameter.')

	# TRUE OR FALSE statement about whether the fasta file is in codons or in amino acids
	# Currently there is not a method in place, to use a nucleotide sequence fasta file as an input.
	code_format = lc_params['code_format'][0]
	linear_gradient = lc_params['linear'][0]

	# if not a linear gradient, a gradient file must be supplied.
	if not linear_gradient:
	    gradient_file = pd.read_csv(custom_gradient)

	# which type of model to use for prediction (from TFA or FA)
	model_type = lc_params['model'][0]

	if model_type == 'FA':
	    print('formic acid')
	elif model_type == 'TFA':
	    print('tri')

	# Initialize empty dictionary of contig names and sequences:
	seq_df = pd.DataFrame(columns = ['contigs', 'seq'])

	# Initialize empty lists of sequences and contigs:
	seq_vec = []
	contig_vec = []

	# Initalize variable that will contain the name of each sequence:
	last_seq = None

	# Reading in fasta file
	fasta_in = open(fasta_file_name, 'r')

	for line in fasta_in:
	    # Strip the line:
	    line = line.strip()
	    # If the line is blank, move on.
	    if len(line) == 0: # blank line
	        continue
	    # If the line is a header, record the header as last_seq
	    elif line[0] == ">": # header-line
	        last_seq = line[1:]
	    # If the line is a sequence, record the sequence:
	    else: # sequence line
	        # separate if statements for if the fasta file was input as amino acids or as genes or as mrna. Note that code_format == 'genes' and code_format == 'rna' are not functional yet. 
	        if(code_format == 'genes'):
	            aa_line = Codon_to_Aminoacid(line)
	            cleaved_line = pyteomics.parser.cleave(str(aa_line), pyteomics.parser.expasy_rules['trypsin'])
	            cleaved_line = list(cleaved_line)
	        elif(code_format == 'rna'):
	            removed_u = line.relace('U', 'T')
	            aa_line = Codon_to_Aminoacid(removed_u)
	            cleaved_line = pyteomics.parser.cleave(str(aa_line), pyteomics.parser.expasy_rules['trypsin'])
	            cleaved_line = list(cleaved_line)
	        elif(code_format == 'aas'):
	            # Digest with trypsin:
	            cleaved_line = pyteomics.parser.cleave(str(line), pyteomics.parser.expasy_rules['trypsin'])
	            cleaved_line = list(cleaved_line)
	        # If the peptide is shorter than 5 amino acids long, then we remove it fromt the dataset:
	        for tryp_pep in cleaved_line:
	            if len(tryp_pep) < 5:
	                continue
	            seq_vec.append(tryp_pep)
	            contig_vec.append(last_seq)

	# Close the fasta file:
	fasta_in.close()

	print('Removing xs and *s from seqs...')

	contig_vec_pd = pd.Series(contig_vec, name = 'contig')
	# Adding in the modification terms for the termini:
	seq_vec_terms = [central_pep + '-OH' for central_pep in seq_vec]
	# Removing contigs with unknown amino acid (X) or selenocysteine (U):
	stars_removed_peps = []
	for starred_peptide in seq_vec_terms:
            line_new = starred_peptide
            if '*' in line_new:
                continue
	    #some peptides have unknown amino acids, remove them.
	    if 'X' in line_new:
	        continue
	    if 'U' in line_new:
	        continue
	    stars_removed_peps.append(line_new)

	# Changing B to asparagine
	b_removed_peps = []
	for b_peptide in stars_removed_peps:
	    line_new = re.sub('B', 'N', b_peptide)
	    b_removed_peps.append(line_new)

	# Changing Z to glutamine
	z_removed_peps = []
	for z_peptide in b_removed_peps:
	    line_new = re.sub('Z', 'Q', z_peptide)
	    z_removed_peps.append(line_new)

	# Removing contigs that have an unknown amino acid (X), or selenocysteine ('U')
	contig_vec_no_x = []
	for contig_name in range(len(contig_vec)):
	    if 'X' in seq_vec_terms[contig_name]:
	        continue
	    if 'U' in seq_vec_terms[contig_name]:
	        continue
            if '*' in seq_vec_terms[contig_name]:
                continue

	    temp_contig = contig_vec[contig_name]
	    contig_vec_no_x.append(temp_contig)

	# Modifying peptides: oxidation of methionine, carbamidomethylation of cysteine, acetylation of N terminal (this one was done upstream)

	print('Modifying peptides...')
	mod_pep = []
	for tryp_pep in z_removed_peps:
	    test_iso = pyteomics.parser.isoforms(tryp_pep,
	                                         fixed_mods = {'ox':['M'], 'cam':['C']},
	                                         show_unmodified_termini = True)
	    for blah in test_iso:
	        mod_pep.append(blah)

	# Modified amino acid dictionary for mass calculation:

	aa_comp = dict(mass.std_aa_comp)
	aa_comp['Ac-'] = mass.Composition({'C': 2, 'H': 3, 'N': 0, 'O': 1, 'P': 0})
	aa_comp['cam'] = mass.Composition({'C': 2, 'H': 3, 'N': 1, 'O': 1, 'P': 0})
	aa_comp['ox'] = mass.Composition({'O':1})

	# Calculate peptide isoelectric points, masses, and charge at pH = 7. Note that we do not use the isoelectric point or charge from this point on, but used it for examining other predictive components of apparent cofragmentation bias.

	print('Calculating peptide physicochemical properties...')
	iso_electric_points = []
	pep_charges = []
	pep_mass = []
	i = 0

	for peptide in mod_pep:
	    peptide_isoelectric_point = electrochem.pI(peptide)
	    peptide_charge = electrochem.charge(peptide, 7)
	    peptide_mass = mass.calculate_mass(sequence = peptide, aa_comp = aa_comp)
	    pep_charges.append(peptide_charge)
	    iso_electric_points.append(peptide_isoelectric_point)
	    pep_mass.append(peptide_mass)
	    i += 1

	print('LC-retention time prediction with the following parameters:')

	print(lc_params)

	# Column length:
	column_length = lc_params['column_length'][0]
	if isinstance(column_length, numbers.Number) != True:
	    raise NameError('Error in parameter input file, column_length takes only Numeric.')

	# Column diameter:
	column_diameter = lc_params['column_diameter'][0]
	if isinstance(column_diameter, numbers.Number) != True:
	    raise NameError('Error in parameter input file, column_diameter takes only Numeric.')

	# Column pore size
	column_pore_size = lc_params['column_pore_size'][0] # 0.11 minutes
	if isinstance(column_pore_size, numbers.Number) != True:
	    raise NameError('Error in parameter input file, column_pore_size takes only Numeric.')

	second_solvent_concentration_a = lc_params['second_solvent_concentration_a'][0]
	if isinstance(second_solvent_concentration_a, numbers.Number) != True:
	    raise NameError('Error in parameter input file, second_solvent_concentration_a takes only Numeric.')

	second_solvent_concentration_b = lc_params['second_solvent_concentration_b'][0]
	if isinstance(second_solvent_concentration_b, numbers.Number) != True:
	    raise NameError('Error in parameter input file, second_solvent_concentration_b takes only Numeric.')

	gradient_0 = lc_params['gradient_0'][0]
	if isinstance(gradient_0, numbers.Number) != True:
	    raise NameError('Error in parameter input file, gradient_0 takes only Numeric.')

	gradient_1 = lc_params['gradient_1'][0]
	if isinstance(gradient_1, numbers.Number) != True:
	    raise NameError('Error in parameter input file, gradient_1 takes only Numeric.')

	gradient_2 = lc_params['gradient_2'][0]
	if isinstance(gradient_2, numbers.Number) != True:
	    raise NameError('Error in parameter input file, gradient_2 takes only Numeric.')

	flow_rate = lc_params['flow_rate'][0]
	if isinstance(flow_rate, numbers.Number) != True:
	    raise NameError('Error in parameter input file, flow_rate takes only Numeric')

	# biolccc predicting RT times
	myChromoConditions = biolccc.ChromoConditions()

	# The column length in mm.
	myChromoConditions.setColumnLength(column_length)

	# The internal column diameter in mm.
	myChromoConditions.setColumnDiameter(column_diameter)

	# The average pore size in A.
	myChromoConditions.setColumnPoreSize(column_pore_size)

	# The concentration of the eluting solvent (ACN for the reversed
	# phase) in component A in %.
	myChromoConditions.setSecondSolventConcentrationA(second_solvent_concentration_a)

	# The concentration of the eluting solvent (ACN for the reversed
	# phase) in component B in %.
	myChromoConditions.setSecondSolventConcentrationB(second_solvent_concentration_b)

	# The shape of the gradient. The example is a linear gradient
	# from gradient_0% to gradient_1% of component B over gradient_2 minutes.

	if linear_gradient:
	    myChromoConditions.setGradient(biolccc.Gradient(gradient_0, gradient_1, gradient_2))
	else:
	    # loop that goes through and sets a custom gradient. another gradient file is required as the argv[4] file.
	    myGradient = biolccc.Gradient()
	    # An older version of this was more static, and left in the comments below to demonstrate what this loop is doing:
	    for set_point in range(len(gradient_file.columns)):
	        myGradient.addPoint(gradient_file.iloc[0, set_point], 
	                            gradient_file.iloc[1, set_point])
	    myChromoConditions.setGradient(myGradient)

	    # The following gradient is an exponential function increasing from gradient_0
	    # to 100, specifically for the Aylward testing datasetself.
	    # def exp_function(x):
	    #   x1 = math.pow(x, 2)//100
	    # this is the function used to compute these setpoints.
	    #myGradient = biolccc.Gradient()
	    #myGradient.addPoint(0.0, gradient_0)
	    #myGradient.addPoint(15.0, 2.0)
	    #myGradient.addPoint(30.0, 9.0)
	    #myGradient.addPoint(45.0, 20.0)
	    #myGradient.addPoint(60.0, 36.0)
	    #myGradient.addPoint(75.0, 56.0)
	    #myGradient.addPoint(90.0, 81.0)
	    #myGradient.addPoint(gradient_2, gradient_1)
	    #myChromoConditions.setGradient(myGradient)

	# The flow rate in ml/min.
	myChromoConditions.setFlowRate(flow_rate)

	print('Calculating retention times...')

	# Designating BioLCCC model to use:
	if model_type == 'TFA':
		model_to_use = biolccc.rpAcnTfaChain
	elif model_type == 'FA':
		model_to_use = biolccc.rpAcnFaRod

	peptide_rts = []
	i = 0

	print('Calculating retention times...')
	for tryp_pep in mod_pep:
	    rt_temp = biolccc.calculateRT(tryp_pep,
	                                  model_to_use,
	                                  myChromoConditions)

	    peptide_rts.append(rt_temp)
	    i += 1

	# Combining the sequences, times, and physicochemical characteristics.
	peptides_pd = pd.Series(z_removed_peps, name = 'peptide_sequence')
	peptide_rts = pd.Series(peptide_rts, name = 'rts')
	iso_electric_points_pd = pd.Series(iso_electric_points, name = 'iso_point')
	pep_charges_pd = pd.Series(pep_charges, name = 'charge')
	pep_mass_pd = pd.Series(pep_mass, name = 'mass')
	contig_pd = pd.Series(contig_vec_no_x, name = 'contig')

	peptide_dataframe = pd.concat([peptides_pd,
	                               peptide_rts,
	                               iso_electric_points_pd,
	                               pep_charges_pd,
	                               pep_mass_pd,
	                               contig_pd], axis = 1)

	current_date = time.strftime("%Y-%m-%d")

	custom_name = output_name
	file_name = custom_name + '_lc-retention-times.csv'
	peptide_dataframe.to_csv(file_name)
