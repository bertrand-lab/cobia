"""
This script is for:

Inputting an idXML file, and subsetting the peptides, and producing a txt file. This is so that when an
SVM model is being trained for predicting retention times you don't need to use so many observations, which
would then slow down model training.

"""

from pyopenms import *
import pandas as pd
import csv
from argparse import ArgumentParser

#inputting command line arguments

#parser = ArgumentParser()
#parser.add_argument("-f", "--file", dest = "filename", help = "Input idXML file for subsetting data")
#parser.add_argument("-n", "--name", dest = "output_name", help = "Character string of .txt file output")
#parser.add_argument("-s", "--subset", dest = "subsetnth", help = "The degree of subsetting, every s-th peptide from the idXML is selected", type = int)

#args = parser.parse_args()

def subsample_idxml(filename, output_name, subset):

# Goal of script:
## Read in idXML files
## Subsample idXML files every 1/1000 peptide
## write subsampled idXML file

    # Reading in file
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(filename, protein_ids, peptide_ids)

    # rt vector and seq vector
    pep_rts = []
    pep_seqs = []

    for peptide_id in peptide_ids:
        pep_rts.append(peptide_id.getRT())
        for hit in peptide_id.getHits():
            pep_seq_i = hit.getSequence()
            pep_seqs.append(pep_seq_i)

    pep_rts_sub = pep_rts[::subsetnth]
    pep_seqs_sub = pep_seqs[::subsetnth]

    # combining subsampled pep seq and rt lists
    zseq = zip(pep_seqs_sub, pep_rts_sub)

    with open(output_name, 'w') as f:
        writer = csv.writer(f, delimiter = '\t')
        writer.writerows(zseq)

    quit()

if __name__ == "__main__":
    subsample_idxml()

