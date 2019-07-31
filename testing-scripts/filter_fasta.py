#!/usr/bin/env python3

# python script for removing sequences with a given header

from Bio import SeqIO
import sys

fasta_file = sys.argv[1]
remove_string = sys.argv[2]
result_file = sys.argv[3]

with open(result_file, 'w') as f:
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if remove_string not in seq.id:
            SeqIO.write(seq, f, 'fasta')


