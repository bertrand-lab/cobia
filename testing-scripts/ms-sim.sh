#!/bin/bash

foo () {
	local file_name='../data/kleiner_data/Mock_Comm_RefDB_V3_abundance'$1
	echo $file_name
# assigning random values to proteins drawn from a weibull distribution
	python FASTAProteinAbundanceSampling.py ../data/kleiner_data/Mock_Comm_RefDB_V3.fasta $file_name'.fasta' --family weibull --mu 0.01 --sigma 5
# digestion, LC, MSMS production of mzML file
	MSSimulator -in $file_name'.fasta' -out $file_name'.mzML'
}

for (( c=1; c<11; c++ )); do foo "$c" & done

