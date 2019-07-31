#!/bin/bash
# Database searching for Aylward data using MSGF+

DIR='../data/aylward_data/'
for FILE in ls "$DIR"*mzML
do
	echo "Processing $FILE file..."
	java -Xmx16G -jar /software/MSGFPlus/MSGFPlus.jar -s $FILE -d ../data/aylward_data/combined_metagenomic_aylward.faa -tda 1 -thread 32
done

#DIR='../data/broberg_data/'

