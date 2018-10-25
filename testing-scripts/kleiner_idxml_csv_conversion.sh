#!/bin/bash
# converting idxml to csv

DIR='../data/kleiner_data/'
for FILE in "$DIR"*_IDF.idXML
do
        echo "Processing $FILE file..."
	temp_string=${FILE/.idXML}
	TextExporter -in $FILE -out $temp_string'.csv' -id:peptides_only -id:first_dim_rt -separator ','
done
