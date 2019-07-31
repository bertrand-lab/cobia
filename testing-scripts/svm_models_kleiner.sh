#!/bin/bash
# testing RT Model SVM kernal types

DIR='../data/kleiner_data/'

#conda activate pyteo_27

# tryptic digest of database file

cobia database_trypsin -f "$DIR"Mock_Comm_RefDB_V3.fasta -n "$DIR"Mock_Comm_RefDB_V3_trypsin.txt -c no-write

for FILE in "$DIR"*ng_FDR.idXML
do
        echo "Processing $FILE file..."
        temp_string=${FILE/.idXML/}
        echo $temp_string

	# setting retention time total for model fitting and predictions

	retention_time=15600 #260 mins

	if [ "${FILE}" == "../data/kleiner_data/Run1_C4_832ng_FDR.idXML" ]
	then
		retention_time=27600 #460 mins
	fi

# subset RTs for fitting a RT model, there are too many straight from the idXML files
       # python open_ms_rt.py -f $FILE  -n $temp_string'subset_peps.txt' -s 40
	cobia subsample_idxml -f $FILE  -n $temp_string'subset_peps.txt' -s 40

# fit the different RTModels with different kernels

	RTModel -in $temp_string'subset_peps.txt' -out $temp_string'model_oligo.txt' -kernel_type 'OLIGO' -total_gradient_time $retention_time &
	RTModel -in $temp_string'subset_peps.txt' -out $temp_string'model_linear.txt' -kernel_type 'LINEAR' -total_gradient_time $retention_time &
	RTModel -in $temp_string'subset_peps.txt' -out $temp_string'model_rbf.txt' -kernel_type 'RBF' -total_gradient_time $retention_time &
	RTModel -in $temp_string'subset_peps.txt' -out $temp_string'model_poly.txt' -kernel_type 'POLY' -total_gradient_time $retention_time &

	wait

# predict the RTs for the database for each SVM model

	RTPredict -in_text ../data/kleiner_data/Mock_Comm_RefDB_V3_trypsin.txt -out_text:file $temp_string'_tryptic_peptide_rt_oligo.csv' -svm_model $temp_string'model_oligo.txt' -total_gradient_time $retention_time &
        RTPredict -in_text ../data/kleiner_data/Mock_Comm_RefDB_V3_trypsin.txt -out_text:file $temp_string'_tryptic_peptide_rt_linear.csv' -svm_model $temp_string'model_linear.txt' -total_gradient_time $retention_time &
        RTPredict -in_text ../data/kleiner_data/Mock_Comm_RefDB_V3_trypsin.txt -out_text:file $temp_string'_tryptic_peptide_rt_rbf.csv' -svm_model $temp_string'model_rbf.txt' -total_gradient_time $retention_time &
        RTPredict -in_text ../data/kleiner_data/Mock_Comm_RefDB_V3_trypsin.txt -out_text:file $temp_string'_tryptic_peptide_rt_poly.csv' -svm_model $temp_string'model_poly.txt' -total_gradient_time $retention_time &

done

