#!/bin/bash
# testing RT Model SVM kernal types

DIR='../data/aylward_data/'

#conda activate pyteo_27

# tryptic digest of database file

retention_time=6000 #100 mins

cobia database_trypsin -f "$DIR"combined_metagenomic_aylward.fasta -n "$DIR"combined_metagenomic_aylward_trypsin.txt

# merging together all ID'd peptides from all idXML files.

IDMerger -in "$DIR"*[::digit::]_FDR.idXML -out "$DIR"aylward_merged_FDR.idXML

# no subsampling, so -s flag is every peptide ( -s 1 )
python open_ms_rt.py -f "$DIR"aylward_merged_FDR.idXML -n "$DIR"aylward_subset_peps.txt -s 1

RTModel -in "$DIR"aylward_subset_peps.txt -out "$DIR"aylward_model_oligo.txt -kernel_type 'OLIGO' -total_gradient_time $retention_time &
RTModel -in "$DIR"aylward_subset_peps.txt -out "$DIR"aylward_model_linear.txt -kernel_type 'LINEAR' -total_gradient_time $retention_time &
RTModel -in "$DIR"aylward_subset_peps.txt -out "$DIR"aylward_model_rbf.txt -kernel_type 'RBF' -total_gradient_time $retention_time &
RTModel -in "$DIR"aylward_subset_peps.txt -out "$DIR"aylward_model_poly.txt -kernel_type 'POLY' -total_gradient_time $retention_time &

wait

# predict the RTs for the database for each SVM model

RTPredict -in_text "$DIR"combined_metagenomic_aylward_trypsin.txt -out_text:file "$DIR"aylward_tryptic_peptide_rt_oligo.csv -svm_model "$DIR"aylward_model_oligo.txt -total_gradient_time $retention_time &
RTPredict -in_text "$DIR"combined_metagenomic_aylward_trypsin.txt -out_text:file "$DIR"aylward_tryptic_peptide_rt_linear.csv -svm_model "$DIR"aylward_model_linear.txt -total_gradient_time $retention_time &
RTPredict -in_text "$DIR"combined_metagenomic_aylward_trypsin.txt -out_text:file "$DIR"aylward_tryptic_peptide_rt_rbf.csv -svm_model "$DIR"aylward_model_rbf.txt -total_gradient_time $retention_time &
RTPredict -in_text "$DIR"combined_metagenomic_aylward_trypsin.txt -out_text:file "$DIR"aylward_tryptic_peptide_rt_poly.csv -svm_model "$DIR"aylward_model_poly.txt -total_gradient_time $retention_time &


