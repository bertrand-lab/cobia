DIR='../data/fungal_pathogen_data/'

# take protein coding genome (MS search database) and digest it

cobia database_trypsin -f "$DIR"coccidioides_posadasii__rmscc_3488_1_proteins.fasta -n "$DIR"coccidioides_posadasii__rmscc_3488_1_proteins_trypsin.txt -c no-write

# fit retention time predictor model with RTModel

RTModel -in "$DIR"fungal_pathogen_peptides_formatted.txt -out "$DIR"fungal_pathogen_model_oligo.txt -kernel_type 'OLIGO' -total_gradient_time 6900

# predict retention times with model

RTPredict -in_text "$DIR"coccidioides_posadasii__rmscc_3488_1_proteins_trypsin.txt -out_text:file "$DIR"coccidioides_posadasii__rmscc_3488_1_proteins_trypsin_rt_oligo.csv -svm_model "$DIR"fungal_pathogen_model_oligo.txt -total_gradient_time 6900

