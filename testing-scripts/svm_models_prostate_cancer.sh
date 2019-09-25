DIR='../data/prostate_cancer_data/'

# take protein coding genome (MS search database) and digest it

cobia database_trypsin -f "$DIR"Human_uniprot_09082016_all_single.fasta -n "$DIR"Human_uniprot_09082016_all_trypsin_single.txt -c no-write

# fit retention time predictor model with RTModel

#RTModel -in "$DIR"prostate_cancer_peptides_formatted.txt -out "$DIR"prostate_cancer_model_oligo.txt -kernel_type 'OLIGO' -total_gradient_time 7200

# predict retention times with model

RTPredict -in_text "$DIR"Human_uniprot_09082016_all_trypsin_single.txt -out_text:file "$DIR"Human_uniprot_09082016_rt_oligo.csv -svm_model "$DIR"prostate_cancer_model_oligo.txt -total_gradient_time 7200

