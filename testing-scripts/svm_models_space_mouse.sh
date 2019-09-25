DIR='../data/space_mouse_data/'

# take protein coding genome (MS search database) and digest it

cobia database_trypsin -f "$DIR"UP000000589_10090_single.fasta -n "$DIR"UP000000589_10090_trypsin_single.txt -c no-write

# fit retention time predictor model with RTModel

#RTModel -in "$DIR"space_mouse_peptides_formatted.txt -out "$DIR"space_mouse_model_oligo.txt -kernel_type 'OLIGO' -total_gradient_time 9300

# predict retention times with model

RTPredict -in_text "$DIR"UP000000589_10090_trypsin_single.txt -out_text:file "$DIR"UP000000589_10090_rt_oligo.csv -svm_model "$DIR"space_mouse_model_oligo.txt -total_gradient_time 9300

