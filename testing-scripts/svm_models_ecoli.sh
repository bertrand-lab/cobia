DIR='../data/ecoli_data/'

# modify the fasta file input sequences because they also contain reverse sequences

#python filter_fasta.py "$DIR"uniprot-ecoli.decoy.fasta 'REV_' "$DIR"uniprot-ecoli.no-decoy.fasta

# take protein coding genome (MS search database) and digest it

cobia database_trypsin -f "$DIR"uniprot-ecoli.no-decoy_single.fasta -n "$DIR"uniprot-ecoli.no-decoy_trypsin_single.txt -c no-write

# fit retention time predictor model with RTModel

# RTModel -in "$DIR"ecoli_peptides_formatted.txt -out "$DIR"ecoli_model_oligo.txt -kernel_type 'OLIGO' -total_gradient_time 10800

# predict retention times with model

RTPredict -in_text "$DIR"uniprot-ecoli.no-decoy_trypsin_single.txt -out_text:file "$DIR"uniprot-ecoli.no-decoy_rt_oligo.csv -svm_model "$DIR"ecoli_model_oligo.txt -total_gradient_time 10800

