
#!/bin/bash
# using openms modelled RTs for cofragmentation risk prediction

# digesting the a cephalotes database
cobia database_trypsin -f ../data/aylward_data/2029527004.a.fasta -n ../data/aylward_data/2029527004.a_trypsin.txt -c no-write

retention_time=6000 #100 mins

# predicting retention times with three different kernels
RTPredict -in_text ../data/aylward_data/2029527004.a_trypsin.txt -out_text:file ../data/aylward_data/2029527004.a_rt_oligo.csv -svm_model ../data/aylward_data/aylward_model_oligo.txt -total_gradient_time $retention_time &
RTPredict -in_text ../data/aylward_data/2029527004.a_trypsin.txt -out_text:file ../data/aylward_data/2029527004.a_rt_linear.csv -svm_model ../data/aylward_data/aylward_model_linear.txt -total_gradient_time $retention_time &
RTPredict -in_text ../data/aylward_data/2029527004.a_trypsin.txt -out_text:file ../data/aylward_data/2029527004.a_rt_rbf.csv -svm_model ../data/aylward_data/aylward_model_rbf.txt -total_gradient_time $retention_time 

cobia openms_modelled_rt -f ../data/aylward_data/2029527004.a_rt_oligo.csv -n ../data/aylward_data/2029527004_tryptic_peptide_linear &
cobia openms_modelled_rt -f ../data/aylward_data/2029527004.a_rt_linear.csv -n ../data/aylward_data/2029527004_tryptic_peptide_oligo &
cobia openms_modelled_rt -f ../data/aylward_data/2029527004.a_rt_rbf.csv -n ../data/aylward_data/2029527004_tryptic_peptide_rbf

cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/2029527004_tryptic_peptide_linear_lc-retention-times.csv -n ../data/aylward_data/2029527004_linear_cofrag --global global
cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/2029527004_tryptic_peptide_oligo_lc-retention-times.csv -n ../data/aylward_data/2029527004_oligo_cofrag --global global
cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/2029527004_tryptic_peptide_rbf_lc-retention-times.csv -n ../data/aylward_data/2029527004_rbf_cofrag --global global

