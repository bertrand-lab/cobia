
#!/bin/bash
# using openms modelled RTs for cofragmentation risk prediction

cobia openms_modelled_rt -f ../data/aylward_data/aylward_tryptic_peptide_rt_linear.csv -n ../data/aylward_data/aylward_tryptic_peptide_linear &
cobia openms_modelled_rt -f ../data/aylward_data/aylward_tryptic_peptide_rt_oligo.csv -n ../data/aylward_data/aylward_tryptic_peptide_oligo &
cobia openms_modelled_rt -f ../data/aylward_data/aylward_tryptic_peptide_rt_rbf.csv -n ../data/aylward_data/aylward_tryptic_peptide_rbf

cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/aylward_tryptic_peptide_linear_lc-retention-times.csv -n ../data/aylward_data/aylward_linear_cofrag --global global
cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/aylward_tryptic_peptide_oligo_lc-retention-times.csv -n ../data/aylward_data/aylward_oligo_cofrag --global global
cobia cofrag_prediction -l ../data/aylward_data/dda_params_aylward.csv -f ../data/aylward_data/aylward_tryptic_peptide_rbf_lc-retention-times.csv -n ../data/aylward_data/aylward_rbf_cofrag --global global

