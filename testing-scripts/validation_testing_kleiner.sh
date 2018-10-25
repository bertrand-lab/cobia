#!/bin/bash
# running cofragmentation prediction and LC for kleiner data (260 LC and 460 LC runs of same mock community sample)

#cobia openms_modelled_rt -f ../data/kleiner_data/Run1_C4_832ng_FDR_tryptic_peptide_rt_rbf.csv -n ../data/kleiner_data/Run1_C4_832ng_rbf &
#cobia openms_modelled_rt -f ../data/kleiner_data/Run1_C4_832ng_FDR_tryptic_peptide_rt_oligo.csv -n ../data/kleiner_data/Run1_C4_832ng_oligo &
#cobia openms_modelled_rt -f ../data/kleiner_data/Run1_C4_832ng_FDR_tryptic_peptide_rt_linear.csv -n ../data/kleiner_data/Run1_C4_832ng_linear &

#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_rbf.csv -n ../data/kleiner_data/Run4_C4_2000ng_rbf &
#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_oligo.csv -n ../data/kleiner_data/Run4_C4_2000ng_oligo &
#cobia openms_modelled_rt -f ../data/kleiner_data/Run4_C4_2000ng_FDR_tryptic_peptide_rt_linear.csv -n ../data/kleiner_data/Run4_C4_2000ng_linear


cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_460.csv -f ../data/kleiner_data/Run1_C4_832ng_rbf_lc-retention-times.csv -n ../data/kleiner_data/Run1_C4_832ng_rbf_cofrag --global global
#cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_460.csv -f ../data/kleiner_data/Run1_C4_832ng_oligo_lc-retention-times.csv -n ../data/kleiner_data/Run1_C4_832ng_oligo_cofrag --global global
#cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_460.csv -f ../data/kleiner_data/Run1_C4_832ng_linear_lc-retention-times.csv -n ../data/kleiner_data/Run1_C4_832ng_linear_cofrag --global global

cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_260.csv -f ../data/kleiner_data/Run4_C4_2000ng_rbf_lc-retention-times.csv  -n ../data/kleiner_data/Run4_C4_2000ng_rbf_cofrag --global global
#cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_260.csv -f ../data/kleiner_data/Run4_C4_2000ng_oligo_lc-retention-times.csv -n ../data/kleiner_data/Run4_C4_2000ng_oligo_cofrag --global global
cobia cofrag_prediction -l ../data/kleiner_data/dda_params_kleiner_260.csv -f ../data/kleiner_data/Run4_C4_2000ng_linear_lc-retention-times.csv -n ../data/kleiner_data/Run4_C4_2000ng_linear_cofrag --global global
