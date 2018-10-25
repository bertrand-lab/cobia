#!/bash/bin

# read in

#cobia database_trypsin -f ../data/bertrand_data/orfs.filtered.pep.fasta -n ../data/bertrand_data/orfs.filtered.pep.trypsin.txt -c no-write

cobia database_trypsin -f ../data/bertrand_data/orfs.filtered.pep.fasta -n ../data/bertrand_data/orfs.filtered.pep.trypsin_wcontigs.txt -c write


# determine cofragmentation of a few targeted peptides we are interested in:

#cobia cofrag_prediction

#RTPredict -in_text ../data/bertrand_data/orfs.filtered.pep.trypsin.txt -out_text:file ../data/bertrand_data/orfs.filtered.pep.trypsin_rt_times.csv -svm_model ../data/kleiner_data/Run1_C4_832ng_FDRmodel_oligo.txt -total_gradient_time 15600

# examine the cofragmentation scores and manually look into other cofragmenting peptides

#cobia cofrag_prediction -f assembly.orf_lc-retention-times.csv -l ../broberg_data/dda_params_broberg.csv --global targeted -n assembly.orf_cofrag_global -t unique_candidate_peps_pseudo_kog_mnsod.csv

#cobia openms_modelled_rt -f ../data/bertrand_data/orfs.filtered.pep.trypsin_rt_times.csv -n ../data/bertrand_data/orfs.filtered.pep.trypsin_rt_times_mass.csv

#cobia cofrag_prediction -f ../data/bertrand_data/orfs.filtered.pep.trypsin_rt_times_mass.csv_lc-retention-times.csv -l ../data/kleiner_data/dda_params_kleiner_260.csv --global global -n ../data/bertrand_data/orfs.filtered.pep.trypsin_cofrag

