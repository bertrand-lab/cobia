#!/bash/bin

# read in fasta database and write out a trypsin digested database, both annotated with contigs and without contigs. 
# this step is not required for cofragmentation prediction, but is used to figure out tryptic peptides that would be useful
# for examining taxon-specific proteins.

cobia database_trypsin -f ../data/bertrand_data/orfs.filtered.pep.fasta -n ../data/bertrand_data/orfs.filtered.pep.trypsin_n.txt -c no-write

cobia database_trypsin -f ../data/bertrand_data/orfs.filtered.pep.fasta -n ../data/bertrand_data/orfs.filtered.pep.trypsin_wcontigs_n.txt -c write

# retention time prediction with BioLCCC. The same gradient and LC characteristics are used from Broberg et al 2018.

#cobia peptide_mod_biolccc_rt_prediction -f ../data/bertrand_data/orfs.filtered.pep.fasta -l ../data/broberg_data/lc_params_broberg2.csv -n ../data/bertrand_data/orfs.filtered.pep.trypsin -g ../data/broberg_data/custom_gradient_broberg2.csv

# examine targeted cofragmentation scores

#cobia cofrag_prediction -f ../data/bertrand_data/orfs.filtered.pep.trypsin_lc-retention-times.csv -l ../data/broberg_data/dda_params_broberg.csv --global targeted -n ../data/bertrand_data/orfs.filtered.pep.trypsin_targeted_frag_cyl_metE -t ../data/bertrand_data/frag_cyl_peps_metE.csv

#cobia cofrag_prediction -f ../data/bertrand_data/orfs.filtered.pep.trypsin_lc-retention-times.csv -l ../data/broberg_data/dda_params_broberg.csv --global global -n ../data/bertrand_data/orfs.filtered.pep.trypsin_global


