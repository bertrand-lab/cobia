#!/bash/bin

# read in

#python ../cobia/cobia.py database_trypsin -f ../data/bertrand_tfg_data/assembly.orf.fasta -n ../data/bertrand_tfg_data/assembly.orf.trypsin_digest.txt

#cobia database_trypsin -f ../data/bertrand_data/orfs.filtered.pep.fa -n ../data/bertrand_data/orfs.filtered.pep.trypsin_digest.txt

# predict retention times of database with BioLCCC

#cobia peptide_mod_biolccc_rt_prediction -f ../data/bertrand_data/orfs.filtered.pep.fasta -l ../data/broberg_data/lc_params_broberg2.csv -n ../data/bertrand_data/assembly.orf -g ../data/broberg_data/custom_gradient_broberg2.csv

# determine cofragmentation of a few targeted peptides we are interested in:

cobia cofrag_prediction -f ../data/bertrand_data/assembly.orf_lc-retention-times.csv -l ../data/broberg_data/dda_params_broberg.csv --global global -n ../data/bertrand_data/assembly.orf.fasta_cofrag

# examine the cofragmentation scores and manually look into other cofragmenting peptides


