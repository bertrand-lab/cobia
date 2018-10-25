
#!/bin/bash

# retention time prediction for broberg et al 2018 study metagenome, metatranscriptome, and metagenome with oak genome appended.

# retention time prediction with BIOLCCC
cobia peptide_mod_biolccc_rt_prediction -f ../data/broberg_data/mg_AT4.fasta -n ../data/broberg_data/broberg_mg_tryptic_peptide -l ../data/broberg_data/lc_params_broberg2.csv -g ../data/broberg_data/custom_gradient_broberg2.csv &
cobia peptide_mod_biolccc_rt_prediction -f ../data/broberg_data/mt_AT4.fasta -n ../data/broberg_data/broberg_mt_tryptic_peptide -l ../data/broberg_data/lc_params_broberg2.csv -g ../data/broberg_data/custom_gradient_broberg2.csv &
cobia peptide_mod_biolccc_rt_prediction -f ../data/broberg_data/combined_mg_qrob_oak_genome.fasta -n ../data/broberg_data/broberg_mg_qrob_oak_genome_tryptic_peptide -l ../data/broberg_data/lc_params_broberg2.csv -g ../data/broberg_data/custom_gradient_broberg2.csv

# cofragmentation prediction with cobia
cobia cofrag_prediction -l ../data/broberg_data/dda_params_broberg.csv -f ../data/broberg_data/broberg_mg_qrob_oak_genome_tryptic_peptide_lc-retention-times.csv -n ../data/broberg_data/broberg_mg_qrob_oak_genome_cofrag --global global
cobia cofrag_prediction -l ../data/broberg_data/dda_params_broberg.csv -f ../data/broberg_data/broberg_mg_tryptic_peptide_lc-retention-times.csv -n ../data/broberg_data/broberg_mg_cofrag --global global
cobia cofrag_prediction -l ../data/broberg_data/dda_params_broberg.csv -f ../data/broberg_data/broberg_mt_tryptic_peptide_lc-retention-times.csv -n ../data/broberg_data/broberg_mt_cofrag --global global


