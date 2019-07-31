# modify peptides and get peptide physicochemical characteristics

DIR='../data/ecoli_data/'

cobia openms_modelled_rt -f "$DIR"uniprot-ecoli.no-decoy_rt_oligo.csv -n "$DIR"uniprot-ecoli.no-decoy_rt_oligo

cobia cofrag_prediction -l "$DIR"dda_params_ecoli.csv -f "$DIR"uniprot-ecoli.no-decoy_rt_oligo_lc-retention-times.csv -n "$DIR"uniprot-ecoli.no-decoy_rt_oligo_cofrag --global global
