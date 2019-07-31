
# modify peptides and get peptide physicochemical characteristics

DIR='../data/prostate_cancer_data/'

cobia openms_modelled_rt -f "$DIR"Human_uniprot_09082016_rt_oligo.csv -n "$DIR"Human_uniprot_09082016_rt_oligo

cobia cofrag_prediction -l "$DIR"dda_params_prostate_cancer.csv -f "$DIR"Human_uniprot_09082016_rt_oligo_lc-retention-times.csv -n "$DIR"Human_uniprot_09082016_rt_oligo_cofrag --global global
