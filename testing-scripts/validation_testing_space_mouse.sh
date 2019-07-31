# modify peptides and get peptide physicochemical characteristics

DIR='../data/space_mouse_data/'

cobia openms_modelled_rt -f "$DIR"UP000000589_10090_rt_oligo.csv -n "$DIR"UP000000589_10090_rt_oligo

cobia cofrag_prediction -l "$DIR"dda_params_space_mouse.csv -f "$DIR"UP000000589_10090_rt_oligo_lc-retention-times.csv -n "$DIR"UP000000589_10090_rt_oligo_cofrag --global global

