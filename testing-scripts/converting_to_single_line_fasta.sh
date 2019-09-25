# converting aylward fasta files
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../data/aylward_data/combined_metagenomic_aylward.fasta > ../data/aylward_data/combined_metagenomic_aylward_single.fasta

# converting kleiner fasta files
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../data/kleiner_data/Mock_Comm_RefDB_V3.fasta > ../data/kleiner_data/Mock_Comm_RefDB_V3_single.fasta

# converting prostate cancer files
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../data/prostate_cancer_data/Human_uniprot_09082016_all.fasta > ../data/prostate_cancer_data/Human_uniprot_09082016_all_single.fasta

# converting e coli file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../data/ecoli_data/uniprot-ecoli.no-decoy.fasta > ../data/ecoli_data/uniprot-ecoli.no-decoy_single.fasta

# converting space mouse file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ../data/space_mouse_data/UP000000589_10090.fasta > ../data/space_mouse_data/UP000000589_10090_single.fasta

