
# script to go in an download all the mass spec data from Schmidt et al 2015, the quantitative e coli proteomics paper.
# ms_data_file_names_ecoli.txt was from copying the FTP from the PRIDE repo, then scraping just the names using cat ms_data_file_names.txt | awk '{ print $2 }' > ms_data_file_names_ecoli.txt

for word in $(cat ../data/ecoli_data/ms_data_file_names_ecoli.txt); 
    do echo $word;
    tempname='prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2015/10/PXD000498/'$word
    echo $tempname
    ../../../.aspera/connect/bin/ascp -TQ -l200m -P 33001 -i "../../../.aspera/connect/etc/asperaweb_id_dsa.openssh" $tempname ../data/ecoli_data/mass_spec_data/
done

