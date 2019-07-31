for word in $(cat ../data/prostate_cancer_data/ms_data_file_names_prostate_cancer.txt);
    do echo $word;
    tempname='prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2018/01/PXD008407/'$word
    echo $tempname
    ../../../.aspera/connect/bin/ascp -TQ -l200m -P 33001 -i "../../../.aspera/connect/etc/asperaweb_id_dsa.openssh" $tempname ../data/prostate_cancer_data/mass_spec_data/
done

