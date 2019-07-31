

for word in $(cat ../data/space_mouse_data/ms_data_file_names_space_mouse_clean_final.txt);
    do echo $word;
    tempname='prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2017/08/PXD005102/'$word
    echo $tempname
    ../../../.aspera/connect/bin/ascp -TQ -l200m -P 33001 -i "../../../.aspera/connect/etc/asperaweb_id_dsa.openssh" $tempname ../data/space_mouse_data/mass_spec_data/
done

