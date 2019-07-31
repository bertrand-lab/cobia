declare -a StringArray=("../data/space_mouse_data/mass_spec_data/" "../data/ecoli_data/mass_spec_data/")

# for each directory
for val in ${StringArray[@]}; do
    echo $val

    for FILE in $val*raw; do
        echo $FILE
        mono /var/www/sfolder/general/ThermoRawFileParser/ThermoRawFileParser.exe -i $FILE -o $val -f 1
    done

done
