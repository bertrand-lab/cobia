for FILE in ls ../data/ms1_data/*txt; do
    echo $FILE
    tempfile=${FILE/.txt}
    echo $tempfile
    sed 's/)//g' $FILE | sed 's/(//g' | sed 's/,//g' > $tempfile'_formatted.txt'
done
