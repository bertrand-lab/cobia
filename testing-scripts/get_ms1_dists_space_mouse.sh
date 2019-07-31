
for FILE in ls ../data/space_mouse_data/mass_spec_data/*mzML; do
     echo $FILE
     python get_ms1_dist.py $FILE
done

