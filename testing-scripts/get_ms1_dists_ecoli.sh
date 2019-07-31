#script to go through each of the raw mzML files from Aylward metaproteomic paper, and print out the ms1 and the
# file name

for FILE in ls ../data/ecoli_data/mass_spec_data/*mzML; do
     echo $FILE
     python get_ms1_dist.py $FILE
done

