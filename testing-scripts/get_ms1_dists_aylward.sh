# script to go through each of the raw mzML files from Aylward metaproteomic paper, and print out the ms1 and the file name
DIR='../data/aylward_data/'

for i in ls $DIR*[0-9].mzML
do
    echo $i
    python get_ms1_dist.py $i
done

#python get_ms1_dist.py ../data/kleiner_data/Run1_C4_832ng.mzML
#python get_ms1_dist.py ../data/kleiner_data/Run4_C4_2000ng.mzML
