
#!/bin/bash
# Database searching

DIR='../data/kleiner_data/'
for FILE in "$DIR"*ng.mzML
do
        echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}
        echo $temp_string
        NoiseFilterSGolay -in $FILE -out $temp_string'_SG.mzML'
        BaselineFilter -in $temp_string'_SG.mzML' -out $temp_string'_BF.mzML'
        PeakPickerHiRes -in $temp_string'_BF.mzML' -out $temp_string'_PP.mzML'
        MSGFPlusAdapter -in $temp_string'.mzML' -executable /software/MSGFPlus/MSGFPlus.jar -database ../data/kleiner_data/Mock_Comm_RefDB_V3.fasta -out $temp_string'.idXML' -add_decoys -threads 32
        PeptideIndexer -in $temp_string.idXML -fasta ../data/kleiner_data/Mock_Comm_RefDB_V3.revCat.fasta -out $temp_string'_PI.idXML' -threads 32 -decoy_string 'XXX_' -enzyme:specificity 'none'
        FalseDiscoveryRate -in $temp_string'_PI.idXML' -out $temp_string'_FDR.idXML' -FDR:PSM 0.01 -threads 32 -PSM 'true'
#       IDFilter -in $temp_string'_FDR.idXML' -out $temp_string'_IDF.idXML' -threads 32
done

