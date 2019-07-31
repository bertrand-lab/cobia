
#!/bin/bash
# Database searching

DIR='../data/aylward_data/'
for FILE in "$DIR"*[[:digit:]].mzML
do
        echo "Processing $FILE file..."
        temp_string=${FILE/.mzML/}
        echo $temp_string
        NoiseFilterSGolay -in $FILE -out $temp_string'_SG.mzML'
        BaselineFilter -in $temp_string'_SG.mzML' -out $temp_string'_BF.mzML'
        PeakPickerHiRes -in $temp_string'_BF.mzML' -out $temp_string'_PP.mzML'
#	CompNovoCID -in $temp_string'_PP.mzML' -out $temp_string'_DN.idXML' &
        TextExporter -in $temp_string'_DN.idXML' -out $temp_string'_DN.csv' -id:peptides_only -id:first_dim_rt -separator ','
        MSGFPlusAdapter -in $temp_string'.mzML' -executable /software/MSGFPlus/MSGFPlus.jar -database ../data/aylward_data/combined_metagenomic_aylward.fasta -out $temp_string'.idXML' -add_decoys -threads 32
        PeptideIndexer -in $temp_string.idXML -fasta ../data/aylward_data/combined_metagenomic_aylward.revCat.fasta -out $temp_string'_PI.idXML' -threads 32 -decoy_string 'XXX_' -enzyme:specificity 'none'
        FalseDiscoveryRate -in $temp_string'_PI.idXML' -out $temp_string'_FDR.idXML' -FDR:PSM 0.05 -threads 32 -PSM 'true'
#       IDFilter -in $temp_string'_FDR.idXML' -out $temp_string'_IDF.idXML' -threads 32
done

