This details how to run all of the 'testing-scripts', which were used for the model validation and the case study.

######## CASE STUDY SCRIPTS

case_study.sh
    - Uses data from the Antarctic metatranscriptome from the manuscript (Bertrand et al 2013), and then predicts cofragmentation scores for targeted peptides from Fragilariopsis cylindrus

case_study_get_group_peps.sh
    - this bash script runs the two R scripts below to find taxon- and KOG class-specific peptides from the metatranscriptome from Bertrand et al 2013

case_study_get_tax_peps.R 
    - this scripts goes into taxonomic groups from the metatranscriptome to determine if a peptide is taxon specific

case_study_get_group_peps.R 
    - this script goes into the metatranscriptome and finds KOG family unique peptides, and returns a csv file

case_study_data_munging.R
    - this script is for examining the Frag-specific peptides and doing that data processing

case_study_tax_func_bias.R
    - this script does the KL divergence analysis, inputs require the case_study_get_tax_peps.R and case_study_get_group_peps.R scripts outputs

######## MS1 DATA COLLECTING AND PROCESSING SCRIPTS

convert_raw_to_mzml.sh
    - Goes through all the raw files downloaded for the Space Mouse dataset and the e coli dataset and converts them to mzML using ThermoRawFileParser

download_ms_data_ecoli.sh
    - script that goes through the FTP from PRIDE and downloads the raw files for the e coli dataset

download_ms_data_prostate_cancer.sh
    - script that goes through the FTP from PRIDE and downloads the raw files for the prostate cancer dataset

download_ms_data_space_mouse.sh
    - script that goes through the FTP from PRIDE and downloads the raw files for the space mouse dataset

get_ms1_dists.sh 
    - gets the MS1 values for Table 1 in the manuscript using a python script (this gets values for the Kleiner data)
get_ms1_dists_aylward.sh
    - script to go through each of the raw mzML files from Aylward metaproteomic paper, and print out the ms1 and the file name

get_ms1_dists_ecoli.sh
    - script to go through each of the raw mzML files for the e coli dataset, and print out the ms1 and the file name
get_ms1_dists_kleiner.sh
    - script to go through each of the raw mzML files for the kleiner dataset, and print out the ms1 and the file name

kleiner_idxml_csv_conversion.sh
    - converting the idXML files from the Kleiner database searching into csvs

process_all_ms1.sh
    - goes through the ms1 data and formats the output files so they can be more easily read

get_ms1_dist.py
    - python script that takes on argument (the mzML file), loads the file, gets the total number of ms1 peaks, 
and the total number of ms1 spectra

######## MAIN VALIDATION APPROACH

validation_database-searching-aylward-openms.sh
    - Mass spec preprocessing (smoothing, baseline filter, peak picking) and database searching of the Aylward dataset

validation_database-searching-kleiner-260-openms.sh
    - Mass spec preprocessing (smoothing, baseline filter, peak picking) and database searching of the Kleiner dataset

svm_models_aylward.sh
    - Fit four retention time models (different kernels) for the Aylward dataset

svm_models_ecoli.sh
    - Filters out reverse sequences from the E. coli database file with a python script
    - Digest the protein sequences with cobia
    - Fit a retention time model (RTModel)
    - Predict the retention times of all the digested database

filter_fasta.py
    - filters out sequences from the e coli database that are reversed (from DB search)

svm_models_kleiner.sh
    - Fit four retention time models (different kernels) for the Kleiner dataset

svm_models_prostate_cancer.sh
    - Fit a retention time model and predict a trypsin digested database using RTModel and cobia for prostate cancer dataset

svm_models_space_mouse.sh
    - Fit a retention time model and predict a trypsin digested database using RTModel and cobia for space mouse data set

validation_testing_aylward.sh
   - Take the RTPredict model output from different kernels, and then run the cofragmentation prediction for the Aylward dataset

validation_testing_broberg.sh
   - Runs the peptide retention time prediction and modification with BioLCCC
   - Predicts cofragmentation scores from these mod/rt predictions

validation_testing_ecoli.sh
   - Inputs the retention time prediction from RTModel and modifies the peptides 
   - Predicts cofragmentation scores

validation_testing_kleiner.sh
   - Take the RTPredict model output from different kernels, and then run the cofragmentation prediction for the two Kleiner datasets

validation_testing_space_mouse.sh
   - Inputs the retention time prediction from RTModel and modifies the peptides
   - Predicts cofragmentation scores

validation_checking_all.R
    - This script is rather large, it takes in all the validation*.sh scripts from above (all the output), and then processes it. The processing requires comparing cofragmentation scores from observed and not observed peptides, fitting a GLM to determine if the cofragmentation score has explanatory value, and then making the figures associated with this.

validation_inferred_peptides.R
    - This is specifically for the second type of validation, where we only look at inferred peptides. This script requires validation_checking_all.R to be run first, so that the functions from that script are available. It also produces a figure.

####### Misc. scripts

model-of-ion-peak-width.R
    - This is a linear regression from a paper that we used to estimate the ion peak width from column characteristics

validation_broberg_fractionation_plots.R
    - The Broberg dataset prefractionated their samples into four. We were thinking of then doing four separate cofragmentation predictions, but wanted to get a sense of their separation. This script looks at the four different fractions, and which peptides were identified, and their hydrophobicity. Most peptides were in the first fraction, and there was a wide range of hydrophobicities, so it didn't make sense to calculate cofragmentation scores for every fraction.

plotting-schematic.R
    - makes the schematic that was used for the TOC graphic and also for the first figure

data_munging_space_mouse.R
    - formatting the space mouse dataset so it can be easily read for RTModel and RTPredict

data_munging_prostate_cancer.sh
    - formatting the prostate cancer dataset so it can be easily read for RTModel and RTPredict

data_munging_ecoli.R
    - formatting the ecoli cancer dataset so it can be easily read for RTModel and RTPredict

plotting-retention-times.R
    - makes the plots for the supplementary material that looks at the distribution of retention times

