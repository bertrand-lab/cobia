# Simple case - running BioLCCC prediction
bash validation_testing_broberg.sh

# RTModel Case - 
## database searching
bash validation_database-searching-aylward-openms.sh
bash validation_database-searching-kleiner-260-openms.sh

## training SVM
bash svm_models_aylward.sh
bash svm_models_kleiner.sh

## predicting cofragmentation
bash validation_testing_aylward.sh
bash validation_testing_kleiner.sh
