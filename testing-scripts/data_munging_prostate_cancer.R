# data munging for single organism proteomics

library(dplyr)
library(magrittr)

## go through each dataset 

## format the peptides and retention times to be a txt file tab separated with no headers

pro_cancer <- read.csv("data/prostate_cancer_data/Peptide list 16012017.csv", skip = 2)

pro_cancer_seq <- pro_cancer$Sequence
pro_cancer_sec <- pro_cancer$Retention.time..min.*60

df1 <- data.frame(as.factor(pro_cancer_seq), pro_cancer_sec)

df1_sub <- sample_n(df1, 1000)

write.table(df1_sub, file = 'data/prostate_cancer_data/prostate_cancer_peptides_formatted.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(df1, file = 'data/prostate_cancer_data/prostate_cancer_peptides_formatted_all.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)