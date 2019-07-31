library(dplyr)
library(magrittr)
library(ggplot2)
library(mzID)

# data munging for e coli

# read in data
ecoli1 <- mzID('data/ecoli_data/F020490.mzid')
ecoli2 <- mzID('data/ecoli_data/F020549.mzid')
#ecoli3 <- mzID('data/ecoli_data/F023569.mzid') this dataset has a offline prefractiontation step, we only want the online chromatography for RT prediction purposes

# convert to dataframe
ecoli_f_all1 <- flatten(ecoli1)
ecoli_f_all2 <- flatten(ecoli2)
# ecoli_f_all3 <- flatten(ecoli3)

# combine all ms files
ecoli_f_all <- rbind(ecoli_f_all1, ecoli_f_all2)

# filter out peptides that did not meet threshold
ecoli_f <- ecoli_f_all %>% filter(passthreshold == TRUE)

# confirm that pass threshold refers to the peptide spectrum match column
ecoli_f_all %>% ggplot(aes(`mascot:score`, fill = passthreshold)) + geom_histogram(bins = 1000, alpha = 0.3)

# get all the retention times and peptide sequences
ecoli_rt <- ecoli_f$`retention time(s)` %>% as.numeric()
ecoli_peps <- ecoli_f$pepseq

# combine into one dataframe
ecoli_df <- data.frame(peptide = ecoli_peps, rts = ecoli_rt)

# randomly sample data because there are too many peptides found to train the model

ecoli_df_sub <- sample_n(ecoli_df, 1000)

# write into a table
write.table(ecoli_df_sub, file = 'data/ecoli_data/ecoli_peptides_formatted.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(ecoli_df, file = 'data/ecoli_data/ecoli_peptides_formatted_all.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)


