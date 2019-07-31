library(dplyr)
library(magrittr)
library(ggplot2)

# data munging for fungal pathogen data

fungal_pathogen <- read.csv("data/fungal_pathogen_data/pr5b01054_si_001_peptides_retention_times.csv")

dim(fungal_pathogen)

fungal_pathogen_re <- data.frame(peptides = fungal_pathogen$sequence, rt_time = fungal_pathogen$scan_time*60)

fungal_pathogen_re_sub <- sample_n(fungal_pathogen_re, 1000)

write.table(fungal_pathogen_re_sub, file = 'data/fungal_pathogen_data/fungal_pathogen_peptides_formatted.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
