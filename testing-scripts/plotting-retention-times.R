# plotting the distribution of retention times

library(ggplot2)
library(dplyr)
library(magrittr)
library(gridExtra)

broberg_mt_lc <- read.csv("data/broberg_data/broberg_mt_tryptic_peptide_lc-retention-times.csv")
broberg_mg_lc <- read.csv("data/broberg_data/broberg_mg_tryptic_peptide_lc-retention-times.csv")
broberg_mg_oak_lc <- read.csv("data/broberg_data/broberg_mg_qrob_oak_genome_tryptic_peptide_lc-retention-times.csv")

kleiner_run1_lc <- read.csv("data/kleiner_data/Run1_C4_832ng_lc_retention_times.csv_lc-retention-times.csv")
kleiner_run4_lc <- read.csv("data/kleiner_data/Run4_C4_2000ng_oligo_lc-retention-times.csv")
aylward_lc <- read.csv("data/aylward_data/aylward_tryptic_peptide_oligo_lc-retention-times.csv")

broberg_mt_lc_p <- broberg_mt_lc %>% ggplot(aes(rts)) + geom_histogram() + 
  theme_bw() + ylab('') + xlab('Predicted Retention Time') + ggtitle('D) Diseased Oak Tree Metatranscriptome') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

broberg_mg_lc_p <- broberg_mg_lc %>% ggplot(aes(rts)) + geom_histogram() + 
  theme_bw() + ylab('Count') + xlab('') + ggtitle('C) Diseased Oak Tree Metagenome') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

broberg_mg_oak_lc_p <- broberg_mg_oak_lc %>% ggplot(aes(rts)) + geom_histogram() + 
  theme_bw() + ylab('Count') + xlab('') + ggtitle('E) Diseased Oak Tree Genome + Metagenome') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

kleiner_run1_lc_p <- kleiner_run1_lc %>% ggplot(aes(rts/60)) + geom_histogram() + 
  theme_bw() + ylab('') + xlab('') + ggtitle('B) Mock Community (LC Run = 460 mins)') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

kleiner_run4_lc_p <- kleiner_run4_lc %>% ggplot(aes(rts/60)) + geom_histogram() + 
  theme_bw() + ylab('Count') + xlab('') + ggtitle('A) Mock Community (LC Run = 260 mins)') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

aylward_lc_p <- aylward_lc %>% ggplot(aes(rts/60)) + geom_histogram() + 
  theme_bw() + ylab('Count') + xlab('Predicted Retention Time') + ggtitle('F) Ant Fungus Garden Metagenome') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5))

jpeg("retention_times_plots.jpeg", width=170, height=210, units="mm", res=850)

grid.arrange(kleiner_run4_lc_p, kleiner_run1_lc_p, broberg_mg_lc_p, broberg_mt_lc_p, broberg_mg_oak_lc_p, aylward_lc_p, nrow = 3)

dev.off()


