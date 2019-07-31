# all validation collapsed into one script

# loading libraries

library("mzR")
library("mzID")
library("MSnID")
library("MSnbase")
library(magrittr)
library(ggplot2)
library(dplyr)
library(stringi)
library(broom)
library(readxl)
library(seqinr)
library(gridExtra)

# data processing scripts

# formatting cofragmentation output for validation
cofrag_processor <- function(cofrag_output_df, observed_peps){
  
  # cofrag_output_df <- aylward_rt_cofrag_linear
  # observed_peps <- aylward_prot_data_orbi$pep_seq
  
  # removing peptide identifiers that are common across contigs
  cofrag_output_df$X <- NULL
  cofrag_output_df$'Unnamed..0' <- NULL
  cofrag_output_df$contig <- NULL
  
  # removing the termini ModX formatting
  cofrag_output_df$pep_seq <- gsub(pattern = '-OH', replacement = "", x = cofrag_output_df$peptide_sequence, fixed = TRUE) 
  cofrag_output_df$pep_seq <- gsub(pattern = 'Ac-', replacement = "", x = cofrag_output_df$pep_seq, fixed = TRUE) 
  
  cofrag_output_df$rescale_cofrag_score <- scales::rescale(cofrag_output_df$mean_cofrag_score, to = c(0, 100))
  
  # multiple instances of the same peptide is within the dataframe, removing those by taking only columns with everything in it
  cofrag_output_df2 <- cofrag_output_df[cofrag_output_df %>% complete.cases(), ]
  
  # making a a column that identifies if the peptide was observed
  cofrag_output_df2$prote_d <- cofrag_output_df2$pep_seq %in% observed_peps
  
  return(cofrag_output_df2)
}

# functino for number annotation
number_of_unseen <- function(cofrag_output_file){
  temp_df <- cofrag_output_file[!is.na(cofrag_output_file$mean_cofrag_score), ]
  number_observed <- temp_df$prote_d %>% sum()
  number_not_observed <- nrow(temp_df) - number_observed
  return(list(number_observed, number_not_observed))
} 

'%!in%' <- function(x,y)!('%in%'(x,y))

extract_coef <- function(glm_model, char_id){
  glm_model_summary <- summary(glm_model)
  glm_coef <- coefficients(glm_model_summary) %>% as.data.frame()
  X0_coef <- glm_coef$Estimate[2]
  X0_se <- glm_coef$`Std. Error`[2]
  return(data.frame(coef = X0_coef, se = X0_se, char_id = rep(char_id, 1)))
}
extract_coef2 <- function(glm_model, char_id){
  glm_model_summary <- summary(glm_model)
  glm_coef <- coefficients(glm_model_summary) %>% as.data.frame()
  dataset_name <- rep(char_id, nrow(glm_coef))
  glm_coef$dataset_name <- dataset_name
  glm_coef$Variable <- rownames(glm_coef)
  rownames(glm_coef) <- NULL
  # X0_coef <- glm_coef$Estimate[2]
  # X0_se <- glm_coef$`Std. Error`[2]
  # return(data.frame(coef = X0_coef, se = X0_se, char_id = rep(char_id, 1)))
  return(glm_coef)
}

## single organism proteomics

# space mouse

space_mouse_observed <- read.table("data/space_mouse_data/space_mouse_peptides_formatted_all.txt", 
                                   col.names = c('sequence', 'retention_time'))
space_mouse_cofrag_oligo <- read.csv("data/space_mouse_data/UP000000589_10090_rt_oligo_cofrag_mi-0.00833333_ipw-0.95_para-15_co-sim.csv")

space_mouse_oligo2 <- cofrag_processor(cofrag_output_df = space_mouse_cofrag_oligo, observed_peps = space_mouse_observed$sequence)

length(space_mouse_oligo2$prote_d)

number_of_unseen(space_mouse_oligo2)
space_mouse_n <- "italic(N)[italic(unobserved)] == 701591"
space_mouse_n_sub <- "italic(N)[italic(observed)] == 9801"
space_mouse_oligo2_p8 <- space_mouse_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.12, "(0.009)"), y = 1e8, x = 13, parse = TRUE, size = 7) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  annotate("text", x = 4.5, y = Inf, vjust = 3, label = space_mouse_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 4.5, y = Inf, vjust = 4.5, label = space_mouse_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  ggtitle('B) Space Mouse') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none", title = element_text(size = 7.5));space_mouse_oligo2_p8

glm_space_mouse <- glm(data = space_mouse_oligo2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_space_mouse_rescale <- glm(data = space_mouse_oligo2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_space_mouse_exp <- glm(data = space_mouse_oligo2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")

# prostate cancer

prostate_cancer_observed <- read.table("data/prostate_cancer_data/prostate_cancer_peptides_formatted_all.txt", 
                                   col.names = c('sequence', 'retention_time'))
prostate_cancer_cofrag_oligo <- read.csv("data/prostate_cancer_data/Human_uniprot_09082016_rt_oligo_cofrag_mi-0.00833333_ipw-0.722_para-15_co-sim.csv")

prostate_cancer_oligo2 <- cofrag_processor(cofrag_output_df = prostate_cancer_cofrag_oligo, observed_peps = prostate_cancer_observed$sequence)

length(prostate_cancer_oligo2$prote_d)

number_of_unseen(prostate_cancer_oligo2)
prostate_cancer_n <- "italic(N)[italic(unobserved)] == 696945"
prostate_cancer_n_sub <- "italic(N)[italic(observed)] == 903"
prostate_cancer_oligo2_p7 <- prostate_cancer_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.12, "(0.009)"), y = 1e8, x = 13, parse = TRUE, size = 7) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  annotate("text", x = 4.5, y = Inf, vjust = 3, label = prostate_cancer_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 4.5, y = Inf, vjust = 4.5, label = prostate_cancer_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  ggtitle('A) Prostate Cancer Biomarkers') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none", title = element_text(size = 7.5));prostate_cancer_oligo2_p7

glm_prostate_cancer <- glm(data = prostate_cancer_oligo2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_prostate_cancer_rescale <- glm(data = prostate_cancer_oligo2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_prostate_cancer_exp <- glm(data = prostate_cancer_oligo2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")


# e coli

ecoli_observed <- read.table("data/ecoli_data/ecoli_peptides_formatted_all.txt", 
                                       col.names = c('sequence', 'retention_time'))
ecoli_cofrag_oligo <- read.csv("data/ecoli_data/uniprot-ecoli.no-decoy_rt_oligo_cofrag_mi-0.00833333_ipw-0.9482396999999999_para-15_co-sim.csv")

ecoli_oligo2 <- cofrag_processor(cofrag_output_df = ecoli_cofrag_oligo, observed_peps = ecoli_observed$sequence)

length(ecoli_oligo2$prote_d)

number_of_unseen(ecoli_oligo2)
ecoli_n <- "italic(N)[italic(unobserved)] == 87162"
ecoli_n_sub <- "italic(N)[italic(observed)] == 13296"
ecoli_oligo2_p9 <- ecoli_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 0.1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.12, "(0.009)"), y = 1e8, x = 13, parse = TRUE, size = 7) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  annotate("text", x = 2.35, y = Inf, vjust = 3, label = ecoli_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 2.35, y = Inf, vjust = 4.5, label = ecoli_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  ggtitle('C) E. coli') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none", title = element_text(size = 7.5));ecoli_oligo2_p9

glm_ecoli <- glm(data = ecoli_oligo2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_ecoli_rescale <- glm(data = ecoli_oligo2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_ecoli_exp <- glm(data = ecoli_oligo2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")

# Kleiner Run 1 Dataset ---------------------------------------------------

# openms idXML file filtered

run1_oms <- read.csv("data/kleiner_data/Run1_C4_832ng_IDF.csv")

# reading in different cofrag scores from different lc predictions
kleiner_cofrag_460_rbf <- read.csv("data/kleiner_data/Run1_C4_832ng_rbf_cofrag_mi-0.00833333_ipw-2.75_para-15_co-sim.csv")
kleiner_cofrag_460_oligo <-  read.csv("data/kleiner_data/Run1_C4_832ng_oligo_cofrag_mi-0.00833333_ipw-2.75_para-15_co-sim.csv")
kleiner_cofrag_460_linear <-  read.csv("data/kleiner_data/Run1_C4_832ng_linear_cofrag_mi-0.00833333_ipw-2.75_para-15_co-sim.csv")

# removing the modX modifications to N termini and C termini
kleiner_cofrag_460_rbf2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_460_rbf, observed_peps = run1_oms$sequence)
kleiner_cofrag_460_oligo2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_460_oligo, observed_peps = run1_oms$sequence)
kleiner_cofrag_460_linear2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_460_linear, observed_peps = run1_oms$sequence)

number_of_unseen(kleiner_cofrag_460_oligo2)
kleiner460_n <- "italic(N)[italic(unobserved)] == 2201822"
kleiner460_n_sub <- "italic(N)[italic(observed)] == 5921"
kleiner_460_p2 <- kleiner_cofrag_460_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.12, "(0.009)"), y = 1e8, x = 13, parse = TRUE, size = 7) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  annotate("text", x = 6, y = Inf, vjust = 3, label = kleiner460_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 6, y = Inf, vjust = 4.5, label = kleiner460_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  ggtitle('B) Mock Community (LC Run = 460 mins)') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none", title = element_text(size = 7.5));kleiner_460_p2

glm_kleiner460_rbf <- glm(data = kleiner_cofrag_460_rbf2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_kleiner460_rbf_rescale <- glm(data = kleiner_cofrag_460_rbf2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_kleiner460_rbf_exp <- glm(data = kleiner_cofrag_460_rbf2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")

glm_kleiner460_oligo <- glm(data = kleiner_cofrag_460_oligo2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_kleiner460_oligo_rescale <- glm(data = kleiner_cofrag_460_oligo2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_kleiner460_oligo_exp <- glm(data = kleiner_cofrag_460_oligo2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")

glm_kleiner460_linear <- glm(data = kleiner_cofrag_460_linear2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_kleiner460_linear_rescale <- glm(data = kleiner_cofrag_460_linear2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_kleiner460_linear_exp <- glm(data = kleiner_cofrag_460_linear2, prote_d ~ mean_cofrag_score + rts + mz, family = "binomial")

glm_kleiner460_rbf_summary <- summary(glm_kleiner460_rbf)
glm_kleiner460_rbf_rescale_summary <- summary(glm_kleiner460_rbf_rescale)
glm_kleiner460_rbf_exp_summary <- summary(glm_kleiner460_rbf_exp)

glm_kleiner460_oligo_summary <- summary(glm_kleiner460_oligo)
glm_kleiner460_oligo_rescale_summary <- summary(glm_kleiner460_oligo_rescale)
glm_kleiner460_oligo_exp_summary <- summary(glm_kleiner460_oligo_exp)

glm_kleiner460_linear_summary <- summary(glm_kleiner460_linear)
glm_kleiner460_linear_rescale_summary <- summary(glm_kleiner460_linear_rescale)
glm_kleiner460_linear_exp_summary <- summary(glm_kleiner460_linear_exp)

# kleiner run 4 -----------------------------------------------------------

# peptides detected using OpenMS pipeline
run4_oms <- read.csv("data/kleiner_data/Run4_C4_2000ng_IDF.csv")

kleiner_cofrag_260_rbf <- read.csv("data/kleiner_data/Run4_C4_2000ng_rbf_cofrag_mi-0.00833333_ipw-1.44_para-15_co-sim.csv")
kleiner_cofrag_260_oligo <- read.csv("data/kleiner_data/Run4_C4_2000ng_oligo_cofrag_mi-0.00833333_ipw-1.44_para-15_co-sim.csv")
kleiner_cofrag_260_linear <- read.csv("data/kleiner_data/Run4_C4_2000ng_linear_cofrag_mi-0.00833333_ipw-1.44_para-15_co-sim.csv")

kleiner_cofrag_260_rbf2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_260_rbf, observed_peps = run4_oms$sequence)
kleiner_cofrag_260_oligo2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_260_oligo, observed_peps = run4_oms$sequence)
kleiner_cofrag_260_linear2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_260_linear, observed_peps = run4_oms$sequence)

# modelling whether cofragmentation risk has an influence on peptide presence/absence
glm_kleiner260_rbf <- glm(data = kleiner_cofrag_260_rbf2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_kleiner260_rbf_rescale <- glm(data = kleiner_cofrag_260_rbf2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_kleiner260_rbf_exp <- glm(data = kleiner_cofrag_260_rbf2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

glm_kleiner260_oligo <- glm(data = kleiner_cofrag_260_oligo2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_kleiner260_oligo_rescale <- glm(data = kleiner_cofrag_260_oligo2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_kleiner260_oligo_exp <- glm(data = kleiner_cofrag_260_oligo2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

glm_kleiner260_linear <- glm(data = kleiner_cofrag_260_linear2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_kleiner260_linear_rescale <- glm(data = kleiner_cofrag_260_linear2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_kleiner260_linear_exp <- glm(data = kleiner_cofrag_260_linear2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

glm_kleiner260_rbf_summary <- summary(glm_kleiner260_rbf)
glm_kleiner260_rbf_rescale_summary <- summary(glm_kleiner260_rbf_rescale)
glm_kleiner260_rbf_exp_summary <- summary(glm_kleiner260_rbf_exp)

glm_kleiner260_oligo_summary <- summary(glm_kleiner260_oligo)
glm_kleiner260_oligo_rescale_summary <- summary(glm_kleiner260_oligo_rescale)
glm_kleiner260_oligo_exp_summary <- summary(glm_kleiner260_oligo_exp)

glm_kleiner260_linear_summary <- summary(glm_kleiner260_linear)
glm_kleiner260_linear_rescale_summary <- summary(glm_kleiner260_linear_rescale)
glm_kleiner260_linear_exp_summary <- summary(glm_kleiner260_linear_exp)

number_of_unseen(kleiner_cofrag_260_oligo2)
kleiner260_n <- "italic(N)[italic(unobserved)] == 2196645"
kleiner260_n_sub <- "italic(N)[italic(observed)] == 11098"
kleiner_260_p1 <- kleiner_cofrag_260_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 0.2) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('log(Count)') +
  # annotate("text", label = paste("beta == ", -0.071, "(0.013)"), 
  # y = 1e9, x = 7, parse = TRUE, size = 6) +
  # xlab('Predicted # Cofragmenting Peptides') +
  annotate("text", x = 7, y = Inf, vjust = 3, label = kleiner260_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 7, y = Inf, vjust = 4.5, label = kleiner260_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  xlab('') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  ggtitle('A) Mock Community (LC Run = 260 mins)') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));kleiner_260_p1

# Broberg MT and MG  ------------------------------------------------------

# read in broberg MS data
bro_ms <- readxl::read_xlsx("data/broberg_data/Supplementary_Table_33.xlsx", skip = 4)

bro_co_mt <- read.csv("data/broberg_data/broberg_mt_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")
bro_co_mg <- read.csv("data/broberg_data/broberg_mg_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")
bro_co_mg_oak <- read.csv("data/broberg_data/broberg_mg_qrob_oak_genome_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")

# function to drop X, 'Unnamed..O', and 'Contig', adjust the peptide sequence, and just get unique peptides as a df

# formatting the peptide file (Supplementary File 33)
bro_peps <- bro_ms[!grepl(pattern = 'Sequence', bro_ms$Sequence), ]
bro_peps_no_numeric <- bro_peps[!grepl(pattern = '[0-9]', bro_peps$Sequence), ]

bro_peps2 <- dplyr::rename(bro_peps_no_numeric, iso_int = 'Isolation Interference [%]', rt = 'RT [min]', pep_seq = "Sequence")
bro_peps2$iso_int <- as.numeric(bro_peps2$iso_int)
bro_peps2$rt <- as.numeric(bro_peps2$rt)
bro_peps2$pep_seq <- as.character(bro_peps2$pep_seq)
bro_peps2$pep_seq <- toupper(bro_peps2$pep_seq)

bro_peps2a <- bro_peps2[!is.na(bro_peps2$AT4),]
bro_peps2a$AT4 <- bro_peps2a$AT4 %>% as.numeric()

bro_peps3 <- aggregate(data = bro_peps2a, iso_int ~ pep_seq, FUN = mean, na.rm = TRUE)
bro_peps3_rt <- aggregate(data = bro_peps2a, rt ~ pep_seq, FUN = mean, na.rm = TRUE)
bro_peps3_abun <- aggregate(data = bro_peps2a, AT4 ~ pep_seq, FUN = mean, na.rm = TRUE)
bro_peps4 <- data.frame(bro_peps3, bro_peps3_rt$rt, bro_peps3_abun$AT4)
bro_peps4$pep_seq <- toupper(bro_peps4$pep_seq)

bro_co_mg2 <- cofrag_processor(cofrag_output_df = bro_co_mg, observed_peps = bro_peps4$pep_seq)
bro_co_mt2 <- cofrag_processor(cofrag_output_df = bro_co_mt, observed_peps = bro_peps4$pep_seq)
bro_co_mg_oak2 <- cofrag_processor(cofrag_output_df = bro_co_mg_oak, observed_peps = bro_peps4$pep_seq)

number_of_unseen(bro_co_mg2)
broberg_mg_n <- "italic(N)[italic(unobserved)] == 107815"
broberg_mg_n_sub <- "italic(N)[italic(observed)] == 178"
broberg_mg_p3 <- bro_co_mg2 %>% 
  # filter(peplen > 8) %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() +
  theme_bw() +
  ylab('log(Count)') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  ggtitle('C) Diseased Oak Tree Metagenome') +
  annotate("text", x = 80, y = Inf, vjust = 3, label = broberg_mg_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 80, y = Inf, vjust = 4.5, label = broberg_mg_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  # annotate("text", label = paste("beta == ", -0.054, "(0.01)"), 
  # y = 1e6, x = 35, parse = TRUE, size = 6) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));broberg_mg_p3

number_of_unseen(bro_co_mt2)
broberg_mt_n <- "italic(N)[italic(unobserved)] == 341461"
broberg_mt_n_sub <- "italic(N)[italic(observed)] == 473"
broberg_mt_p4 <- bro_co_mt2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  ggtitle('D) Diseased Oak Tree Metatranscriptome') +
  annotate("text", x = 250, y = Inf, vjust = 3, label = broberg_mt_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 250, y = Inf, vjust = 4.5, label = broberg_mt_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  # annotate("text", label = paste("beta == ", -0.021, "(0.0023)"), 
  # y = 1e7, x = 80, parse = TRUE, size = 6) +
  xlab('') + theme(legend.position = "none", title = element_text(size = 7.5));broberg_mt_p4
# geom_density(alpha = 0.4)


number_of_unseen(bro_co_mg_oak2)
broberg_mg_oak_n <- "italic(N)[italic(unobserved)] == 682956"
broberg_mg_oak_n_sub <- "italic(N)[italic(observed)] == 439"
broberg_mg_p3_oak <- bro_co_mg_oak2 %>% 
  # filter(peplen > 8) %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('log(Count)') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  ggtitle('E) Diseased Oak Tree Genome + Metagenome') +
  annotate("text", x = 375, y = Inf, vjust = 3, label = broberg_mg_oak_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 375, y = Inf, vjust = 4.5, label = broberg_mg_oak_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  xlab('Cofragmentation Score') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));broberg_mg_p3_oak


# glm_broberg_mg <- glm(data = bro_mt_joined, prote_d ~ mean_cofrag_score, family = 'binomial')

glm_broberg_mt <- glm(data = bro_co_mt2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mt_rescale <- glm(data = bro_co_mt2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_broberg_mt_exp <- glm(data = bro_co_mt2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

glm_broberg_mg <- glm(data = bro_co_mg2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mg_rescale <- glm(data = bro_co_mg2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_broberg_mg_exp <- glm(data = bro_co_mg2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

glm_broberg_mg_oak <- glm(data = bro_co_mg_oak2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mg_oak_rescale <- glm(data = bro_co_mg_oak2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_broberg_mg_oak_exp <- glm(data = bro_co_mg_oak2, prote_d ~ mean_cofrag_score + rts + mz, family = 'binomial')

summary(glm_broberg_mt)
summary(glm_broberg_mg)
summary(glm_broberg_mg_oak)
# summary(glm_broberg_mt)

# glm_mt_2 <- glm(data = bro_mt_joined, prote_d ~ mean_cofrag_score + mz + rts, family = 'binomial')
# glm_mg_2 <- glm(data = bro_mg_joined, prote_d ~ mean_cofrag_score + mz + rts, family = 'binomial')

# summary(glm_mt_2)
# summary(glm_mg_2)


# aylward  ----------------------------------------------------------------

# read in proteomic results (peptide list)
aylward_prot_data1 <- read_xls("data/aylward_data/ismej201210x5.xls", sheet = 1)
aylward_prot_data2 <- read_xls("data/aylward_data/ismej201210x5.xls", sheet = 2)
aylward_prot_data3 <- read_xls("data/aylward_data/ismej201210x5.xls", sheet = 3)
aylward_prot_data <- rbind(aylward_prot_data1, aylward_prot_data2, aylward_prot_data3)

test1 <- aylward_prot_data1[grepl(pattern = 'Orbi', x = aylward_prot_data1$Samples, ignore.case = TRUE),]$peptide
length(unique(test1))
test2 <- aylward_prot_data2[grepl(pattern = 'Orbi', x = aylward_prot_data2$Samples, ignore.case = TRUE),]$peptide
length(unique(test2))
test3 <- aylward_prot_data3[grepl(pattern = 'Orbi', x = aylward_prot_data3$Samples, ignore.case = TRUE),]$peptide
length(unique(test3))


# data formatting for peptide sequences from proteomic data
aylward_prot_data$peptide <- stri_sub(aylward_prot_data$peptide, 3, -3)
aylward_prot_data <- aylward_prot_data %>% dplyr::rename(pep_seq = peptide)
aylward_prot_data_orbi <- aylward_prot_data[grepl(pattern = 'Orbi', x = aylward_prot_data$Samples),]

#read in cofragmentation simulation output
aylward_rt_cofrag_linear <- read.csv("data/aylward_data/aylward_linear_cofrag_mi-0.00833333_ipw-0.341_para-20_co-sim.csv")
aylward_rt_cofrag_oligo <- read.csv("data/aylward_data/aylward_oligo_cofrag_mi-0.00833333_ipw-0.341_para-20_co-sim.csv")
aylward_rt_cofrag_rbf <- read.csv("data/aylward_data/aylward_rbf_cofrag_mi-0.00833333_ipw-0.341_para-20_co-sim.csv")

# removing the modX modifications to N termini and C termini
aylward_rt_cofrag_rbf2 <- cofrag_processor(cofrag_output_df = aylward_rt_cofrag_rbf, observed_peps = aylward_prot_data_orbi$pep_seq)
aylward_rt_cofrag_oligo2 <- cofrag_processor(cofrag_output_df = aylward_rt_cofrag_oligo, observed_peps = aylward_prot_data_orbi$pep_seq)
aylward_rt_cofrag_linear2 <- cofrag_processor(cofrag_output_df = aylward_rt_cofrag_linear, observed_peps = aylward_prot_data_orbi$pep_seq)

## testing something out

## aylward_rt_cofrag_linear2$pep_seq %!in% aylward_prot_data_orbi$pep_seq
## aylward_prot_data_orbi$pep_seq[aylward_prot_data_orbi$pep_seq %!in% aylward_rt_cofrag_linear2$pep_seq]
##

number_of_unseen(aylward_rt_cofrag_oligo2)
aylward_oligo_n <- "italic(N)[italic(unobserved)] == 2799709"
aylward_oligo_n_sub <- "italic(N)[italic(observed)] == 63"
aylward_p5_oligo <- aylward_rt_cofrag_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Cofragmentation Score') + 
  annotate("text", x = 97, y = Inf, vjust = 3, label = aylward_oligo_n, parse = TRUE, colour = 'grey62', size = 3) +
  annotate("text", x = 97, y = Inf, vjust = 4.5, label = aylward_oligo_n_sub, parse = TRUE, colour = 'grey24', size = 3) +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('E) Ant Fungus Garden Metagenome');aylward_p5_oligo

number_of_unseen(aylward_rt_cofrag_linear2)
aylward_linear_n <- "italic(N)[italic(unobserved)] == 2824883"
aylward_linear_n_sub <- "italic(N)[italic(observed)] == 63"
aylward_p5_linear <- aylward_rt_cofrag_linear2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Predicted # Cofragmenting Peptides') + 
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('E) Ant Fungus Garden Metagenome');aylward_p5_linear

number_of_unseen(aylward_rt_cofrag_rbf2)
aylward_oligo_n <- "italic(N)[italic(unobserved)] == 2824883"
aylward_oligo_n_sub <- "italic(N)[italic(observed)] == 63"
aylward_p5_rbf <- aylward_rt_cofrag_rbf2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Predicted # Cofragmenting Peptides') + 
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('E) Ant Fungus Garden Metagenome');aylward_p5_rbf


# glm.fit fitted probabilities numerically 0 or 1 occured
glm_aylward_oligo <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_oligo2)
glm_aylward_oligo_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_oligo2)
glm_aylward_oligo_exp <- glm(prote_d ~ mean_cofrag_score + rts + mz, family = "binomial", data = aylward_rt_cofrag_oligo2)

glm_aylward_linear <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_linear2)
glm_aylward_linear_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_linear2)
glm_aylward_linear_exp <- glm(prote_d ~ mean_cofrag_score + rts + mz, family = "binomial", data = aylward_rt_cofrag_linear2)

glm_aylward_rbf <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_rbf2)
glm_aylward_rbf_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_rbf2)
glm_aylward_rbf_exp <- glm(prote_d ~ mean_cofrag_score + rts + mz, family = "binomial", data = aylward_rt_cofrag_rbf2)

# aggregated plot ---------------------------------------------------------

# plotting aggregated validation results

# to do

# change font size of coef plot
# add in number of peps
# add ggtitle for plot F

# Making 


coef_diagram <- rbind(extract_coef(glm_kleiner260_oligo, 'Mock Community, LC 260'),
                      extract_coef(glm_kleiner460_oligo, 'Mock Community, LC 460'),
                      extract_coef(glm_broberg_mg, 'Diseased Oak Tree Metagenome'),
                      extract_coef(glm_broberg_mt, 'Diseased Oak Tree Metatranscriptome'),
                      # extract_coef(glm_broberg_mg_oak, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                      extract_coef(glm_aylward_oligo, 'Ant Fungus Garden'),
                      extract_coef(glm_prostate_cancer, 'Prostate Cancer Biomarkers'),
                      extract_coef(glm_space_mouse, 'Space Mouse'),
                      extract_coef(glm_ecoli, 'E. coli'))
coef_diagram$scaled <- rep('Cofragmentation Score', nrow(coef_diagram))

coef_diagram_rescale <- rbind(extract_coef(glm_kleiner260_oligo_rescale, 'Mock Community, LC 260'),
                      extract_coef(glm_kleiner460_oligo_rescale, 'Mock Community, LC 460'),
                      extract_coef(glm_broberg_mg_rescale, 'Diseased Oak Tree Metagenome'),
                      extract_coef(glm_broberg_mt_rescale, 'Diseased Oak Tree Metatranscriptome'),
                      # extract_coef(glm_broberg_mg_oak_rescale, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                      extract_coef(glm_aylward_oligo_rescale, 'Ant Fungus Garden'),
                      extract_coef(glm_prostate_cancer_rescale, 'Prostate Cancer Biomarkers'),
                      extract_coef(glm_space_mouse_rescale, 'Space Mouse'),
                      extract_coef(glm_ecoli_rescale, 'E. coli'))
coef_diagram_rescale$scaled <- rep('Scaled Cofragmentation Score', nrow(coef_diagram_rescale))

coef_diagram_pooled <- rbind(coef_diagram, coef_diagram_rescale)

coef_graph <- ggplot(coef_diagram_pooled, aes(x = char_id, y = coef)) + 
  geom_point(alpha = 0.5) +
  coord_flip() +
  geom_linerange(aes(ymin = coef - se, ymax = coef + se, x = char_id)) +
  geom_hline(yintercept = 0, lty = 2) + xlab("") + ylab("Coefficient Estimate") + 
  # ylim(c(-0.25, 0.25)) + 
  theme_bw() + 
  facet_grid(~scaled) +
  # ggtitle('F) Summaries of Validation by Dataset') + 
  theme(legend.position = "none") +
  scale_x_discrete(limits = rev(levels(coef_diagram$char_id)));coef_graph

dev.off()

jpeg("figures/main_validation_plots_metaproteome.jpeg", width=170, height=170, units="mm", res=850)

grid.arrange(kleiner_260_p1, kleiner_460_p2, broberg_mg_p3, broberg_mt_p4, 
             # broberg_mg_p3_oak, 
             aylward_p5_oligo, nrow = 3)

dev.off()

jpeg("figures/main_validation_plots_single_proteome.jpeg", width=170/2, height=170/5, units="mm", res=850)

grid.arrange(prostate_cancer_oligo2_p7, space_mouse_oligo2_p8, ecoli_oligo2_p9, nrow = 3)

dev.off()

ggsave(coef_graph, filename = "figures/main_coefficient_plot.jpeg", width = 170, height = 170/5, units = "mm")


# kernel plot -------------------------------------------------------------


### coefficient plot of different kernels for Aylward and Kleiner predictions

coef_diagram_supplement <- rbind(extract_coef(glm_kleiner260_oligo, 'Mock Community, LC 260, OLIGO kernel'),
                      extract_coef(glm_kleiner260_oligo_exp, 'Mock Community, LC 260, Additional Explanatory, OLIGO kernel'),
                      extract_coef(glm_kleiner260_oligo_rescale, 'Mock Community, LC 260, Rescaled, OLIGO kernel'),
                      extract_coef(glm_kleiner260_linear, 'Mock Community, LC 260, Linear kernel'),
                      extract_coef(glm_kleiner260_linear_exp, 'Mock Community, LC 260, Additional Explanatory, Linear kernel'),
                      extract_coef(glm_kleiner260_linear_rescale, 'Mock Community, LC 260, Rescaled, linear kernel'),
                      extract_coef(glm_kleiner260_rbf, 'Mock Community, LC 260, RBF kernel'),
                      extract_coef(glm_kleiner260_rbf_exp, 'Mock Community, LC 260, Additional Explanatory, RBF kernel'),
                      extract_coef(glm_kleiner260_rbf_rescale, 'Mock Community, LC 260, Rescaled, RBF kernel'),
                      extract_coef(glm_kleiner460_oligo, 'Mock Community, LC 460, OLIGO kernel'),
                      extract_coef(glm_kleiner460_oligo_exp, 'Mock Community, LC 460, Additional Explanatory, OLIGO kernel'),
                      extract_coef(glm_kleiner460_oligo_rescale, 'Mock Community, LC 460, Rescaled, OLIGO kernel'),
                      extract_coef(glm_kleiner460_linear, 'Mock Community, LC 460, Linear kernel'),
                      extract_coef(glm_kleiner460_linear_exp, 'Mock Community, LC 460, Additional Explanatory, Linear kernel'),
                      extract_coef(glm_kleiner460_linear_rescale, 'Mock Community, LC 460, Rescaled, linear kernel'),
                      extract_coef(glm_kleiner460_rbf, 'Mock Community, LC 460, RBF kernel'),
                      extract_coef(glm_kleiner460_rbf_exp, 'Mock Community, LC 460, Additional Explanatory, RBF kernel'),
                      extract_coef(glm_kleiner460_rbf_rescale, 'Mock Community, LC 460, Rescaled, RBF kernel'),
                      extract_coef(glm_broberg_mg, 'Diseased Oak Tree Metagenome'),
                      extract_coef(glm_broberg_mg_exp, 'Diseased Oak Tree Metagenome, Additional Explanatory'),
                      extract_coef(glm_broberg_mg_rescale, 'Diseased Oak Tree Metagenome, Rescaled'),
                      extract_coef(glm_broberg_mt, 'Diseased Oak Tree Metatranscriptome'),
                      extract_coef(glm_broberg_mt_exp, 'Diseased Oak Tree Metatranscriptome, Additional Explanatory'),
                      extract_coef(glm_broberg_mt_rescale, 'Diseased Oak Tree Metatranscriptome, Rescaled'),
                      # extract_coef(glm_broberg_mg_oak, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                      # extract_coef(glm_broberg_mg_oak_exp, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended, Additional Explanatory'),
                      # extract_coef(glm_broberg_mg_oak_rescale, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended, Rescaled'),
                      extract_coef(glm_aylward_oligo, 'Ant Fungus Garden, OLIGO kernel'),
                      extract_coef(glm_aylward_oligo_exp, 'Ant Fungus Garden, OLIGO kernel, Additional Explanatory'),
                      extract_coef(glm_aylward_oligo_rescale, 'Ant Fungus Garden, OLIGO Kernel, Rescaled'),
                      extract_coef(glm_aylward_linear, 'Ant Fungus Garden, Linear kernel'),
                      extract_coef(glm_aylward_linear_exp, 'Ant Fungus Garden, Linear kernel, Additional Explanatory'),
                      extract_coef(glm_aylward_linear_rescale, 'Ant Fungus Garden, Linear Kernel, Rescaled'),
                      extract_coef(glm_aylward_rbf, 'Ant Fungus Garden, RBF kernel'),
                      extract_coef(glm_aylward_rbf_exp, 'Ant Fungus Garden, RBF kernel, Additional Explanatory'),
                      extract_coef(glm_aylward_rbf_rescale, 'Ant Fungus Garden, RBF Kernel, Rescaled'),
                      extract_coef(glm_prostate_cancer, 'Prostate Cancer Biomarkers, OLIGO kernel'),
                      extract_coef(glm_prostate_cancer_exp, 'Prostate Cancer Biomarkers, OLIGO kernel, Additional Explanatory'),
                      extract_coef(glm_prostate_cancer_rescale, 'Prostate Cancer Biomarkers, OLIGO kernel, Rescaled'),
                      extract_coef(glm_space_mouse, 'Space Mouse, OLIGO kernel'),
                      extract_coef(glm_space_mouse_exp, 'Space Mouse, OLIGO kernel, Additional Explanatory'),
                      extract_coef(glm_space_mouse_rescale, 'Space Mouse, OLIGO kernel, Rescaled'),
                      extract_coef(glm_ecoli, 'E. coli, OLIGO kernel'),
                      extract_coef(glm_ecoli_exp, 'E. coli, OLIGO kernel, Additional Explanatory'),
                      extract_coef(glm_ecoli_rescale, 'E. coli, OLIGO kernel, Rescaled'))

coef_diagram_supplement$dataset <- c(rep('Mock Community, LC 260, OLIGO kernel', 3), 
                                     rep('Mock Community, LC 260, Linear kernel', 3),
                                     rep('Mock Community, LC 260, RBF kernel', 3),
                                     rep('Mock Community, LC 460, OLIGO kernel', 3), 
                                     rep('Mock Community, LC 460, Linear kernel', 3),
                                     rep('Mock Community, LC 460, RBF kernel', 3),
                                     rep('Diseased Oak Metagenome', 3),
                                     rep('Diseased Oak Metatranscriptome', 3),
                                     # rep('Diseased Oak Metagenome, \nOak tree genome appended', 3),
                                     rep('Ant Fungus Garden, OLIGO kernel', 3),
                                     rep('Ant Fungus Garden, Linear kernel', 3),
                                     rep('Ant Fungus Garden, RBF kernel', 3),
                                     rep('Prostate Cancer Biomarkers, OLIGO kernel', 3),
                                     rep('Space Mouse, OLIGO kernel', 3),
                                     rep('E. coli, OLIGO kernel', 3))
coef_diagram_supplement$model_structure <- rep(c('Cofrag. Score Alone', 'Additional Explanatory Variables', 'Rescaled Cofrag. Score'), 14)

coef_plot_suppl <- ggplot(coef_diagram_supplement, aes(x = dataset, y = coef)) + 
  coord_flip() +
  geom_point(position = position_dodge(0.2), aes(x = dataset, shape = model_structure, colour = model_structure), alpha = 0.4) +
  geom_linerange(aes(ymin = coef - se, ymax = coef + se, x = dataset, colour = model_structure), alpha = 0.4, 
                 position = position_dodge(0.2)) +
  theme_bw() +
  theme(legend.position = 'top', legend.title = element_blank()) +
  scale_colour_manual(values = c('blue', 'red', 'orange')) +
  scale_x_discrete(limits = rev(levels(coef_diagram_supplement$dataset))) +
  ylab("Coefficient Estimate") +
  xlab("");coef_plot_suppl

dev.off()

jpeg("figures/supp_coefficient_plot_all.jpeg", width=170, height=150, units="mm", res=850)

print(coef_plot_suppl)

dev.off()

###

#  supplementary table ----------------------------------------------------

coef_diagram_supplement_table <- rbind(extract_coef2(glm_kleiner260_oligo, 'Mock Community, LC 260, OLIGO kernel'),
                                 extract_coef2(glm_kleiner260_oligo_exp, 'Mock Community, LC 260, Additional Explanatory, OLIGO kernel'),
                                 extract_coef2(glm_kleiner260_oligo_rescale, 'Mock Community, LC 260, Rescaled, OLIGO kernel'),
                                 extract_coef2(glm_kleiner260_linear, 'Mock Community, LC 260, Linear kernel'),
                                 extract_coef2(glm_kleiner260_linear_exp, 'Mock Community, LC 260, Additional Explanatory, Linear kernel'),
                                 extract_coef2(glm_kleiner260_linear_rescale, 'Mock Community, LC 260, Rescaled, linear kernel'),
                                 extract_coef2(glm_kleiner260_rbf, 'Mock Community, LC 260, RBF kernel'),
                                 extract_coef2(glm_kleiner260_rbf_exp, 'Mock Community, LC 260, Additional Explanatory, RBF kernel'),
                                 extract_coef2(glm_kleiner260_rbf_rescale, 'Mock Community, LC 260, Rescaled, RBF kernel'),
                                 extract_coef2(glm_kleiner460_oligo, 'Mock Community, LC 460, OLIGO kernel'),
                                 extract_coef2(glm_kleiner460_oligo_exp, 'Mock Community, LC 460, Additional Explanatory, OLIGO kernel'),
                                 extract_coef2(glm_kleiner460_oligo_rescale, 'Mock Community, LC 460, Rescaled, OLIGO kernel'),
                                 extract_coef2(glm_kleiner460_linear, 'Mock Community, LC 460, Linear kernel'),
                                 extract_coef2(glm_kleiner460_linear_exp, 'Mock Community, LC 460, Additional Explanatory, Linear kernel'),
                                 extract_coef2(glm_kleiner460_linear_rescale, 'Mock Community, LC 460, Rescaled, linear kernel'),
                                 extract_coef2(glm_kleiner460_rbf, 'Mock Community, LC 460, RBF kernel'),
                                 extract_coef2(glm_kleiner460_rbf_exp, 'Mock Community, LC 460, Additional Explanatory, RBF kernel'),
                                 extract_coef2(glm_kleiner460_rbf_rescale, 'Mock Community, LC 460, Rescaled, RBF kernel'),
                                 extract_coef2(glm_broberg_mg, 'Diseased Oak Tree Metagenome'),
                                 extract_coef2(glm_broberg_mg_exp, 'Diseased Oak Tree Metagenome, Additional Explanatory'),
                                 extract_coef2(glm_broberg_mg_rescale, 'Diseased Oak Tree Metagenome, Rescaled'),
                                 extract_coef2(glm_broberg_mt, 'Diseased Oak Tree Metatranscriptome'),
                                 extract_coef2(glm_broberg_mt_exp, 'Diseased Oak Tree Metatranscriptome, Additional Explanatory'),
                                 extract_coef2(glm_broberg_mt_rescale, 'Diseased Oak Tree Metatranscriptome, Rescaled'),
                                 # extract_coef2(glm_broberg_mg_oak, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                                 # extract_coef2(glm_broberg_mg_oak_exp, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended, Additional Explanatory'),
                                 # extract_coef2(glm_broberg_mg_oak_rescale, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended, Rescaled'),
                                 extract_coef2(glm_aylward_oligo, 'Ant Fungus Garden, OLIGO kernel'),
                                 extract_coef2(glm_aylward_oligo_exp, 'Ant Fungus Garden, OLIGO kernel, Additional Explanatory'),
                                 extract_coef2(glm_aylward_oligo_rescale, 'Ant Fungus Garden, OLIGO Kernel, Rescaled'),
                                 extract_coef2(glm_aylward_linear, 'Ant Fungus Garden, Linear kernel'),
                                 extract_coef2(glm_aylward_linear_exp, 'Ant Fungus Garden, Linear kernel, Additional Explanatory'),
                                 extract_coef2(glm_aylward_linear_rescale, 'Ant Fungus Garden, Linear Kernel, Rescaled'),
                                 extract_coef2(glm_aylward_rbf, 'Ant Fungus Garden, RBF kernel'),
                                 extract_coef2(glm_aylward_rbf_exp, 'Ant Fungus Garden, RBF kernel, Additional Explanatory'),
                                 extract_coef2(glm_aylward_rbf_rescale, 'Ant Fungus Garden, RBF Kernel, Rescaled'),
                                 extract_coef2(glm_prostate_cancer, 'Prostate Cancer Biomarkers, OLIGO kernel'),
                                 extract_coef2(glm_prostate_cancer_exp, 'Prostate Cancer Biomarkers, OLIGO kernel, Additional Explanatory'),
                                 extract_coef2(glm_prostate_cancer_rescale, 'Prostate Cancer Biomarkers, OLIGO kernel, Rescaled'),
                                 extract_coef2(glm_space_mouse, 'Space Mouse, OLIGO kernel'),
                                 extract_coef2(glm_space_mouse_exp, 'Space Mouse, OLIGO kernel, Additional Explanatory'),
                                 extract_coef2(glm_space_mouse_rescale, 'Space Mouse, OLIGO kernel, Rescaled'),
                                 extract_coef2(glm_ecoli, 'E. coli, OLIGO kernel'),
                                 extract_coef2(glm_ecoli_exp, 'E. coli, OLIGO kernel, Additional Explanatory'),
                                 extract_coef2(glm_ecoli_rescale, 'E. coli, OLIGO kernel, Rescaled'))

write.csv(coef_diagram_supplement_table, "tables/supplementary-table-1_all-glms.csv")

