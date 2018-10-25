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
library(seqinr)

# data processing scripts

# formatting cofragmentation output for validation
cofrag_processor <- function(cofrag_output_df, observed_peps){
  
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

# Kleiner Run 1 Dataset ---------------------------------------------------

# openms idXML file filtered

run1_oms <- read.csv("data/kleiner_data/Run1_C4_832ng_IDF.csv")


# reading in different cofrag scores from different lc predictions
kleiner_cofrag_460_rbf <- read.csv("data/kleiner_data/Run1_C4_832ng_rbf_cofrag_mi-0.00833333_ipw-2.75_para-15_co-sim.csv")
kleiner_cofrag_460_oligo <-  read.csv("data/kleiner_data/Run1_C4_832ng_oligo_cofrag_mi-0.00833333_ipw-2.75_para-15_co-sim.csv")

# removing the modX modifications to N termini and C termini
kleiner_cofrag_460_rbf2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_460_rbf, observed_peps = run1_oms$sequence)
kleiner_cofrag_460_oligo2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_460_oligo, observed_peps = run1_oms$sequence)

number_of_unseen(kleiner_cofrag_460_oligo2)
kleiner460_n <- "italic(N)[italic(unobserved)] == 2164345"
kleiner460_n_sub <- "italic(N)[italic(observed)] == 5920"
kleiner_460_p2 <- kleiner_cofrag_460_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.12, "(0.009)"), y = 1e8, x = 13, parse = TRUE, size = 7) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  annotate("text", x = 6, y = Inf, vjust = 3, label = kleiner460_n, parse = TRUE, size = 3, colour = 'grey39') +
  annotate("text", x = 6, y = Inf, vjust = 4.5, label = kleiner460_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  ggtitle('B) Mock Community (LC Run = 460 mins)') +
  scale_fill_manual(values = c("grey39", "grey24")) +
  theme(legend.position = "none", title = element_text(size = 7.5));kleiner_460_p2

glm_kleiner460_rbf <- glm(data = kleiner_cofrag_460_rbf2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_kleiner460_rbf_rescale <- glm(data = kleiner_cofrag_460_rbf2, prote_d ~ rescale_cofrag_score, family = "binomial")
glm_kleiner460_oligo <- glm(data = kleiner_cofrag_460_oligo2, prote_d ~ mean_cofrag_score, family = "binomial")
glm_kleiner460_oligo_rescale <- glm(data = kleiner_cofrag_460_oligo2, prote_d ~ rescale_cofrag_score, family = "binomial")

glm_kleiner460_rbf_summary <- summary(glm_kleiner460_rbf)
glm_kleiner460_oligo_summary <- summary(glm_kleiner460_oligo)


# kleiner run 4 -----------------------------------------------------------

# peptides detected using OpenMS pipeline
run4_oms <- read.csv("data/kleiner_data/Run4_C4_2000ng_IDF.csv")

kleiner_cofrag_260_rbf <- read.csv("data/kleiner_data/Run4_C4_2000ng_rbf_cofrag_mi-0.00833333_ipw-1.44_para-15_co-sim.csv")
kleiner_cofrag_260_oligo <- read.csv("data/kleiner_data/Run4_C4_2000ng_oligo_cofrag_mi-0.00833333_ipw-1.44_para-15_co-sim.csv")

kleiner_cofrag_260_rbf2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_260_rbf, observed_peps = run4_oms$sequence)
kleiner_cofrag_260_oligo2 <- cofrag_processor(cofrag_output_df = kleiner_cofrag_260_oligo, observed_peps = run4_oms$sequence)

# modelling whether cofragmentation risk has an influence on peptide presence/absence
glm_kleiner260_rbf <- glm(data = kleiner_cofrag_260_rbf2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_kleiner260_rbf_rescale <- glm(data = kleiner_cofrag_260_rbf2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_kleiner260_oligo <- glm(data = kleiner_cofrag_260_oligo2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_kleiner260_oligo_rescale <- glm(data = kleiner_cofrag_260_oligo2, prote_d ~ rescale_cofrag_score, family = 'binomial')

glm_kleiner260_rbf_summary <- summary(glm_kleiner260_rbf)
glm_kleiner260_oligo_summary <- summary(glm_kleiner260_oligo)

number_of_unseen(kleiner_cofrag_260_oligo2)
kleiner260_n <- "italic(N)[italic(unobserved)] == 2153169"
kleiner260_n_sub <- "italic(N)[italic(observed)] == 11096"
kleiner_260_p1 <- kleiner_cofrag_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 0.2) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('log(Count)') +
  # annotate("text", label = paste("beta == ", -0.071, "(0.013)"), 
  # y = 1e9, x = 7, parse = TRUE, size = 6) +
  # xlab('Predicted # Cofragmenting Peptides') +
  annotate("text", x = 7, y = Inf, vjust = 3, label = kleiner260_n, parse = TRUE, size = 3, colour = 'grey39') +
  annotate("text", x = 7, y = Inf, vjust = 4.5, label = kleiner260_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  xlab('') +
  scale_fill_manual(values = c("grey39", "grey24")) +
  ggtitle('A) Mock Community (LC Run = 260 mins)') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));kleiner_260_p1


glm_kleiner260_rbf_exp <- glm(data = kleiner_cofrag_rbf, prote_d ~ mean_cofrag_score + rts + peplen, family = 'binomial')
glm_kleiner260_oligo_exp <- glm(data = kleiner_cofrag_oligo, prote_d ~ mean_cofrag_score + rts + peplen, family = 'binomial')

glm_kleiner260_rbf_exp_summary <- summary(glm_kleiner260_rbf_exp)
glm_kleiner260_oligo_exp_summary <- summary(glm_kleiner260_oligo_exp) 



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
broberg_mg_n <- "italic(N)[italic(unobserved)] == 106452"
broberg_mg_n_sub <- "italic(N)[italic(observed)] == 178"
broberg_mg_p3 <- bro_co_mg2 %>% 
  # filter(peplen > 8) %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() +
  theme_bw() +
  ylab('log(Count)') +
  scale_fill_manual(values = c("grey39", "grey24")) +
  ggtitle('C) Diseased Oak Tree Metagenome') +
  annotate("text", x = 80, y = Inf, vjust = 3, label = broberg_mg_n, parse = TRUE, size = 3, colour = 'grey39') +
  annotate("text", x = 80, y = Inf, vjust = 4.5, label = broberg_mg_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  # annotate("text", label = paste("beta == ", -0.054, "(0.01)"), 
  # y = 1e6, x = 35, parse = TRUE, size = 6) +
  # xlab('Predicted # Cofragmenting Peptides') +
  xlab('') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));broberg_mg_p3

number_of_unseen(bro_co_mt2)
broberg_mt_n <- "italic(N)[italic(unobserved)] == 337301"
broberg_mt_n_sub <- "italic(N)[italic(observed)] == 473"
broberg_mt_p4 <- bro_co_mt2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  scale_fill_manual(values = c("grey39", "grey24")) +
  ggtitle('D) Diseased Oak Tree Metatranscriptome') +
  annotate("text", x = 250, y = Inf, vjust = 3, label = broberg_mt_n, parse = TRUE, size = 3, colour = 'grey39') +
  annotate("text", x = 250, y = Inf, vjust = 4.5, label = broberg_mt_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  # annotate("text", label = paste("beta == ", -0.021, "(0.0023)"), 
  # y = 1e7, x = 80, parse = TRUE, size = 6) +
  xlab('') + theme(legend.position = "none", title = element_text(size = 7.5));broberg_mt_p4
# geom_density(alpha = 0.4)


number_of_unseen(bro_co_mg_oak2)
broberg_mg_oak_n <- "italic(N)[italic(unobserved)] == 672692"
broberg_mg_oak_n_sub <- "italic(N)[italic(observed)] == 439"
broberg_mg_p3_oak <- bro_co_mg_oak2 %>% 
  # filter(peplen > 8) %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('log(Count)') +
  scale_fill_manual(values = c("grey39", "grey24")) +
  ggtitle('E) Diseased Oak Tree Genome + Metagenome') +
  annotate("text", x = 375, y = Inf, vjust = 3, label = broberg_mg_oak_n, parse = TRUE, size = 3, colour = 'grey39') +
  annotate("text", x = 375, y = Inf, vjust = 4.5, label = broberg_mg_oak_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  xlab('Cofragmentation Score') +
  theme(legend.position = "none", 
        title = element_text(size = 7.5));broberg_mg_p3_oak


# glm_broberg_mg <- glm(data = bro_mt_joined, prote_d ~ mean_cofrag_score, family = 'binomial')

glm_broberg_mt <- glm(data = bro_co_mt2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mt_rescale <- glm(data = bro_co_mt2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_broberg_mg <- glm(data = bro_co_mg2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mg_rescale <- glm(data = bro_co_mg2, prote_d ~ rescale_cofrag_score, family = 'binomial')
glm_broberg_mg_oak <- glm(data = bro_co_mg_oak2, prote_d ~ mean_cofrag_score, family = 'binomial')
glm_broberg_mg_oak_rescale <- glm(data = bro_co_mg_oak2, prote_d ~ rescale_cofrag_score, family = 'binomial')

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

# aylward_rt_cofrag_linear2$pep_seq %!in% aylward_prot_data_orbi$pep_seq
# aylward_prot_data_orbi$pep_seq[aylward_prot_data_orbi$pep_seq %!in% aylward_rt_cofrag_linear2$pep_seq]
##

number_of_unseen(aylward_rt_cofrag_oligo2)
aylward_oligo_n <- "italic(N)[italic(unobserved)] == 2750772"
aylward_oligo_n_sub <- "italic(N)[italic(observed)] == 61"
aylward_p5_oligo <- aylward_rt_cofrag_oligo2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d), binwidth = 1) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Cofragmentation Score') + 
  annotate("text", x = 97, y = Inf, vjust = 3, label = aylward_oligo_n, parse = TRUE, colour = 'grey39', size = 3) +
  annotate("text", x = 97, y = Inf, vjust = 4.5, label = aylward_oligo_n_sub, parse = TRUE, colour = 'grey24', size = 3) +
  scale_fill_manual(values = c("grey39", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('F) Ant Fungus Garden Metagenome');aylward_p5_oligo


number_of_unseen(aylward_rt_cofrag_linear2)
aylward_oligo_n <- "italic(N)[italic(unobserved)] == 2750772"
aylward_oligo_n_sub <- "italic(N)[italic(observed)] == 61"
aylward_p5_linear <- aylward_rt_cofrag_linear2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Predicted # Cofragmenting Peptides') + 
  scale_fill_manual(values = c("grey39", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('F) Ant Fungus Garden Metagenome');aylward_p5_linear

aylward_p5_rbf <- aylward_rt_cofrag_rbf2 %>% 
  ggplot(aes(mean_cofrag_score)) + 
  geom_histogram(aes(fill = prote_d)) + 
  scale_y_log10() + 
  theme_bw() +
  ylab('') +
  # annotate("text", label = paste("beta == ", -0.043, "(0.16)"), 
  # y = 1e7, x = 7, parse = TRUE, size = 6) +
  xlab('Predicted # Cofragmenting Peptides') + 
  scale_fill_manual(values = c("grey39", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + 
  ggtitle('F) Ant Fungus Garden Metagenome');aylward_p5_rbf

glm_aylward_oligo <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_oligo2)
glm_aylward_oligo_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_oligo2)
glm_aylward_linear <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_linear2)
glm_aylward_linear_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_linear2)
glm_aylward_rbf <- glm(prote_d ~ mean_cofrag_score, family = "binomial", data = aylward_rt_cofrag_rbf2)
glm_aylward_rbf_rescale <- glm(prote_d ~ rescale_cofrag_score, family = "binomial", data = aylward_rt_cofrag_rbf2)

summary(glm_aylward_oligo)
summary(glm_aylward_linear)
summary(glm_aylward_rbf)

glm_aylward_oligo_exp <- glm(prote_d ~ mean_cofrag_score + rts + iso_point + mz + pep_len, 
                             family = "binomial", data = aylward_rt_cofrag_oligo2)
glm_aylward_linear_exp <- glm(prote_d ~ mean_cofrag_score + rts + iso_point + mz + pep_len, 
                              family = "binomial", data = aylward_rt_cofrag_linear2)
glm_aylward_rbf_exp <- glm(prote_d ~ mean_cofrag_score + rts + iso_point + mz + pep_len, 
                           family = "binomial", data = aylward_rt_cofrag_rbf2)

summary(glm_aylward_linear_exp)
summary(glm_aylward_oligo_exp)
summary(glm_aylward_rbf_exp)




# aggregated plot ---------------------------------------------------------

# plotting aggregated validation results

library(gridExtra)
# to do

# change font size of coef plot
# add in number of peps
# add ggtitle for plot F

# Making 


coef_diagram <- rbind(extract_coef(glm_kleiner260_oligo, 'Mock Community, LC 260'),
                      extract_coef(glm_kleiner460_oligo, 'Mock Community, LC 460'),
                      extract_coef(glm_broberg_mg, 'Diseased Oak Tree Metagenome'),
                      extract_coef(glm_broberg_mt, 'Diseased Oak Tree Metatranscriptome'),
                      extract_coef(glm_broberg_mg_oak, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                      extract_coef(glm_aylward_oligo, 'Ant Fungus Garden'))
coef_diagram$scaled <- rep('Cofragmentation Score', nrow(coef_diagram))

coef_diagram_rescale <- rbind(extract_coef(glm_kleiner260_oligo_rescale, 'Mock Community, LC 260'),
                      extract_coef(glm_kleiner460_oligo_rescale, 'Mock Community, LC 460'),
                      extract_coef(glm_broberg_mg_rescale, 'Diseased Oak Tree Metagenome'),
                      extract_coef(glm_broberg_mt_rescale, 'Diseased Oak Tree Metatranscriptome'),
                      extract_coef(glm_broberg_mg_oak_rescale, 'Diseased Oak Tree Metagenome, \nOak Tree genome appended'),
                      extract_coef(glm_aylward_oligo_rescale, 'Ant Fungus Garden'))
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


jpeg("validation_plots.jpeg", width=170, height=210, units="mm", res=850)

grid.arrange(kleiner_260_p1, kleiner_460_p2, broberg_mg_p3, broberg_mt_p4, broberg_mg_p3_oak, aylward_p5_oligo, nrow = 3)

dev.off()


ggsave(coef_graph, filename = "coefficient_plot.jpeg")

