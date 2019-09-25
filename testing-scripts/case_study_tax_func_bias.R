# case study: taxonomic and functional biases

library(ggplot2)
library(dplyr)
library(readxl)
library(seqinr)
library(reshape2)
library(cleaver)
library(Peptides)
library(Biostrings)
library(entropy)
library(gridExtra)

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2)) 
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# read in file with ORF ids and annotations

# annot_contigs <- read_excel("data/bertrand_tfg_data/annotation_allTFG.mmetsp_fc_pn_reclassified.edgeR.xlsx")

annot_contigs <- read_excel("data/bertrand_data/antarctica_2013_MCM_FeVit_annotations.xlsx", 
                            skip = 1)

bertrand_cofrag <- read.csv("data/bertrand_data/orfs.filtered.pep.trypsin_global_mi-0.00833333_ipw-0.725_para-30_co-sim.csv",
                            stringsAsFactors = FALSE)
bertrand_cofrag2 <- bertrand_cofrag[complete.cases(bertrand_cofrag), ]

get_group_specific_peptides <- function(cofrag_df, subset_grouping){
  
  # tester <- bertrand_cofrag2
  # function to filter out peptides that are not group specific
  # subset_grouping <- "group"
  recorded_peptides <- vector()
  
  tally_peps <- as.data.frame(table(cofrag_df$peptide_sequence))
  unique_peps_tally <- tally_peps[tally_peps$Freq == 1, ]$Var1 %>% as.character()
  
  # make a list of unique peptides
  unique_peptides_obs <- tally_peps[tally_peps$Freq > 1, ]$Var1 %>% as.character()

  # loop through each unique peptide
  for(i in 1:length(unique_peptides_obs)){
    # get matching contigs for each unique peptide
    print(i)
    # print(unique_peptides_obs)
    contig_for_peptide <- cofrag_df[cofrag_df$peptide_sequence == unique_peptides_obs[i], ]$contig
    if(length(contig_for_peptide) == 1){
      recorded_peptides <- c(recorded_peptides, unique_peptides_obs[i])
      next
    }
    if(subset_grouping == "group"){
      subsample_annot_file <- annot_contigs[annot_contigs$contig %in% contig_for_peptide, ]$group %>% unique()
    }
    if(subset_grouping == "KOG_class"){
      subsample_annot_file <- annot_contigs[annot_contigs$contig %in% contig_for_peptide, ]$KOG_class %>% unique()
    }
    if(length(subsample_annot_file) == 1){
      recorded_peptides <- c(recorded_peptides, unique_peptides_obs[i])
    } else {
      next
    }
  }
  return(c(recorded_peptides, unique_peps_tally))
}

# blah <- get_group_specific_peptides(cofrag_df = bertrand_cofrag2[c(1:200),], subset_grouping = "group")

# tax_group_peps <- get_group_specific_peptides(cofrag_df = bertrand_cofrag2, subset_grouping = "group")
# tax_group_peps_df <- data.frame(peptide_index_tax_group = tax_group_peps)
# write.csv(tax_group_peps_df, "data/bertrand_data/tax_group_pep_indexes.csv")

tax_group_peps <- read.csv("data/bertrand_data/tax_group_pep_indexes.csv")
func_group_peps <- read.csv("data/bertrand_data/func_group_pep_indexes.csv")

dim(tax_group_peps)
dim(func_group_peps)

bertrand_cofrag_tax_specific <- bertrand_cofrag2[bertrand_cofrag2$peptide_sequence %in% tax_group_peps$peptide_index_tax_group, ]
bertrand_cofrag_func_specific <- bertrand_cofrag2[bertrand_cofrag2$peptide_sequence %in% func_group_peps$peptide_index_func_group, ]

# func_group_peps <- get_group_specific_peptides(cofrag_df = bertrand_cofrag2, subset_grouping = "KOG_class")

# filter the annotation file based on these matched contigs

# determine the number of unique groups within this filtered annotation file

# if greater == 1, add peptide to a list

# if not, continue



# merge annotations with peptide specific scores
# head(annot_contigs)
# head(bertrand_cofrag)

# consider only peptides that are contig-specific


# map those peptides to contigs to get the average, median, and range of cofragmentation scores
annot_contigs <- annot_contigs %>% dplyr::rename(contig = orf_id)

# summarize the bertrand_cofrag file by contig

format_contig_file <- function(dataframe_input){
  
  bertrand_cofrag_contig <- dataframe_input %>% 
    group_by(contig) %>% 
    summarise(contig_mean = mean(mean_cofrag_score, na.rm = TRUE),
              contig_median = median(mean_cofrag_score, na.rm = TRUE),
              contig_min = min(mean_cofrag_score, na.rm = TRUE),
              contig_max = max(mean_cofrag_score, na.rm = TRUE),
              number_peps_per_contig = n())
  
  return(bertrand_cofrag_contig)
}

bertrand_cofrag_tax_specific_form <- format_contig_file(dataframe_input = bertrand_cofrag_tax_specific)
bertrand_cofrag_func_specific_form <- format_contig_file(dataframe_input = bertrand_cofrag_func_specific)

bertrand_cofrag_annot_tax <- inner_join(bertrand_cofrag_tax_specific_form, annot_contigs, by = 'contig')
bertrand_cofrag_annot_func <- inner_join(bertrand_cofrag_func_specific_form, annot_contigs, by = 'contig')

bertrand_cofrag_annot_no_group <- bertrand_cofrag_annot_tax[,-27]
bertrand_cofrag_annot_no_kog_class <- bertrand_cofrag_annot_func[,-19]

# note that there are some contigs that are in the cofragmentation output but NOT in the annotation output.
# this is because there was an abundance threshold for whether or not to annotate a transcript from the 
# metatranscriptomic output, but this threshold was not put for whether or not to put in the fasta
# sequence file.

# convert NA to string

bertrand_cofrag_annot_tax[is.na(bertrand_cofrag_annot_tax$group), ]$group <- rep("NA", length(bertrand_cofrag_annot_tax[is.na(bertrand_cofrag_annot_tax$group), ]$group))
bertrand_cofrag_annot_func[is.na(bertrand_cofrag_annot_func$KOG_class), ]$KOG_class <- rep("NA", length(bertrand_cofrag_annot_func[is.na(bertrand_cofrag_annot_func$KOG_class), ]$KOG_class))

# calculate the KL divergence between all groupings together and each subset

calc_tax_kl_div <- function(target_tax, subset_grouping = "group", annot_df){
  
  # target_tax <- 'Rhodophyta'
  # subsetting dataframe for the target taxonomy and removing NAs
  if(subset_grouping == "group"){
    subset_1 <- annot_df[annot_df$group == target_tax, ]$contig_min[!is.na(annot_df[annot_df$group == target_tax, ]$contig_min)]
  }
  if(subset_grouping == "KOG_class"){
    subset_1 <- annot_df[annot_df$KOG_class == target_tax, ]$contig_min[!is.na(annot_df[annot_df$KOG_class == target_tax, ]$contig_min)]
    
  }
  # adding the maximum which slightly skews the results by a very small amount but makes the bins identical
  subset_1_app <- c(subset_1, max(annot_df$contig_min))
  
  # the total distribution
  subset_2 <- annot_df$contig_min[!is.na(annot_df$contig_min)] 
  
  # discretized distributions of subsets
  dis_1 <- entropy::discretize(x = subset_1_app, numBins = 100)
  dis_1_mod <- dis_1 + 0.000001 # small values are added to avoid a KL divergence that fails
  dis_2 <- entropy::discretize(x = subset_2, numBins = 100)
  dis_2_mod <- dis_2 + 0.000001
  
  kl_div <- entropy::KL.plugin(freqs1 = dis_2_mod, freqs2 = dis_1_mod, unit = 'log2')
  
  return(kl_div)
}

calc_kl_div_bootstrap <- function(target_tax, subset_grouping = "group", number_boot, annot_df){

    # target_tax <- 'Rhodophyta'
    # subsetting dataframe for the target taxonomy and removing NAs
    if(subset_grouping == "group"){
      subset_1 <- annot_df[annot_df$group == target_tax, ]$contig_min[!is.na(annot_df[annot_df$group == target_tax, ]$contig_min)] %>% length()
    }
    if(subset_grouping == "KOG_class"){
      subset_1 <- annot_df[annot_df$KOG_class == target_tax, ]$contig_min[!is.na(annot_df[annot_df$KOG_class == target_tax, ]$contig_min)] %>% length()
      
    }
    print(subset_1)  
      # the total distribution
    subset_2 <- annot_df$contig_min[!is.na(annot_df$contig_min)] 
    
    kl_div_vector <- vector(length = number_boot)
    
    for(i in 1:number_boot){
      random_sub_df <- sample_n(tbl = annot_df, size = subset_1)
      print(nrow(random_sub_df))
      
      random_sub_contig_min <- random_sub_df$contig_min
      # adding the maximum which slightly skews the results by a very small amount but makes the bins identical
      subset_1_app <- c(random_sub_contig_min, max(annot_df$contig_min))
      
      # discretized distributions of subsets
      dis_1 <- entropy::discretize(x = subset_1_app, numBins = 100)
      dis_1_mod <- dis_1 + 0.000001 # small values are added to avoid a KL divergence that fails
      dis_2 <- entropy::discretize(x = subset_2, numBins = 100)
      dis_2_mod <- dis_2 + 0.000001

      kl_div <- entropy::KL.plugin(freqs1 = dis_2_mod, freqs2 = dis_1_mod, unit = 'log2')
      print(kl_div)
      kl_div_vector[i] <- kl_div
    }
    
    quantile95 <- quantile(kl_div_vector, 0.95) %>% as.numeric()
    quantile05 <- quantile(kl_div_vector, 0.05) %>% as.numeric()
    
    if(subset_grouping == "group"){
      df1 <- data.frame(group = target_tax, 
                        KL_divergence95 = quantile95, 
                        KL_divergence05 = quantile05)
    }
    if(subset_grouping == "KOG_class"){
      df1 <- data.frame(KOG_class = target_tax, 
                        KL_divergence95 = quantile95, 
                        KL_divergence05 = quantile05)
    }
    
    return(df1)
}

# diatom_test <- calc_kl_div_bootstrap(target_tax = 'Centric Diatom', subset_grouping = 'group', number_boot = 1)
# diatom_test_noboot <- calc_tax_kl_div(target_tax = 'Centric Diatom', subset_grouping = 'group')

calc_multiple_kl_div_boot <- function(list_of_tax, subset_grouping = "group", number_boot, annot_df){
  
  if(subset_grouping == "group"){
    master_df <- data.frame(group = character(), 
                            KL_divergence95 = numeric(), 
                            KL_divergence05 = numeric())
  }
  if(subset_grouping == "KOG_class"){
    master_df <- data.frame(KOG_class = character(), 
                      KL_divergence95 = numeric(), 
                      KL_divergence05 = numeric())
  }  
  
  for(individual_tax in 1:length(list_of_tax)){
    temp_tax <- calc_kl_div_bootstrap(target_tax = list_of_tax[individual_tax], 
                                      subset_grouping = subset_grouping, 
                                      number_boot = number_boot, annot_df = annot_df)
    master_df <- rbind(temp_tax, master_df)
    # print(temp_tax)
    # print(master_df)
  }
  
  return(master_df)

}

calc_multiple_kl_div <- function(list_of_tax, subset_grouping = "group", number_boot, annot_df){
  
  list_of_kl <- vector(length = length(list_of_tax))
  
  for(individual_tax in 1:length(list_of_tax)){
    temp_tax <- calc_tax_kl_div(target_tax = list_of_tax[individual_tax], 
                                subset_grouping = subset_grouping,
                                annot_df = annot_df)
    list_of_kl[individual_tax] <- temp_tax
  }
  
  return(list_of_kl)
  
}

# kl diverences for all groups
kl_div_group <- calc_multiple_kl_div(list_of_tax = unique(bertrand_cofrag_annot_tax$group)[complete.cases(unique(bertrand_cofrag_annot_tax$group))],
                              subset_grouping = "group", annot_df = bertrand_cofrag_annot_tax)
kl_div_kog_class <- calc_multiple_kl_div(list_of_tax = unique(bertrand_cofrag_annot_func$KOG_class)[complete.cases(unique(bertrand_cofrag_annot_func$KOG_class))],
                              subset_grouping = "KOG_class", annot_df = bertrand_cofrag_annot_func)

set.seed(10101)
# bootstrapped confidence intervals
kl_div_group_boot <- calc_multiple_kl_div_boot(list_of_tax = unique(bertrand_cofrag_annot_tax$group)[complete.cases(unique(bertrand_cofrag_annot_tax$group))],
                                                   subset_grouping = "group", number_boot = 1000, 
                                               annot_df = bertrand_cofrag_annot_tax)
kl_div_kog_class_boot <- calc_multiple_kl_div_boot(list_of_tax = unique(bertrand_cofrag_annot_func$KOG_class)[complete.cases(unique(bertrand_cofrag_annot_func$KOG_class))],
                                              subset_grouping = "KOG_class", number_boot = 1000, 
                                              annot_df = bertrand_cofrag_annot_func)

df1_group_kl_div <- data.frame(KL_divergence = round(kl_div_group, 3), 
           group = unique(bertrand_cofrag_annot_tax$group)[complete.cases(unique(bertrand_cofrag_annot_tax$group))])

df1_kog_class_kl_div <- data.frame(KL_divergence = round(kl_div_kog_class, 3), 
           KOG_class = unique(bertrand_cofrag_annot_func$KOG_class)[complete.cases(unique(bertrand_cofrag_annot_func$KOG_class))])



cofrag_tax_plot_supp <- bertrand_cofrag_annot_tax %>% 
  ggplot(aes(x = contig_min)) + 
  geom_density(data = bertrand_cofrag_annot_no_group, 
               fill = 'grey90', lwd = 0.1) +
  geom_density(colour = 'grey10', alpha = 0.8) + 
  facet_wrap(~group, labeller = label_wrap_gen(width=10)) +
  geom_label(data = df1_group_kl_div, aes(label = KL_divergence), 
             x = 45, y = 0.065, 
             hjust=1, vjust=0,
             inherit.aes = FALSE) +
  theme_bw() +
  xlim(0, 50) +
  xlab("Minimum Cofragmentation Score by Open Reading Frame") +
  ylab("Probability Density");cofrag_tax_plot_supp

ggsave(cofrag_tax_plot_supp, 
       filename = "figures/supp_cofrag_tax_plot.jpeg")


tally_groups <- bertrand_cofrag_annot_tax %>% 
  group_by(group) %>% 
  dplyr::tally() %>% as.data.frame()

cofrag_func_plot_supp <- bertrand_cofrag_annot_func %>% 
  ggplot(aes(x = contig_min)) + 
  geom_density(data = bertrand_cofrag_annot_no_kog_class, 
               fill = 'grey90', lwd = 0.1) +
  geom_density(colour = 'grey10', alpha = 0.8) + 
  facet_wrap(~KOG_class, labeller = label_wrap_gen(width=30)) +
  geom_label(data = df1_kog_class_kl_div, aes(label = KL_divergence),
             x = 45, y = 0.05, 
             hjust=1, vjust=0,
             inherit.aes = FALSE) +
  theme_bw() +
  xlim(0, 50) +
  xlab("Minimum Cofragmentation Score by Open Reading Frame") +
  ylab("Probability Density")

ggsave(cofrag_func_plot_supp, 
       filename = "figures/supp_cofrag_func_plot.jpeg")

tally_kogs <- bertrand_cofrag_annot_func %>% 
  group_by(KOG_class) %>% 
  dplyr::tally() %>% as.data.frame()
  
ggplot(aes(y = n, x = KOG_class)) + 
  geom_bar(stat='identity') +
  facet_wrap(~KOG_class) +
  theme(axis.text.x = element_blank()) +
  theme_bw()


kog_class_together <- inner_join(tally_kogs, df1_kog_class_kl_div)
kog_class_together2 <- inner_join(kog_class_together, kl_div_kog_class_boot)
group_together <- inner_join(tally_groups, df1_group_kl_div)
group_together_2 <- inner_join(group_together, kl_div_group_boot)

kog_class_plot <- kog_class_together2 %>% 
  ggplot(aes(y = KL_divergence, x = log(n))) + 
  geom_ribbon(aes(ymin = KL_divergence05, ymax = KL_divergence95), alpha = 0.1) +
  geom_point() +
  theme_bw() +
  ylab("Kullback-Leibler Divergence") +
  xlab("log(Number of assigned ORFs per KOG class)");kog_class_plot

group_plot <- group_together_2 %>% 
  ggplot(aes(y = KL_divergence, x = log(n))) + 
  geom_ribbon(aes(ymin = KL_divergence05, ymax = KL_divergence95), alpha = 0.1) +
  geom_point() +
  theme_bw() +
  ylab("Kullback-Leibler Divergence") +
  xlab("log(Number of assigned ORFs per taxonomic Group)");group_plot

grid.arrange(kog_class_plot, group_plot, nrow = 1)


kog_class_together2$Grouping <- rep("KOG Class", nrow(kog_class_together2))
group_together_2$Grouping <- rep("Taxonomic Group", nrow(group_together_2))
kog_class_together3 <- dplyr::rename(.data = kog_class_together2, 
                                     subgrouping = KOG_class)

group_together3 <- dplyr::rename(.data = group_together_2,
                                 subgrouping = group)
kl_plot_df <- rbind(kog_class_together3, group_together3)

kl_plot_df %>% 
  ggplot(aes(y = KL_divergence, x = log(n))) + 
  geom_ribbon(aes(ymin = KL_divergence05, ymax = KL_divergence95), alpha = 0.1) +
  geom_point() +
  theme_bw() +
  facet_grid(~Grouping) +
  ylab("Kullback-Leibler Divergence") +
  xlab("log(Number of assigned ORFs per Grouping (taxonomic or functional))")
  
