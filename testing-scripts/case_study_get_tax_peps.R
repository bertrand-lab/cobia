# case study: taxonomic and functional biases

# this scripts goes into taxonomic groups from the metatranscriptome to determine if a peptide is taxon specific

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

annot_contigs <- read_excel("../data/bertrand_data/antarctica_2013_MCM_FeVit_annotations.xlsx", skip = 1)
annot_contigs <- annot_contigs %>% dplyr::rename(contig = orf_id)

bertrand_cofrag <- read.csv("../data/bertrand_data/orfs.filtered.pep.trypsin_global_mi-0.00833333_ipw-0.725_para-30_co-sim.csv",
                            stringsAsFactors = FALSE)
bertrand_cofrag2 <- bertrand_cofrag[complete.cases(bertrand_cofrag), ]

get_group_specific_peptides <- function(cofrag_df, subset_grouping){
  
  # tester <- bertrand_cofrag2 function to filter out peptides that are not group specific subset_grouping <- "group"
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
      print("group")
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

tax_group_peps <- get_group_specific_peptides(cofrag_df = bertrand_cofrag2, subset_grouping = "group")
tax_group_peps_df <- data.frame(peptide_index_tax_group = tax_group_peps)
write.csv(tax_group_peps_df, "../data/bertrand_data/tax_group_pep_indexes.csv")

