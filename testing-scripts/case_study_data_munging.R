library(ggplot2)
library(dplyr)
library(readxl)
library(seqinr)
library(reshape2)
library(cleaver)
library(Peptides)
library(Biostrings)

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2)) 
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# read in file with ORF ids and annotations

# annot_contigs <- read_excel("data/bertrand_tfg_data/annotation_allTFG.mmetsp_fc_pn_reclassified.edgeR.xlsx")

annot_contigs <- read_excel("data/bertrand_data/antarctica_2013_MCM_FeVit_annotations.xlsx", skip = 1)

vit_keywords <- c('vitamin-B12 independent', 'Cobalamin-independent')

all_tryptic_peps <- read.table(file = 'data/bertrand_data/orfs.filtered.pep.trypsin_wcontigs.txt', sep = ',')

get_contigs <- function(key_word, annotation_file){
  
  if(is.character(key_word[1]) != TRUE ){
    stop
  }
  
  orf_list_dups <- vector()
  
  for(i in 1:length(key_word)){
    annot_contigs2_orfs_kegg <- annot_contigs[grepl(pattern = key_word, x = annot_contigs$kegg_desc),]$orf_id
    annot_contigs2_orfs_kog <- annot_contigs[grepl(pattern = key_word, x = annot_contigs$KOG_desc),]$orf_id
    annot_contigs2_orfs_ko <- annot_contigs[grepl(pattern = key_word, x = annot_contigs$KO_desc),]$orf_id
    annot_contigs2_orfs_all <- annot_contigs[grepl(pattern = key_word, x = annot_contigs$best_hit_annotation),]$orf_id
    annot_contigs2_orfs_pfams <- annot_contigs[grepl(pattern = key_word_pfams, x = annot_contigs$PFams_desc),]$orf_id
    all_orfs <- c(annot_contigs2_orfs_kegg, annot_contigs2_orfs_kog, annot_contigs2_orfs_ko, annot_contigs2_orfs_all, annot_contigs2_orfs_pfams) %>% unique()
    orf_list_dups <- c(orf_list_dups, all_orfs)
  }
  all_orfs <- unique(orf_list_dups)
  
  return(all_orfs)
}

get_peps <- function(contig_list, tryptic_peptide_file){
  # of all the tryptic peptides, which peptides are in the sub group of tryptic peptides
  target_tryp_peps <- tryptic_peptide_file[which(tryptic_peptide_file$V2 %in% contig_list), ]
  nontarget_tryp_peps <- tryptic_peptide_file[which(tryptic_peptide_file$V2 %!in% contig_list), ]
  
  proteotypic_informative_tryp_peps <- target_tryp_peps[target_tryp_peps$V1 %!in% nontarget_tryp_peps$V1,]$V1 %>% as.character()
  proteotypic_informative_tryp_pep_df <- target_tryp_peps[target_tryp_peps$V1 %!in% nontarget_tryp_peps$V1,]
  
  return_list <- list(proteotypic_informative_tryp_peps, proteotypic_informative_tryp_pep_df)
  
  return(return_list)
}

write_targeted_cobia <- function(get_peps_out){
  # write file for targeted cobia
  pep_targeted <- get_peps_out
  pep_targeted$mz_nomod <- mw(pep_targeted$pep_seq)/2
  pep_targeted$len_prot <- nchar(pep_targeted$pep_seq %>% as.character())
  pep_targeted$num_m <- str_count(pep_targeted$pep_seq %>% as.character(), pattern = "M")
  pep_targeted$num_c <- str_count(pep_targeted$pep_seq %>% as.character(), pattern = "C")
  pep_targeted$mz <- pep_targeted$mz_nomod + pep_targeted$num_c*28.5 + pep_targeted$num_m*8
  pep_targeted2 <- pep_targeted %>% filter(mz < 2000, len_prot > 4)
  return(pep_targeted2)
}

get_tax_specific_peps <- function(contig_annot_file, proteotypic_peps, taxonomy_id){
  # names(proteotypic_informative_tryp_pep_df) <- c('pep_seq', 'orf_id')
  good_peptide_candidates <- inner_join(proteotypic_peps, contig_annot_file[, c(1:37)], by = 'orf_id')
  
  # subset good peptide candidates by taxonomy
  tax_specific_peps <- good_peptide_candidates[good_peptide_candidates$best_LPI_species == taxonomy_id, ]$pep_seq %>% as.character()
  # tax specific proteins
  not_tax_specific_peps <- good_peptide_candidates[good_peptide_candidates$best_LPI_species != taxonomy_id, ]$pep_seq %>% as.character()
  
  good_tax_peps <- tax_specific_peps[tax_specific_peps %!in% not_tax_specific_peps] %>% as.character()
  
  return(good_tax_peps)
}

find_tax_peps <- function(tryptic_peptide_file, key_word, annotation_file, target_tax){
  
  # tryptic_peptide_file <- all_tryptic_peps 
  # key_word <- vit_keywords 
  # annotation_file <- annot_contigs
  # target_tax <- "Fragilariopsis cylindrus"
  
  target_contigs <- get_contigs(key_word = key_word, annotation_file = annotation_file)
  candidate_peps <- get_peps(contig_list = target_contigs, tryptic_peptide_file = tryptic_peptide_file)
  names(candidate_peps[[2]]) <- c('pep_seq', 'orf_id')
  tax_specific_peps <- get_tax_specific_peps(contig_annot_file = annotation_file, 
                                             proteotypic_peps = candidate_peps[[2]], 
                                             taxonomy_id = target_tax)
  return(tax_specific_peps)
}

find_tax_peps(tryptic_peptide_file = all_tryptic_peps, 
              key_word = vit_keywords, 
              annotation_file = annot_contigs, 
              target_tax = "Fragilariopsis cylindrus")
# write.csv(pep_targeted2['pep_seq'], row.names = FALSE, file = "data/bertrand_data/candidate_peps_targeted.csv")

write.fasta(as.list(good_frag_peps), names = seq(from = 1, to = length(good_frag_peps)), file.out = 'data/bertrand_data/good_frag_peps.fasta')

write.csv(good_frag_peps['pep_seq'], row.names = FALSE, file = "data/bertrand_data/frag_cyl_peps_metE.csv")
write.csv(atth_peps_candidates['pep_seq'], row.names = FALSE, file = "data/bertrand_data/atth_kerg_peps_metE.csv")

# # align and check the contigs
# frag_contigs <- good_peptide_candidates[good_peptide_candidates$best_LPI_species == 'Fragilariopsis cylindrus', ]$orf_id %>% unique()
# frag_fasta <- contig_seqs[which(names(contig_seqs) %in% frag_contigs)]
# write.fasta(frag_fasta, names = names(frag_fasta), file.out = 'data/bertrand_data/frag_mete.fasta')
# 
# # align and check the contigs
# frag_contigs <- good_peptide_candidates[good_peptide_candidates$best_LPI_species == 'Fragilariopsis cylindrus', ]$orf_id %>% unique()
# frag_fasta <- contig_seqs[which(names(contig_seqs) %in% frag_contigs)]
# write.fasta(frag_fasta, names = names(frag_fasta), file.out = 'data/bertrand_data/frag_mete.fasta')

# reading in cofragmentation data
test <- read.csv("data/bertrand_data/assembly.orf.fasta_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")
targ <- read.csv("data/bertrand_data/asembly.orf_frag_peps_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")

# subset the CONSEQUENCE scores of four
really_good_peps <- c("HSTFAQTEGSIDVQR", "AQAVEELGWSLQLADDK", "WFTTNYHYLPSEVDTK")

pep_lc_file <- read.csv("data/bertrand_data/assembly.orf_lc-retention-times.csv")

dda_params_file <- read.csv("data/broberg_data/dda_params_broberg.csv")

# look for peptides of similar mass and retention time as above

cofrag_buddies <- function(pep_seq_cofrag, lc_file, dda_params_file){
  # lc_file <- pep_lc_file
  # pep_seq_cofrag <- "HSTFAQTEGSIDVQR"

  rt_pep <- lc_file[lc_file$peptide_sequence == paste0(pep_seq_cofrag, '-OH'), ]$rts
  rt_upper_bound <- rt_pep + dda_params_file[1, c('ion_peak_width')]
  rt_lower_bound <- rt_pep - dda_params_file[1, c('ion_peak_width')]
  
  mz_pep <- lc_file[lc_file$peptide_sequence == paste0(pep_seq_cofrag, '-OH'), ]$mass/2
  mz_upper_pep <- mz_pep + 0.5*dda_params_file[1, c('precursor_selection_window')]
  mz_lower_pep <- mz_pep - 0.5*dda_params_file[1, c('precursor_selection_window')]
  
  lc_file$mz <- lc_file$mass/2
  
  other_peps <- lc_file %>% dplyr::filter(rts > rt_lower_bound, 
                                          rts < rt_upper_bound, 
                                          mz > mz_lower_pep, 
                                          mz < mz_upper_pep)
  
  other_peps_seqs_unique <- unique(other_peps$peptide_sequence) %>% as.character()
  other_peps_seqs <- other_peps$peptide_sequence %>% as.character()
  
  other_peps_contigs_unique <- unique(other_peps$contig) %>% as.character()
  other_peps_contigs <- other_peps$contig %>% as.character()
  
  finale_list <- list(other_peps_seqs, other_peps_contigs, other_peps_seqs_unique, other_peps_contigs_unique)
  
  return(finale_list)

  }

test <- cofrag_buddies(pep_seq_cofrag = 'HSTFAQTEGSIDVQR', lc_file = pep_lc_file, dda_params_file = dda_params_file)
cofrag_buddies(pep_seq_cofrag = 'AQAVEELGWSLQLADDK', lc_file = pep_lc_file, dda_params_file = dda_params_file)
cofrag_buddies(pep_seq_cofrag = 'WFTTNYHYLPSEVDTK', lc_file = pep_lc_file, dda_params_file = dda_params_file)

cofrag_buddy_annot <- function(cofrag_buddy_output, annot_file){
  # cofrag_buddy_output <- test 
  # annot_file <- annot_contigs

  cofrag_buddy_output_peps <- data.frame(pep_seq = cofrag_buddy_output[[1]], orf_id = cofrag_buddy_output[[2]])
  
  annot_sub <- annot_file[annot_file$orf_id %in% cofrag_buddy_output_peps[[2]], c('best_hit_annotation', 
                                                                               'kegg_desc', 
                                                                               'KOG_desc', 
                                                                               'KO_desc', 
                                                                               'best_LPI_species', 'orf_id')]
  
  annot_sub_finale <- inner_join(cofrag_buddy_output_peps, annot_sub, by = 'orf_id')
  
  # finale_df <- cbind(cofrag_buddy_output[[1]],
  #                    annot_sub$best_hit_annotation, 
  #                    annot_sub$kegg_desc,
  #                    annot_sub$KOG_desc,
  #                    annot_sub$KO_desc,
  #                    annot_sub$best_LPI_species)
  
  return(annot_sub_finale)
}

cofrag_proc <- function(pep_seq, lc_file_master = pep_lc_file, dda_params_file_master = dda_params_file, annot_file = annot_contigs){
  # pep_seq <- 'HSTFAQTEGSIDVQR'
  
  co_buddies <- cofrag_buddies(pep_seq_cofrag = pep_seq, lc_file = lc_file_master, dda_params_file = dda_params_file_master)
  co_buddies_annot <- cofrag_buddy_annot(cofrag_buddy_output = co_buddies, annot_file = annot_file)
  
  
  return(co_buddies_annot)
}



cofrag_proc(pep_seq = 'HSTFAQTEGSIDVQR')
# lots of unknown proteins. even three separately identified proteins from frag

cofrag_proc(pep_seq = 'AQAVEELGWSLQLADDK')[13,]
cofrag_proc(pep_seq = 'AQAVEELGWSLQLADDK')[7,]
cofrag_proc(pep_seq = 'AQAVEELGWSLQLADDK')[64,]

cofrag_proc(pep_seq = 'WFTTNYHYLPSEVDTK')[15,]
cofrag_proc(pep_seq = 'WFTTNYHYLPSEVDTK')[22,]
# determine which of the contigs are also in that bin

# figure out what they do biologically, and see what would happen if expression patterns changed




targ[targ$pep_seq %in% really_good_peps,]

blah <- inner_join(targ, frag_peps_cofrag, by = 'pep_seq')

plot(blah$mean_cofrag_score.x, blah$mean_cofrag_score.y)

abline(a = 0, b = 1)
plot(frag_peps_cofrag$mean_cofrag_score, targ$mean_cofrag_score)
head(test)
test$pep_seq <- str_sub(test$peptide_sequence, start = 1, end = -4)

frag_peps_cofrag <- test[test$pep_seq %in% good_frag_peps, ]

good_peptide_candidates %>% dplyr::filter(best_LPI_species == 'Karenia brevis Wilson')

# function that inputs a peptide sequence, outputs a vector of contigs

names(all_tryptic_peps) <- c('peptide_seq', 'orf_id')
names(sub_peptide_joined)[1] <- c('peptide_seq')


# for each candidate peptide

# determine which taxonomy it maps to

candidate_peps_matched <- inner_join(x = candidate_peps, y = all_tryptic_peps, by = 'peptide_seq')

dim(candidate_peps)
dim(candidate_peps_matched)

candidate_peps_tax <- left_join(x = candidate_peps_matched, y = annot_contigs3_df, by = 'orf_id')

dim(candidate_peps_matched)
dim(candidate_peps_tax)
head(candidate_peps_tax)

candidate_pennate <- candidate_peps_tax %>% 
  dplyr::filter(group == 'Pennate Diatom')
candidate_frag <- candidate_peps_tax %>% 
  dplyr::filter(best_LPI_tax_string == 'Eukaryota;Stramenopiles;Stramenopiles_X;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;Fragilariopsis;Fragilariopsis cylindrus' |
                best_LPI_tax_string == 'Eukaryota;Stramenopiles;Stramenopiles_X;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;Fragilariopsis;Fragilariopsis kerguelensis L26_C5')
candidate_pseudo <- candidate_peps_tax %>% 
  dplyr::filter(taxgrp == 'Pseudo-nitzschia')
candidate_hapto <- candidate_peps_tax %>% 
  dplyr::filter(taxgrp == 'Other Haptophyta')

candidate_not_pennate <- candidate_peps_tax %>% 
  dplyr::filter(group != 'Pennate Diatom')
candidate_not_frag <- candidate_peps_tax %>% 
  dplyr::filter(best_LPI_tax_string != 'Eukaryota;Stramenopiles;Stramenopiles_X;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;Fragilariopsis;Fragilariopsis cylindrus',
                  best_LPI_tax_string != 'Eukaryota;Stramenopiles;Stramenopiles_X;Bacillariophyta;Bacillariophyta_X;Raphid-pennate;Fragilariopsis;Fragilariopsis kerguelensis L26_C5')
candidate_not_pseudo <- candidate_peps_tax %>% 
  dplyr::filter(taxgrp != 'Pseudo-nitzschia')

candidate_not_hapto <- candidate_peps_tax %>% 
  dplyr::filter(taxgrp != 'Other Haptophyta')


unique_candidate_pennate <- candidate_pennate[candidate_pennate$peptide_seq %!in% candidate_not_pennate$peptide_seq,]$peptide_seq %>% unique()

unique_candidate_frag <- candidate_frag[candidate_frag$peptide_seq %!in% candidate_not_frag$peptide_seq,]$peptide_seq %>% unique()
unique_candidate_pseudo <- candidate_pseudo[candidate_pseudo$peptide_seq %!in% candidate_not_pseudo$peptide_seq, ]$peptide_seq %>% unique()


unique_candidate_hapto <- candidate_hapto[candidate_hapto$peptide_seq %!in% candidate_not_hapto$peptide_seq, ]$peptide_seq %>% unique()

write.csv(x = unique_candidate_frag, 'data/bertrand_tfg_data/unique_candidate_peps_frag_kog_mnsod.csv')
write.csv(x = unique_candidate_pseudo, 'data/bertrand_tfg_data/unique_candidate_peps_pseudo_kog_mnsod.csv')









