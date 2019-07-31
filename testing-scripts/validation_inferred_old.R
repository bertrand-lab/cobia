# second validation

library(dplyr)
library(magrittr)
library(ggplot2)

bro_co_mg <- read.csv("data/broberg_data/broberg_mg_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")
bro_co_mt <- read.csv("data/broberg_data/broberg_mt_cofrag_mi-0.00833333_ipw-0.725_para-15_co-sim.csv")

bro_ms <- readxl::read_xlsx("data/broberg_data/Supplementary_Table_33.xlsx", skip = 4)

# observed proteins (from protein inference) digested into peptides
obs_protein <- read.table('data/broberg_data/mp_AT4_trypsin_digest.txt')

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

bro_co_mg2 <- cofrag_processor(cofrag_output_df = bro_co_mg, 
                               observed_peps = bro_peps4$pep_seq)
bro_co_mt2 <- cofrag_processor(cofrag_output_df = bro_co_mt, 
                               observed_peps = bro_peps4$pep_seq)

####
names(obs_protein) <- 'pep_seq'

obs_protein$ms_observed <- obs_protein$pep_seq %in% bro_peps4$pep_seq

obs_protein2 <- inner_join(x = obs_protein, y = bro_co_mg2, by = 'pep_seq')
obs_protein3 <- inner_join(x = obs_protein, y = bro_co_mt2, by = 'pep_seq')

obs_protein2$rescale_cofrag_score <- scales::rescale(obs_protein2$mean_cofrag_score, to = c(0, 100))
obs_protein3$rescale_cofrag_score <- scales::rescale(obs_protein3$mean_cofrag_score, to = c(0, 100))

number_of_unseen(obs_protein2)
mginferred_n <- "italic(N)[italic(unobserved)] == 4773"
mginferred_n_sub <- "italic(N)[italic(observed)] == 250"

mg_inferred <- obs_protein2 %>% ggplot(aes(x = mean_cofrag_score)) + 
  geom_histogram(aes(fill = ms_observed), binwidth = 1) +
  theme_bw() +
  ggtitle('A) Diseased Oak Tree Metagenome') +
  ylab('log(Count)') + xlab('Cofragmentation Score')  +
  annotate("text", x = 75, y = Inf, vjust = 3, label = mginferred_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 75, y = Inf, vjust = 4.5, label = mginferred_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + scale_y_log10();mg_inferred
# obs_protein2 %>% ggplot(aes(x = rescale_cofrag_score)) + 
#   geom_histogram(aes(fill = ms_observed), binwidth = 1) +
#   theme_bw() +
#   ggtitle('Cofragmentation histogram of only peptides \ntheoretically present after protein inference. Scores calculated from metagenome.')




number_of_unseen(obs_protein3)
mtinferred_n <- "italic(N)[italic(unobserved)] == 12612"
mtinferred_n_sub <- "italic(N)[italic(observed)] == 586"
mt_inferred <- obs_protein3 %>% ggplot(aes(x = mean_cofrag_score)) + 
  geom_histogram(aes(fill = ms_observed), binwidth = 1) +
  theme_bw() +
  ggtitle('B) Diseased Oak Tree Metatranscriptome') +
  ylab('log(Count)') + xlab('Cofragmentation Score')  +
  annotate("text", x = 250, y = Inf, vjust = 3, label = mtinferred_n, parse = TRUE, size = 3, colour = 'grey62') +
  annotate("text", x = 250, y = Inf, vjust = 4.5, label = mtinferred_n_sub, parse = TRUE, size = 3, colour = 'grey24') +
  scale_fill_manual(values = c("grey62", "grey24")) +
  theme(legend.position = "none",
        title = element_text(size = 7.5)) + scale_y_log10();mt_inferred
# obs_protein3 %>% ggplot(aes(x = rescale_cofrag_score)) + 
#   geom_histogram(aes(fill = ms_observed), binwidth = 1) +
#   theme_bw() +
#   ggtitle('Cofragmentation histogram of only peptides \ntheoretically present after protein inference. Scores calculated from metatranscriptome.')

glm_peptide_restricted_mg <- glm(ms_observed ~ mean_cofrag_score, data = obs_protein2, family = 'binomial')
glm_peptide_restricted_mt <- glm(ms_observed ~ mean_cofrag_score, data = obs_protein3, family = 'binomial')

summary(glm_peptide_restricted_mg)
summary(glm_broberg_mg)

summary(glm_peptide_restricted_mt)
summary(glm_broberg_mt)

coef_diagram_inferred <- rbind(extract_coef(glm_peptide_restricted_mg, 'Diseased Oak Tree Metagenome, \nInferred Peptides'),
                               extract_coef(glm_broberg_mg, 'Diseased Oak Tree Metagenome, \nAll Potential Peptides'),
                               extract_coef(glm_peptide_restricted_mt, 'Diseased Oak Tree Metatranscriptome, \nInferred Peptides'),
                               extract_coef(glm_broberg_mt, 'Diseased Oak Tree Metatranscriptome, \nAll Potential Peptides'))

coef_diagram_inferred$inferred <- c('Inferred', 'All', 'Inferred', 'All')
coef_diagram_inferred$dataset <- c('Diseased Oak Tree Metagenome', 'Diseased Oak Tree Metagenome', 'Diseased Oak Tree Metatranscriptome', 'Diseased Oak Tree Metatranscriptome')

coef_diagram_inferred$dataset <- factor(coef_diagram_inferred$dataset, c("Diseased Oak Tree Metatranscriptome", "Diseased Oak Tree Metagenome"))

coef_graph_inferred <- ggplot(coef_diagram_inferred, aes(x = dataset, y = coef)) + 
  coord_flip() +
  geom_linerange(aes(ymin = coef - se, ymax = coef + se, x = dataset, colour = inferred), alpha = 0.7, 
                 position = position_dodge(0.2)) +
  geom_point(alpha = 0.5, position = position_dodge(0.2), aes(x = dataset, colour = inferred)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab("") + ylab("Coefficient Estimate") + 
  # ylim(c(-0.25, 0.25)) + 
  theme_bw() + 
  # facet_wrap(~inferred, nrow = 2)+
  # ggtitle('F) Summaries of Validation by Dataset') +
  scale_colour_manual(values = c('grey60', 'grey20')) +
  theme(legend.title = element_blank());coef_graph_inferred
# scale_x_discrete(limits = rev(levels(coef_diagram_inferred$char_id)));coef_graph_inferred


jpeg("figures/main_validation_plot_inferred.jpeg", width=170, height=210, units="mm", res=850)

grid.arrange(mg_inferred, mt_inferred, coef_graph_inferred)

dev.off()

