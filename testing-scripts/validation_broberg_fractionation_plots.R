library(Peptides)

bro_ms_frac1 <- bro_peps2 %>% 
  dplyr::filter(`Spectrum File` == '250516_Fraction_1.raw') %>%
  dplyr::select(pep_seq) %>% unique()

bro_ms_frac2 <- bro_peps2 %>% 
  dplyr::filter(`Spectrum File` == '250516_Fraction_2.raw') %>%
  dplyr::select(pep_seq) %>% unique()

bro_ms_frac3 <- bro_peps2 %>% 
  dplyr::filter(`Spectrum File` == '250516_Fraction_3.raw') %>%
  dplyr::select(pep_seq) %>% unique()

bro_ms_frac4 <- bro_peps2 %>% 
  dplyr::filter(`Spectrum File` == '250516_Fraction_4.raw') %>%
  dplyr::select(pep_seq) %>% unique()


bro_frac1 <- bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_1.raw', ]$pep_seq[complete.cases(bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_1.raw', ]$pep_seq)]
bro_frac2 <- bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_2.raw', ]$pep_seq[complete.cases(bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_2.raw', ]$pep_seq)]
bro_frac3 <- bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_3.raw', ]$pep_seq[complete.cases(bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_3.raw', ]$pep_seq)]
bro_frac4 <- bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_4.raw', ]$pep_seq[complete.cases(bro_peps2[bro_peps2$`Spectrum File` == '250516_Fraction_4.raw', ]$pep_seq)]


par(mfrow = c(4, 1))
hydrophobicity(seq = bro_frac1 %>% unique()) %>% hist(xlim = c(-4, 4), ylim = c(0, 40), breaks = 100, main = 'High pH Fraction 1')
hydrophobicity(seq = bro_frac2 %>% unique()) %>% hist(xlim = c(-4, 4), ylim = c(0, 40), breaks = 100, main = 'High pH Fraction 2')
hydrophobicity(seq = bro_frac3 %>% unique()) %>% hist(xlim = c(-4, 4), ylim = c(0, 40), breaks = 100, main = 'High pH Fraction 3')
hydrophobicity(seq = bro_frac4 %>% unique()) %>% hist(xlim = c(-4, 4), ylim = c(0, 40), breaks = 100, main = 'High pH Fraction 4', xlab = 'Hydrophobicity')





