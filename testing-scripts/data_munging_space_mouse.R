library(dplyr)
library(magrittr)
library(ggplot2)

# data munging for the space mouse proteome

space_mouse <- read.table('data/space_mouse_data/allPeptides.txt', sep = '\t')

space_mouse_peptides <- as.character(space_mouse$V23)
space_mouse_rts <- as.numeric(as.character(space_mouse$V16))*60

space_mouse_df <- data.frame(space_mouse_peptides, space_mouse_rts)
space_mouse_df_2 <- space_mouse_df[-1,]

space_mouse_df_3 <- space_mouse_df_2[which(space_mouse_df_2$space_mouse_peptides != " "), ]

space_mouse_df_sub <- sample_n(space_mouse_df_3, 1000)

write.table(space_mouse_df_sub, file = 'data/space_mouse_data/space_mouse_peptides_formatted.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(space_mouse_df_3, file = 'data/space_mouse_data/space_mouse_peptides_formatted_all.txt', 
            col.names = FALSE, sep = "\t",
            row.names = FALSE, quote = FALSE)
