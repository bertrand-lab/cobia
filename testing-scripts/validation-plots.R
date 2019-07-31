# plotting aggregated validation results

library(gridExtra)
# to do

# change font size of coef plot
# add in number of peps
# add ggtitle for plot F

# Making 

extract_coef <- function(glm_model, char_id){
  glm_model_summary <- summary(glm_model)
  glm_coef <- coefficients(glm_model_summary) %>% as.data.frame()
  X0_coef <- glm_coef$Estimate[2]
  X0_se <- glm_coef$`Std. Error`[2]
  return(data.frame(coef = X0_coef, se = X0_se, char_id = rep(char_id, 1)))
}

coef_diagram <- rbind(extract_coef(glm_kleiner260_oligo, 'A'),
                      extract_coef(glm_kleiner460_oligo, 'B'),
      extract_coef(glm_broberg_mg, 'C'),
      extract_coef(glm_broberg_mt, 'D'),
      extract_coef(glm_aylward_oligo, 'E'))

coef_graph <- ggplot(coef_diagram, aes(x = char_id, y = coef)) + 
  geom_point(alpha = 0.5, size = 2) +
  coord_flip() +
  geom_linerange(aes(ymin = coef - se, ymax = coef + se, x = char_id)) +
  geom_hline(yintercept = 0, lty = 2) + xlab("") + ylab("Coefficient Estimate") + 
  # ylim(c(-0.25, 0.25)) + 
  theme_bw() + 
  ggtitle('F) Summaries of Validation by Dataset') + 
  theme(legend.position = "none", 
        title = element_text(size = 7.5)) +
  scale_x_discrete(limits = rev(levels(coef_diagram$char_id)));coef_graph


jpeg("validation_plots.jpeg", width=170, height=210, units="mm", res=850)

grid.arrange(kleiner_260_p1, kleiner_460_p2, broberg_mg_p3, broberg_mt_p4, aylward_p5_oligo, coef_graph, nrow = 3)

dev.off()
