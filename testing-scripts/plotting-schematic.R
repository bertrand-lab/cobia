library(truncnorm)
library(ggplot2)
library(magrittr)
library(ggforce)

set.seed(33331)

rt_xmin <- runif(500, min = 10, max = 40)
rt_xmax <- rt_xmin + 3
mz_ymin <- rtruncnorm(n = 500, a = 300, b = 1000, mean = 500, sd = 300)
mz_ymax <- mz_ymin + 3

df <- data.frame(rt_xmin, rt_xmax, mz_ymax, mz_ymin)

p1 <- df %>% 
  ggplot() + 
  geom_rect(mapping = aes(xmin = rt_xmin, 
                          xmax = rt_xmax, 
                          ymin = mz_ymin, 
                          ymax = mz_ymax), alpha = 0.4) +
  ylab('m/z') +
  xlab('Retention Time') +
  theme_bw();p1


p <- df %>% 
  ggplot() + 
  geom_rect(mapping = aes(xmin = rt_xmin, 
                          xmax = rt_xmax, 
                          ymin = mz_ymin, 
                          ymax = mz_ymax), alpha = 0.4) +
  ylab('m/z') +
  xlab('Retention Time') +
  theme_bw() + 
  facet_zoom(x = rt_xmin > 29 & rt_xmax < 35, 
             y = mz_ymin > 370 & mz_ymin < 420, zoom.size = 0.5, show.area = FALSE) + 
  annotate('text',label = 'Low Risk \nPeptide', y = 405, x = 33) +
  annotate('text',label = 'High Risk \nPeptide', y = 380, x = 30.5) +
  geom_segment(aes(x = 33, 
                   y = 407, 
                   xend = 33, 
                   yend = 412), 
               colour='black', 
               size = 1,
               arrow = arrow(length = unit(0.3, "cm"), type = 'closed')) +
  geom_segment(aes(x = 31, 
                   y = 377, 
                   xend = 32, 
                   yend = 373),
               lineend = 'butt',
               colour='black', 
               size = 1,
               # arrow = arrow(length = unit(0.3, "cm")));p
               arrow = arrow(length = unit(0.3, "cm"), type = 'closed')) +
  geom_rect(xmin = 29, xmax = 35, ymin = 370, ymax = 420, 
            fill = NA, colour = 'grey20');p



# manually adjusting the zoomed in figure
pb <- ggplot_build(p)
pb$data[[2]][1, 'alpha'] <- 0
pb$data[[3]][1, 'alpha'] <- 0
pb$data[[4]][1:1000, 'alpha'] <- rep(0, 1000)
pb$data[[5]][1:1000, 'alpha'] <- rep(0, 1000)
pb$data[[6]][1001:2000, 'xmin'] <- rep(0, 1000)
pb$data[[6]][1001:2000, 'ymin'] <- rep(0, 1000)
pb$data[[6]][1001:2000, 'ymax'] <- rep(1000, 1000)
pb$data[[6]][1001:2000, 'xmax'] <- rep(1000, 1000)
pg <- ggplot_gtable(pb)
plot(pg)


p2 <- df %>% 
  ggplot() + 
  theme_bw() +
  geom_rect(xmin = 29, xmax = 35, ymin = 370, ymax = 420, 
            fill = NA, colour = 'grey20') +
  geom_rect(mapping = aes(xmin = rt_xmin, 
                          xmax = rt_xmax, 
                          ymin = mz_ymin, 
                          ymax = mz_ymax), alpha = 0.4) +
  ylab('m/z') +
  xlab('Retention Time');p2

p2b <- df %>% 
  ggplot() +
  xlim(c(29, 35)) + ylim(c(370, 420)) +
  geom_rect(mapping = aes(xmin = rt_xmin, 
                          xmax = rt_xmax, 
                          ymin = mz_ymin, 
                          ymax = mz_ymax), alpha = 0.4) +
  theme_bw() +
  annotate('text',label = 'Low Risk \nPeptide', y = 405, x = 33) +
  annotate('text',label = 'High Risk \nPeptide', y = 380, x = 30.5) +
  geom_segment(aes(x = 33, 
                   y = 407, 
                   xend = 33, 
                   yend = 412), 
               colour='black', 
               size = 1,
               arrow = arrow(length = unit(0.3, "cm"), type = 'closed')) +
  geom_segment(aes(x = 31, 
                   y = 377, 
                   xend = 32, 
                   yend = 373),
               lineend = 'butt',
               colour='black', 
               size = 1,
               # arrow = arrow(length = unit(0.3, "cm")));p
               arrow = arrow(length = unit(0.3, "cm"), type = 'closed'));p2b


