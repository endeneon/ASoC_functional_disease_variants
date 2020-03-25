# Siwei 13 Mar 2019
# make volcano plot of CN20

nrow(Glut_20_anno_input_extended)

nrow(ASoC_w_pBinom)

library(ggplot2)
library(dplyr)

ASoC_w_pBinom$type <- "black"
ASoC_w_pBinom[ASoC_w_pBinom$FDR < 0.05, 11] <- "red"

ASoC_w_pBinom$REF_ratio <- ASoC_w_pBinom$REF/ASoC_w_pBinom$DP

sample_n(ASoC_w_pBinom, size = 10)
ASoC_w_pBinom_plot <- sample_n(ASoC_w_pBinom, size = 10000, replace = F)

ggplot(ASoC_w_pBinom_plot, aes(x = REF_ratio, y = ASoC_w_pBinom_plot$`-log10P`)) +

  geom_point(color = ASoC_w_pBinom_plot$type, cex = 0.5, alpha = 1) +
  stat_density_2d(aes(fill = 0-stat(density)), 
                  geom = "raster", alpha = 0.5, 
                  contour = F, 
                  n = c(1000,1000)) +
  scale_fill_gradientn(colours = c("red", "yellow", "green", "cyan", "white"), 
                       values = c(0, 0.25, 0.5,0.75, 1)) +

  theme_classic() +
  ylim(0, 20)


ggplot(ASoC_w_pBinom_plot, aes(x = REF_ratio, y = ASoC_w_pBinom_plot$`-log10P`)) +
  
  geom_point(color = ASoC_w_pBinom_plot$type, cex = 0.2, alpha = 1) +
  theme_classic() +
  ylim(0, 20) +
  labs(x = "Ratio of the reference allele", 
       y = "-log10 P value") +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8))
