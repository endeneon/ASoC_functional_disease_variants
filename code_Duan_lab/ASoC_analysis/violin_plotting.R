# make violin plots

library(ggplot2)
library(grid)
library(gridExtra)

p1 <- 
  ggplot(CN_20_ASoC_bed, aes(x = X3, y = (CN_20_ASoC_bed$X2), 
                             # color = CN_20_ASoC_bed$X3, 
                             fill = CN_20_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  # scale_y_continuous(trans ='log10') +
  labs(x = "Glut_20_non_ASoC", y = "-log2 p Value") +
  # guides(fill=guide_legend(title=NULL)) +
  scale_fill_discrete(name="Legend")
  # theme(legend.title = element_text())

p2 <- 
  ggplot(NSC_20_ASoC_bed, aes(x = X3, y = (NSC_20_ASoC_bed$X2), 
                             # color = NSC_20_ASoC_bed$X3, 
                             fill = NSC_20_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  # scale_y_continuous(trans ='log10') +
  labs(x = "NPC_20_non_ASoC", y = "-log2 p Value") +
  # guides(fill=guide_legend(title=NULL)) +
  scale_fill_discrete(name="Legend2")
  # theme(legend.title = "Legend")

grid_arrange_shared_legend(p1, p2,
                           nrow = 2)


#######


#######
p1 <- 
  ggplot(Glut_ASoC_bed, aes(x = X3, y = (Glut_ASoC_bed$X2), 
                             fill = Glut_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "Glutamatergic_non_ASoC", y = "-log2 p Value") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))

p2 <- 
  ggplot(NPC_ASoC_bed, aes(x = X3, y = (NPC_ASoC_bed$X2), 
                            fill = NPC_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "NPC_non_ASoC", y = "-log2 p Value") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))

p3 <- 
  ggplot(new_GA_ASoC_bed, aes(x = X3, y = (new_GA_ASoC_bed$X2), 
                           fill = new_GA_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "GABAergic_non_ASoC", y = "-log2 p Value") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))
p4 <- 
  ggplot(DN_ASoC_bed, aes(x = X3, y = (DN_ASoC_bed$X2), 
                           fill = DN_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "Dopaminergic_non_ASoC", y = "-log2 p Value") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))
p5 <- 
  ggplot(iPS_ASoC_bed, aes(x = X3, y = (iPS_ASoC_bed$X2), 
                           fill = iPS_ASoC_bed$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "iPS_non_ASoC", y = "-log2 p Value") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))

grid_arrange_shared_legend(p1, p2, p3, p4, p5, 
                           nrow = 5)



######
Glut_matched <- Glut_ASoC_bed[(order(Glut_ASoC_bed$X1) %in% order(cell_type_specific_ASoC$Glut)), ]
Glut_matched <- Glut_ASoC_bed[!((order(Glut_ASoC_bed$X1) %in% order(cell_type_specific_ASoC$Glut))), ]






NPC_matched <- NPC_ASoC_bed[NPC_ASoC_bed$X1 %in% cell_type_specific_ASoC$NPC, ]
DN_matched <- DN_ASoC_bed[DN_ASoC_bed$X1 %in% cell_type_specific_ASoC$DN, ]
GA_matched <- new_GA_ASoC_bed[new_GA_ASoC_bed$X1 %in% cell_type_specific_ASoC$GA, ]
iPS_matched <- iPS_ASoC_bed[iPS_ASoC_bed$X1 %in% cell_type_specific_ASoC$iPS, ]

# p1 <- #
  ggplot(Glut_matched, aes(x = X3, y = X2, 
                            fill = Glut_matched$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "Glutamatergic-specific ASoC SNPs", y = "-log10 p Value") +
  # scale_fill_discrete(name="Legend") +
  scale_fill_manual(name = "", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF")) +
    theme_classic() +
    theme(axis.text = element_text(size = 14), 
          axis.title = element_text(size = 14), 
          legend.text = element_text(size = 12))


p2 <- 
  ggplot(NPC_matched, aes(x = X3, y = (NPC_matched$X2), 
                           fill = NPC_matched$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "NPC_specific_ASoC", y = "-log2 p Value") +
  # scale_fill_discrete(name="Legend") +
  scale_fill_manual(name = "Legend", 
                    values = c("#FFBF00FF", "#00FF40FF",
                               "#00FFFFFF", "#8000FFFF", "#FF00BFFF", 
                               "#FFFFFF", "#FFFFFF", "#FFFFFF"))

p3 <- 
  ggplot(GA_matched, aes(x = X3, y = (GA_matched$X2), 
                              fill = GA_matched$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "GABAergic_specific_ASoC", y = "-log2 p Value") +
  # scale_fill_discrete(name="Legend")
  scale_fill_manual(values = c("#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF",
                               "#FFFFFF"))


p4 <- 
  ggplot(DN_matched, aes(x = X3, y = (DN_matched$X2), 
                          fill = DN_matched$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "Dopaminergic_specific_ASoC", y = "-log2 p Value") +
  # scale_fill_discrete(name="Legend")
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))


p5 <- 
  ggplot(iPS_matched, aes(x = X3, y = (iPS_matched$X2), 
                           fill = iPS_matched$X3)) +
  geom_violin(trim = F) +
  scale_y_continuous(trans ='log10') +
  labs(x = "iPS_specific_ASoC", y = "-log2 p Value") +
  # scale_fill_discrete(name="Legend")
  scale_fill_manual(name = "Legend", 
                    values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF",
                               "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF"))


NPC_matched <- as.data.frame(rbind(NPC_matched, NPC_matched[129, ]))
NPC_matched$X2 <- as.numeric(NPC_matched$X2)
######

grid_arrange_shared_legend(p1, p2, p3, p4, p5, 
                           nrow = 5)

# new_GA_ASoC_bed[1468, ] <- c("rsplot", 1.01, "3'-UTR")
new_GA_ASoC_bed <- as.data.frame(rbind(new_GA_ASoC_bed, new_GA_ASoC_bed[new_GA_ASoC_bed$X3 == "3'-UTR", ]))
new_GA_ASoC_bed$X2 <- as.numeric(new_GA_ASoC_bed$X2)
