# 19 Dec 2019
# Make plots similar to Fig S3 of Nott et al., Science 2019
# plot use cell type specific (not exclusive) peaks

# init
library(readr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(stringr)
library(purrr)
library(egg)
library(plyr)

### load data
main_motifs_list <- list()
main_motifs_list[["iN-Glut"]] <- read_delim("../specific_peaks/CN_specific_summit/knownResults.txt", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)
main_motifs_list[["iN-DN"]] <- read_delim("../specific_peaks/DN_specific_summit//knownResults.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)
main_motifs_list[["iN-GA"]] <- read_delim("../specific_peaks/GA_specific_summit/knownResults.txt", 
                                          "\t", escape_double = FALSE, trim_ws = TRUE)
main_motifs_list[["iPSC"]] <- read_delim("../specific_peaks/iPS_specific_summit/knownResults.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE)
main_motifs_list[["NPC"]] <- read_delim("../specific_peaks/NPC_specific_summit/knownResults.txt", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)

###

# main_motifs_list[["iPSC"]]$`Motif Name`[1] <- "OCT4"
# main_motifs_list[["iPSC"]] <- main_motifs_list[["iPSC"]][-4, ]
# main_motifs_list[["iPSC"]] <- main_motifs_list[["iPSC"]][-12, ]
# main_motifs_list[["iN-GA"]] <- main_motifs_list[["iN-GA"]][-9, ]

### plot top 25 motifs of each cell type
i <- 1
merged_cell_type_specific_plot <- list()
for (i in 1:length(main_motifs_list)) {
  ## enriched P value
  enriched_P_plot <- 
    ggplot() +
    scale_fill_distiller(palette = "Spectral", direction = 1) +
    geom_col(aes(x = main_motifs_list[[i]]$`Motif Name`[1:30], 
                 y = 0 - main_motifs_list[[i]]$`Log P-value`[1:30], 
                 fill = c(1:30))) +
    ylab("-log10P") +
    theme_classic() +
    ggtitle(label = names(main_motifs_list)[i]) +
    scale_x_discrete(limits =  main_motifs_list[[i]]$`Motif Name`[1:30],
                     labels = toupper(unlist(map(str_split(main_motifs_list[[i]]$`Motif Name`[1:30], 
                                                           pattern = '\\('), 1)))) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(), 
          legend.position = "none")
  
  ## % of targeted motifs
  targeted_perc_plot <- 
    ggplot() +
    scale_fill_distiller(palette = "Spectral", direction = 1L) +
    geom_col(aes(x = main_motifs_list[[i]]$`Motif Name`[1:30], 
                 y = as.numeric(str_sub(main_motifs_list[[i]]$`% of Target Sequences with Motif`[1:30], 
                                        end = -2L)), 
                 fill = c(1:30))) +
    ylab("%") +
    theme_classic() +
    scale_x_discrete(limits =  main_motifs_list[[i]]$`Motif Name`[1:30],
                     labels = toupper(unlist(map(str_split(main_motifs_list[[i]]$`Motif Name`[1:30], 
                                                           pattern = '\\('), 1)))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"), 
          legend.position = "none", 
          axis.title.x = element_blank())
  
  ## merge the two plots to one and align
  merged_cell_type_specific_plot[[i]] <- as_ggplot(egg::ggarrange(enriched_P_plot, targeted_perc_plot, 
                                                                  ncol = 1))
}

# ggarrange from egg package can make all plots aligned
egg::ggarrange(merged_cell_type_specific_plot[[4]], 
               merged_cell_type_specific_plot[[5]], 
               merged_cell_type_specific_plot[[1]], 
               merged_cell_type_specific_plot[[2]], 
               merged_cell_type_specific_plot[[3]], 
               ncol = 1, 
               left = 10, 
               heights = c(1.5, 1, 1, 1, 1))
####


test_string <- unlist(map(str_split(main_motifs_list[[i]]$`Motif Name`[1:30], pattern = '\\('), 1))
ggarrange(enriched_P_plot, targeted_perc_plot, ncol = 1)



data_top25motifs_list <- list()
data_top25motifs_list[["iPSC"]] <- main_motifs_list[[4]][1:30, ]
data_top25motifs_list[["NPC"]] <- main_motifs_list[[5]][1:30, ]
data_top25motifs_list[["iN-Glut"]] <- main_motifs_list[[1]][1:30, ]
data_top25motifs_list[["iN-DN"]] <- main_motifs_list[[2]][1:30, ]
data_top25motifs_list[["iN-GA"]] <- main_motifs_list[[3]][1:30, ]


rm(data_4_S_table_top25motifs)
data_4_S_table_top25motifs <- plyr::ldply(data_top25motifs_list, rbind)

save.image()
