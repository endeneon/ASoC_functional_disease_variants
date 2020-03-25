# 23 Jan 2019
# plot with postmortem ATAC_seq results
library(ggplot2)
library(factoextra)


summit_binding_matrix <- DiffBind_read_counted_summit$binding
summit_binding_matrix <- as.data.frame(summit_binding_matrix)

summit_binding_matrix[, 1:3] <- NULL

#############
# 
# hist(unlist(summit_binding_matrix[, 1:8]), xlim = c(0,30), breaks = 3000)
# median(unlist(summit_binding_matrix[, 1:8]))
# median(unlist(summit_binding_matrix[, 9:16]))
# median(unlist(summit_binding_matrix[, 17:24]))
# median(unlist(summit_binding_matrix[, 25:32]))
# median(unlist(summit_binding_matrix[, 33:40]))
# 
# hist(unlist(summit_binding_matrix[, 25:32]), xlim=c(0,30), breaks = 1000)

# summit_binding_matrix <- cbind(summit_binding_matrix, Xin_250_summits_DNAseI_count_18Jan2018[, 7:85])

# summit_binding_matrix <- cbind(summit_binding_matrix, Xin_250_all_hg38_DNAseI_organoid_ATAC_21Jan2019[, 7:40])
summit_binding_matrix <- cbind(summit_binding_matrix, all_hg38_DNAseI_organoid_ATAC_22Jan2019[, 7:44]) # 22Jan, no organoid
summit_binding_matrix$Geneid <- Xin_original_summits_250_05Jan2019$X4

summit_binding_matrix <- merge(summit_binding_matrix, counts, 
                               by.x = "Geneid", by.y = "Geneid_hg38")
rownames(summit_binding_matrix) <- summit_binding_matrix$Geneid
summit_binding_matrix$Geneid <- NULL

rownames(summit_binding_matrix[is.na(rowSums(summit_binding_matrix)) == T], )
colnames(apply(summit_binding_matrix, 2, function(x) any(is.na(x))))


########
k_signal <- 10
k_sum <- 6

########
all_shared_genes <- as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 1:8] > 
                                                                        median(unlist(summit_binding_matrix[, 1:8])))) > k_sum, ]))

# all_shared_genes_temp <- as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16] > 
#                                                                              median(unlist(summit_binding_matrix[, 9:16])))) > k_sum, ]))

all_shared_genes <- as.vector(union(all_shared_genes, 
                                    as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16] > 
                                                                                        median(unlist(summit_binding_matrix[, 9:16])))) 
                                                                             > k_sum, ]))
))

all_shared_genes <- as.vector(union(all_shared_genes, 
                                    as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 17:24] > 
                                                                                        median(unlist(summit_binding_matrix[, 17:24])))) 
                                                                             > k_sum, ]))
))

all_shared_genes <- as.vector(union(all_shared_genes, 
                                    as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 25:32] > 
                                                                                        median(unlist(summit_binding_matrix[, 25:32])))) 
                                                                             > k_sum, ]))
))

all_shared_genes <- as.vector(union(all_shared_genes, 
                                    as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 33:40] > 
                                                                                        median(unlist(summit_binding_matrix[, 33:40])))) 
                                                                             > k_sum, ]))
))

plot_summit_binding_matrix <- summit_binding_matrix[rownames(summit_binding_matrix) %in% all_shared_genes, ]

# plot_summit_binding_matrix <- summit_binding_matrix
# plot_summit_binding_matrix[, colnames(plot_summit_binding_matrix) == "ENCFF982IRZ.bam"] <- NULL


group <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8), 
           rep("Brain_DNAseI", 38), rep("tel_organoid_ATAC", 3), rep("Brain_Neuro", 115))#, rep("tel_organoid_ATAC", 3))
group[grepl(pattern = "_G_", x = colnames(plot_summit_binding_matrix))] <- "Brain_non_Neuro"



main_DGE <- DGEList(counts = plot_summit_binding_matrix, remove.zeros = T, group = group)

main_DGE <- DGEList(counts = plot_summit_binding_matrix, remove.zeros = T)

system.time({
  main_DGE <- calcNormFactors(main_DGE)
})  

# sum(colSums(main_DGE$counts))


system.time({ 
  main_DGE <- estimateCommonDisp(main_DGE)
})

# 
# FPKM_PCA <- rpkm.DGEList(main_DGE, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
#                          log = T)

# FPKM_PCA <- rpkm(main_DGE$pseudo.counts, gene.length = rep(500, length.out = length(all_shared_genes)),
#                          log = T)

FPKM_PCA <- rpkm(main_DGE$pseudo.counts, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
                 log = T)

# test remove part of high/low columns
###################
hist(rowSums(FPKM_PCA), 50)
hist(rowSds(FPKM_PCA), 50)
hist(rowMedians(FPKM_PCA), 50)

FPKM_PCA <- FPKM_PCA[rowSds(FPKM_PCA) > 0.5, ]
FPKM_PCA <- FPKM_PCA[rowSds(FPKM_PCA) < 1.5, ]
# FPKM_PCA <- FPKM_PCA[rowMedians(FPKM_PCA) < 2, ]
# FPKM_PCA <- FPKM_PCA[rowMedians(FPKM_PCA) > -1, ]
###################


system.time({
  p <- prcomp(FPKM_PCA, center = TRUE, scale. = TRUE, tol = 0) # why not all samples??
})


system.time({
  p_alt <- prcomp(as.matrix(t(FPKM_PCA)), center = TRUE, scale. = TRUE, tol = 0) # yifan said do transpose
})

PCA_data <- data.frame(x = 2, y = 1:ncol(plot_summit_binding_matrix))

# PCA_data$x <- p$rotation[, 1]
# PCA_data$y <- p$rotation[, 2]
# PCA_data$z <- p$rotation[, 3]

PCA_data$x <- p$rotation[1, ]
PCA_data$y <- p$rotation[2, ]
PCA_data$z <- p$rotation[3, ]

PCA_data$x <- p$x[, 1]
PCA_data$y <- p$x[, 2]
PCA_data$z <- p$x[, 3]
# 
# PCA_data$x <- p$loadings[, 1]
# PCA_data$y <- p$loadings[, 2]



scr_p <- fviz_eig(p)

rownames(PCA_data) <- colnames(FPKM_PCA)
PCA_data$type <- NA

PCA_data$type <- c(rep("Glutamatergic Neuron", 8), 
                   rep("Dopaminergic Neuron", 8), 
                   rep("GABAergic Neuron", 8), 
                   rep("iPS", 8), 
                   rep("Neural Progenitor Cell", 8), 
                   rep("brain (fetal) DNAseI", 38), 
                   rep("Telecephalic Organoid", 3), 
                   rep("Postmortem Brain (neural)", 115)) #, rep("tel_organoid_ATAC", 3)
PCA_data$type[grepl(pattern = "_G_", x = rownames(PCA_data))] <- "Postmortem Brain (non-neural)"

df_tel_organoid <- PCA_data[79:81, ]
# plot_color <- #
ggplot(PCA_data, aes(x=x, y=y, label = rownames(PCA_data))) +
  geom_point(aes(colour = PCA_data$type, shape = PCA_data$type)) +
  #  geom_text(angle = 45, colour = "gray") + 
  theme_classic() +
  # labs(x = "PC1", y = "PC2") +
  labs(x = paste("PC1 =", round(scr_p$data$eig[1], 2), "%"), y = paste("PC2 =", round(scr_p$data$eig[2], 2), "%")) +
  scale_colour_manual(values = rep(c("seagreen", "darkorchid2", "magenta", "blue", 
                                     "brown", "darkred", "royalblue", "orange", "red"), length.out = 26)) +
  scale_shape_manual(values = rep(c(0,10,13,5,8,12,20,20,11,3,6), length.out = 26)) +
  labs(color = "Cell type", shape = "Cell type") +
  geom_text(data = df_tel_organoid, 
            aes(x = x, y = y, 
                label = c("day_0", "day_11", "day_30")), 
            nudge_x = 0.004, 
            nudge_y = -0.005, size = 4) +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10))


  geom_text(data = NULL, inherit.aes = F,
            aes(x = PCA_data[79:81, 1], 
                y = PCA_data[79:81, 2], 
                label = rownames(PCA_data)[79:81]))
  # guides(color=guide_legend(title="Cell Type")) #+
xlim(0.07, 0.14) +# 
  ylim(0.01, 0.05)#

############

heatmap.2(cor(FPKM_PCA[, 1:40]), 
          density.info = "none", trace = "none", dendrogram = "row",
          col = colorRampPalette(c("blue4", "darkblue", "blue", "green", "yellow", "red", "darkred"))(n=400),
          breaks = c(seq(0, 0.20, length.out = 65), # darkblue
                     seq(0.201, 0.402, length.out = 65), # blue
                     seq(0.403, 0.65, length.out = 65), # green
                     seq(0.651, 0.70, length.out = 65), # yellow
                     seq(0.701, 0.85, length.out = 65), # orange
                     seq(0.851, 0.99, length.out = 65), # red
                     seq(0.991, 1, length.out = 11)), 
          RowSideColors = c(rep("#F8766D", 8),
                            rep("#619CFF", 8),
                            
                            rep("#00BA38", 8), 
                            rep("#B79F00", 8),
                            rep("#00BFC4", 8) 
                            ))#,
                            # rep("#F564E3", 38))) # 
scales::show_col(scales::hue_pal()(6))
#######################

average_peak_size <- cleanup_Xin_original_peaks$X3 - cleanup_Xin_original_peaks$X2
mean(average_peak_size)
