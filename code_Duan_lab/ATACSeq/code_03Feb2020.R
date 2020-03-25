# 03 Feb 2020
# Siwei add microglia samples
# use DN omni 8 set

library(DiffBind)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(factoextra)
library(ggplot2)


summit_binding_matrix <- as.data.frame(cbind(summit_binding_matrix_not_DN_omni, 
                                             Xin_666614_summits_microglia[, 7:10]))

k_signal <- 10
k_sum <- 4

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

all_shared_genes <- as.vector(union(all_shared_genes, 
                                    as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 40:44] > 
                                                                                        median(unlist(summit_binding_matrix[, 40:44])))) 
                                                                             > k_sum, ]))
))

plot_summit_binding_matrix <- summit_binding_matrix[rownames(summit_binding_matrix) %in% all_shared_genes, ]
# plot_summit_binding_matrix[, 1:3] <- NULL
# plot_summit_binding_matrix <- summit_binding_matrix

group <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8), rep("MG", 4))
main_DGE <- DGEList(counts = plot_summit_binding_matrix, remove.zeros = T, group = group)

# FPKM_PCA <- rpkm.DGEList(main_DGE, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
#                          log = T)

FPKM_PCA <- rpkm.DGEList(main_DGE, gene.length = rep(501, length.out = length(all_shared_genes)),
                         log = T)

###################

p <- prcomp(FPKM_PCA, center = TRUE, scale. = TRUE, tol = 0) # why not all samples??


PCA_data <- data.frame(x = 2, y = 1:ncol(plot_summit_binding_matrix))

PCA_data$x <- p$rotation[, 1]
PCA_data$y <- p$rotation[, 2]
PCA_data$z <- p$rotation[, 3]
# 
# PCA_data$x <- p$loadings[, 1]
# PCA_data$y <- p$loadings[, 2]



scr_p <- fviz_eig(p)

rownames(PCA_data) <- colnames(FPKM_PCA)
PCA_data$type <- NA

PCA_data$type <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8), rep("MG", 4))

ggplot(PCA_data, aes(x=x, y=y, label = rownames(PCA_data))) +
  geom_point(aes(colour = PCA_data$type)) +
  #  geom_text(angle = 45, colour = "gray") + 
  theme_minimal() +
  labs(x = "PC1", y = "PC2") +
  labs(x = paste(round(scr_p$data$eig[1], 2), "%"), y = paste(round(scr_p$data$eig[2], 2), "%")) +
  guides(color=guide_legend(title="Cell Type")) #+
# xlim(-0.19, -0.16) + ylim(-0.17, -0.02)#

heatmap.2(cor(FPKM_PCA), 
          density.info = "none", trace = "none", dendrogram = "row",
          col = colorRampPalette(c("blue4", "darkblue", "blue", "green", "yellow", "red", "darkred"))(n=800))#,
breaks = c(seq(0.4, 0.50, length.out = 130),
           seq(0.501, 0.60, length.out = 130),
           seq(0.601, 0.65, length.out = 130),
           seq(0.651, 0.78, length.out = 130),
           seq(0.781, 0.85, length.out = 130),
           seq(0.851, 0.99, length.out = 130),
           seq(0.991, 1, length.out = 21)))

heatmap.2(cor(FPKM_PCA), 
          density.info = "none", trace = "none", dendrogram = "row",
          col = colorRampPalette(c("blue4", "darkblue", "blue", "green", "yellow", "red", "darkred"))(n=800))#,
