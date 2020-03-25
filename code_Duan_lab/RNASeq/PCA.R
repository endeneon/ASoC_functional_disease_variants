library(edgeR)
library(ggplot2)
library(gplots)
library(factoextra)




#############
k_cpm <- 1
k_sum <- 7



group <- factor(c(rep("CN", 8), rep("DN", 8), rep("iPS", 8), rep("NPC", 8), rep("GA", 8)))

# main_DGE <- DGEList(counts = raw_reads, group = group, remove.zeros = T) #?
main_DGE <- DGEList(counts = raw_reads, remove.zeros = T) #?

all_shared_genes <- as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 1:8]) > k_cpm)) > k_sum, , 
                                                    keep.lib.size = FALSE]$counts), stringsAsFactors = F)

all_shared_genes <- as.data.frame(merge(all_shared_genes, as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 9:16]) > k_cpm)) > k_sum, , 
                                                              keep.lib.size = FALSE]$counts), stringsAsFactors = F), 
                          by.x = colnames(all_shared_genes), 
                          by.y = colnames(as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 9:16]) > k_cpm)) > k_sum, , 
                                                                 keep.lib.size = FALSE]$counts), stringsAsFactors = F)), 
                          all.x = T, all.y = T))
all_shared_genes <- as.data.frame(merge(all_shared_genes, as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 17:24]) > k_cpm)) > k_sum, , 
                                                                            keep.lib.size = FALSE]$counts), stringsAsFactors = F), 
                          by.x = colnames(all_shared_genes), 
                          by.y = colnames(as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 17:24]) > k_cpm)) > k_sum, , 
                                                                          keep.lib.size = FALSE]$counts), stringsAsFactors = F)), 
                          all.x = T, all.y = T))
all_shared_genes <- as.data.frame(merge(all_shared_genes, as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 25:32]) > k_cpm)) > k_sum, , 
                                                                            keep.lib.size = FALSE]$counts), stringsAsFactors = F), 
                          by.x = colnames(all_shared_genes), 
                          by.y = colnames(as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 25:32]) > k_cpm)) > k_sum, , 
                                                                          keep.lib.size = FALSE]$counts), stringsAsFactors = F)), 
                          all.x = T, all.y = T))
all_shared_genes <- as.data.frame(merge(all_shared_genes, as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 33:40]) > k_cpm)) > k_sum, , 
                                                                            keep.lib.size = FALSE]$counts), stringsAsFactors = F), 
                          by.x = colnames(all_shared_genes), 
                          by.y = colnames(as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 33:40]) > k_cpm)) > k_sum, , 
                                                                          keep.lib.size = FALSE]$counts), stringsAsFactors = F)), 
                          all.x = T, all.y = T))


main_DGE <- main_DGE[all_shared_genes[[1]], , 
                     keep.lib.size = F]
main_DGE <- calcNormFactors(main_DGE)
main_DGE <- estimateDisp(main_DGE)
design <- model.matrix(~group)
main_DGE_fitted <- glmQLFit(main_DGE, design = design)
# main_DGE_fitted <- main_DGE

####################
gene_length_FPKM <- as.vector(merge(all_shared_genes, read_length, by.x = colnames(all_shared_genes), by.y = "Geneid"))

# FPKM_PCA <- rpkm(main_DGE_fitted$fitted.values, gene.length = gene_length_FPKM$Length,
#                  log = T)
FPKM_PCA <- rpkm.DGEList(main_DGE, gene.length = gene_length_FPKM$Length,
log = T)
# filtered_cpm <- cpm.DGEList(main_DGE_fitted, log = T)
# FPKM_PCA <- rpkm(main_DGE_fitted$counts, gene.length = gene_length_FPKM$Length, 
#                  log = T)

# FPKM_PCA <- FPKM_PCA[ , -(9:16)]
###################

p <- prcomp(FPKM_PCA, center = TRUE, scale. = TRUE, tol = 0) # why not all samples??
# p <- princomp(FPKM_PCA)chr19:49,419,503-49,445,098
# p <- prcomp(FPKM_PCA) # why not all samples??
# p <- prcomp(filtered_cpm, tol = 0)

PCA_data <- data.frame(x = 2, y = 1:ncol(filtered_cpm))

PCA_data$x <- p$rotation[, 1]
PCA_data$y <- p$rotation[, 2]
PCA_data$z <- p$rotation[, 3]
# 
# PCA_data$x <- p$loadings[, 1]
# PCA_data$y <- p$loadings[, 2]



scr_p <- fviz_eig(p)

rownames(PCA_data) <- colnames(FPKM_PCA)
PCA_data$type <- NA
# PCA_data$type <- c(rep("GA_08_Hanwen", 4), rep("GA_RNASeq", 5), rep("NSC", 8), rep("Yan_GA", 2), rep("CN", 8),
#                    rep("iPS", 8), rep("GA_09_Hanwen", 4))
PCA_data$type <- c(rep("Glut", 8), rep("Dopa", 8), rep("iPS", 8), rep("NPC", 8), rep("GA", 8))
PCA_data$type <- c(rep("red", 8), rep("orange", 8), rep("blue", 8), rep("magenta", 8), rep("green", 8))
# PCA_data$type <- c(rep("Glut(n=8)", 8), rep("iPS(n=8)", 8), rep("NPC(n=8)", 8), rep("GA_08(n=4)", 4), rep("GA_09(n=4)", 4), rep("GA_Yan(n=2)", 2))
# plot(PCA_data)
# text(PCA_data, labels = rownames(PCA_data))
# rownames(PCA_data) <- PCA_data$type
# rownames(PCA_data) <- c("NPC_07", "NPC_08", "NPC_09", "NPC_11", "NPC_12", "NPC_16", "NPC_18", "NPC_21", 
# "iPS_07", "iPS_08", "iPS_09", "iPS_11", "iPS_12", "iPS_16", "iPS_18", "iPS_21")
# rownames(PCA_data) <- as.vector(colnames())
ggplot(PCA_data, aes(x=x, y=y, label = rownames(PCA_data))) +
  geom_point(aes(colour = PCA_data$type)) +
  #  geom_text(angle = 45, colour = "gray") + 
  theme_minimal() +
  labs(x = "PC1", y = "PC2") +
 labs(x = paste(round(scr_p$data$eig[1], 2), "%"), y = paste(round(scr_p$data$eig[2], 2), "%")) +
  guides(color=guide_legend(title="Cell Type")) #+
  xlim(-0.19, -0.16) + ylim(-0.17, -0.02)#

#############
  library(rgl)
  #plot3d(x, y, z, col="red")
  plot3d(PCA_data[, 1:3], type = "s", col = PCA_data$type, size = 1, top = T, box = F, alpha = 0.7) +
    plot3d(PCA_data[, 1:3], type = "h", col = PCA_data$type, size = 1, top = T, box = F, 
           add = T, alpha = 0.5, lty = 2) +
    # rgl.texts(PCA_data[, 1:3], text = c("abc", "bcd", "def")) +
    # axis3d('xx', labels = paste("PC1: ", round(scr_p$data$eig[1], 2), "%")) +
    # axis3d('y', labels = paste("PC2:", round(scr_p$data$eig[2], 2), "%")) +
    # axiz3d('z', labels = paste("PC3:", round(scr_p$data$eig[3], 2), "%")) +
    axis3d(edge = (c('x', 'y', 'z'))) +
    title3d(
      xlab = paste("PC1: ", round(scr_p$data$eig[1], 2), "%"),
      ylab = paste("PC2:", round(scr_p$data$eig[2], 2), "%"),
      zlab = paste("PC3:", round(scr_p$data$eig[3], 2), "%")) +
    # text3d()
    planes3d(a = 0, b = 0, c = 1, d = 0, color = "cyan", alpha = 0.3, fog = T, emission = "cyan") #+
  abclines3d(0.15, 0, 0, a = diag(3), col = "black", lty = 4)
  
