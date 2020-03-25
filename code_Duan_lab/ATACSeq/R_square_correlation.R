# Siwei Jan 21 2019
# draw R square correlation

# export the DiffBind binding matrix, use signal threshold

library(edgeR)
library(gplots)
library(RColorBrewer)
library(factoextra)
library(ggplot2)
library(matrixStats) # for rowSds
library(dplyr)


# DiffBind_read_summit_DBA_READS <- dba.count(DiffBind_read_counted_summit, peaks = NULL, score = DBA_SCORE_READS)

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
# summit_binding_matrix <- cbind(summit_binding_matrix, all_hg38_DNAseI_organoid_ATAC_22Jan2019[, 7:44]) # 22Jan

### draw correlation density plot
summit_binding_matrix <- cbind(summit_binding_matrix, all_hg38_DNAseI_organoid_ATAC_22Jan2019[, 7:47]) # 22Jan

########
k_signal <- 10
k_sum <- 6

# all_shared_genes <- as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 1:8] > k_signal)) > k_sum, ]),
#                                   stringsAsFactors = F)
# 
# all_shared_genes <- as.data.frame(merge(all_shared_genes, 
#                                         as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16]) > k_signal) 
#                                                                                                        > k_sum, ]), 
#                                                                         stringsAsFactors = F), 
#                                         by.x = colnames(all_shared_genes), 
#                                         by.y = colnames(as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16]) > k_signal) 
#                                                                                                      > k_sum, ]), 
#                                                                       stringsAsFactors = F)), 
#                                         all.x = T, all.y = T))
###########
all_shared_genes <- as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 1:8] > 
                                                                        median(unlist(summit_binding_matrix[, 1:8])))) > k_sum, ]))

all_shared_genes_temp <- as.vector(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16] > 
                                                                             median(unlist(summit_binding_matrix[, 9:16])))) > k_sum, ]))

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

# plot_summit_binding_matrix <- summit_binding_mgrep "saf" ~/.bash_history 
# atrix
# plot_summit_binding_matrix[, colnames(plot_summit_binding_matrix) == "ENCFF982IRZ.bam"] <- NULL


group <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8), 
           rep("Brain_DNAseI", 38), rep("tel_organoid_ATAC", 3))

main_DGE <- DGEList(counts = plot_summit_binding_matrix, remove.zeros = T, group = group)

system.time({
  main_DGE <- calcNormFactors(main_DGE)
  main_DGE <- estimateCommonDisp(main_DGE)
})

####### Use log(pseudo counts within peaks +1)




# FPKM_PCA <- rpkm.DGEList(main_DGE, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
#                          log = T)

# FPKM_PCA <- rpkm(main_DGE$pseudo.counts, gene.length = rep(500, length.out = length(all_shared_genes)),
#                          log = T)

# FPKM_PCA <- rpkm(main_DGE$pseudo.counts, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
#                  log = T, prior.count = 1)
FPKM_PCA <- main_DGE$pseudo.counts
# FPKM_PCA <- rpkm(main_DGE, gene.length = rep(500, length.out = nrow(plot_summit_binding_matrix)),
#                  log = T)
FPKM_PCA <- as.data.frame(FPKM_PCA+1)
##########################

colnames(FPKM_PCA)

XY_data_frame <- as.data.frame(cbind(FPKM_PCA$GA_09, FPKM_PCA$ENCFF348VND.bam))
# colnames(XY_data_frame) <- c("GA_09", "ENCFF348VND.bam")
colnames(XY_data_frame) <- c("x", "y")

XY_data_frame <- XY_data_frame[XY_data_frame$x > 20, ]
XY_data_frame <- XY_data_frame[XY_data_frame$y > 20, ]
XY_data_frame <- log(sample_n(XY_data_frame, 50000), 2)

regression_data_frame <- sample_n(XY_data_frame, 500)

# 
# ggplot()+
#   geom_point(aes(x = XY_data_frame$GA_09, y = XY_data_frame$ENCFF348VND.bam),
#                  colour = "blue", alpha = 0.5) +
#   theme_minimal()
# 
# ggplot(x = XY_data_frame$GA_09, y = XY_data_frame$ENCFF348VND.bam)+
#   geom_bin2d(bins = c(100, 100)) +
#   theme_minimal()

ggplot(XY_data_frame, aes(x, y)) +

  geom_bin2d(bins = 500) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue4", "green", "yellow", "darkred"))(n=800), 
                       limits = c(0.5, 5)
                       # colours = c("red", "darkblue"), 
                       # breaks = c(seq(0, 100, length.out = 90)))
                       ) +
  geom_point(size = 0.5, colour = "gray", alpha = 0) +
  theme_minimal()

plot_f <-
  ggplot(XY_data_frame, aes(x, y)) +
  geom_bin2d(bins = 500) +
  scale_fill_gradientn(colours = colorRampPalette(c("lightblue", "steelblue", "royalblue", "darkblue"))(n=800), 
                       limits = c(0.5, 5)) +
                       # colours = c("red", "darkblue"), 
                       # breaks = c(seq(0, 100, length.out = 90)))
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  geom_point(size = 0.5, colour = "gray", alpha = 0) +
  theme_minimal() +
  geom_text(aes(x = 1, y = 5), label = lm_eqn(regression_data_frame), parse = TRUE)

  plot_f + 
    geom_text(aes(x = 1, y = 5), label = lm_eqn(regression_data_frame), parse = TRUE)

  ggplot() +
    geom_text(label = lm_eqn(regression_data_frame), parse = T, 
              aes(x = 1, y = 5, colour = "black"))
  
#############
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}




cor.test(XY_data_frame$GA_09, XY_data_frame$ENCFF348VND.bam)

##########################
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

PCA_data$type <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8), 
                   rep("Brain_DNAseI", 38)) #, rep("tel_organoid_ATAC", 3)

ggplot(PCA_data, aes(x=x, y=y, label = rownames(PCA_data))) +
  geom_point(aes(colour = PCA_data$type)) +
  #  geom_text(angle = 45, colour = "gray") + 
  theme_minimal() +
  labs(x = "PC1", y = "PC2") +
  labs(x = paste(round(scr_p$data$eig[1], 2), "%"), y = paste(round(scr_p$data$eig[2], 2), "%")) +
  guides(color=guide_legend(title="Cell Type")) #+
xlim(0.07, 0.14) +# 
  ylim(0.01, 0.05)#

############
DNAseI_to_use <- colnames(plot_summit_binding_matrix)[(PCA_data$y > 0) & (PCA_data$y < 0.05)]
DNAseI_to_use <- colnames(plot_summit_binding_matrix)[(PCA_data$x < 0.08)]  
############

heatmap.2(cor(FPKM_PCA), 
          density.info = "none", trace = "none", dendrogram = "row",
          col = colorRampPalette(c("blue4", "darkblue", "blue", "green", "yellow", "red", "darkred"))(n=800),
          breaks = c(seq(0.4, 0.50, length.out = 130),
                     seq(0.501, 0.60, length.out = 130),
                     seq(0.601, 0.65, length.out = 130),
                     seq(0.651, 0.78, length.out = 130),
                     seq(0.781, 0.85, length.out = 130),
                     seq(0.851, 0.99, length.out = 130),
                     seq(0.991, 1, length.out = 21)))


#######################

nrow(as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 1:8] > k_signal)) > k_sum, ]),
                   stringsAsFactors = F))

nrow(as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 9:16] > k_signal)) > k_sum, ]), 
                   stringsAsFactors = F))

sum(rownames(summit_binding_matrix[((rowSums(summit_binding_matrix[, 1:8] > k_signal)) > k_sum), ]) %in%
      rownames(summit_binding_matrix[((rowSums(summit_binding_matrix[, 9:16] > k_signal)) > k_sum), ]))


nrow(as.data.frame(rownames(summit_binding_matrix[(rowSums(summit_binding_matrix[, 25:32] > k_signal)) > k_sum, ]), 
                   stringsAsFactors = F))
max(rowSds(as.matrix(summit_binding_matrix)))
sum(rowSds(as.matrix(summit_binding_matrix)) == 0)

plot_summit_binding_matrix <- summit_binding_matrix[(rowSds(as.matrix(0-log(summit_binding_matrix))) > 1e-15), ]

heatmap.2(cor(0-log(plot_summit_binding_matrix)), 
          col = colorRampPalette(c("blue4", "darkblue", "blue", "green", "yellow", "red", "darkred"))(n=800))

heatmap(cor(0-log(plot_summit_binding_matrix, 2)))

max((rowSds(as.matrix(0-log(summit_binding_matrix)))))
hist((rowSds(as.matrix(0-log(summit_binding_matrix)))), breaks = 100)
###########

group <- c(rep("CN", 8), rep("DN", 8), rep("GA", 8), rep("iPS", 8), rep("NPC", 8))

# main_DGE <- DGEList(counts = raw_reads, group = group, remove.zeros = T) #?
main_DGE <- DGEList(counts = raw_reads, remove.zeros = T) #?

all_shared_genes <- as.data.frame(rownames(main_DGE[(rowSums(cpm(main_DGE$counts[, 1:8]) > k_cpm)) > k_sum, , 
                                                    keep.lib.size = FALSE]$counts), stringsAsFactors = F)