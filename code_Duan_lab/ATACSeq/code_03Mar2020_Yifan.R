# 3 Mar 2020
# for Yifan
# exclude iPS from PCA calculation

FPKM_PCA <- FPKM_PCA[, -c(25:32)]

library(ggplot2)
library(factoextra)

system.time({
  p <- prcomp(FPKM_PCA, center = TRUE, scale. = TRUE, tol = 0) # why not all samples??
})


system.time({
  p_alt <- prcomp(as.matrix(t(FPKM_PCA)), center = TRUE, scale. = TRUE, tol = 0) # yifan said do transpose
})

PCA_data <- data.frame(x = 2, y = p_alt$x[2, ])

# PCA_data$x <- p$rotation[, 1]
# PCA_data$y <- p$rotation[, 2]
# PCA_data$z <- p$rotation[, 3]

PCA_data$x <- p$rotation[1, ]
PCA_data$y <- p$rotation[2, ]
PCA_data$z <- p$rotation[3, ]

PCA_data$x <- p_alt$x[1, ]
PCA_data$y <- p_alt$x[2, ]
PCA_data$z <- p_alt$x[3, ]

PCA_data$x <- p_alt$x[, 1]
PCA_data$y <- p_alt$x[, 2]
PCA_data$z <- p_alt$x[, 3]
# 
# PCA_data$x <- p$loadings[, 1]
# PCA_data$y <- p$loadings[, 2]



scr_p <- fviz_eig(p_alt)

rownames(PCA_data) <- colnames(FPKM_PCA)
PCA_data$type <- NA

PCA_data$type <- c(rep("Glutamatergic Neuron", 8), 
                   rep("Dopaminergic Neuron", 8), 
                   rep("GABAergic Neuron", 8), 
                   # rep("iPS", 8), 
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
                                     # "brown", 
                                     "darkred", "royalblue", "orange", "red"), 
                                   length.out = 26)) +
  scale_shape_manual(values = rep(c(0,10,13,5,
                                    # 8,
                                    12,20,20,11,3,6), length.out = 26)) +
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