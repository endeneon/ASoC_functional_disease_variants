# Siwei 09 Mar 2020
# Make -logRatio vs logDP plot 
# both axis use log scale

# init
library(ggplot2)
library(dplyr)
library(RColorBrewer)


ASoC_w_pBinom$type <- "black"
ASoC_w_pBinom[ASoC_w_pBinom$FDR < 0.05, 11] <- "red"

ASoC_w_pBinom$REF_ratio <- ASoC_w_pBinom$REF/ASoC_w_pBinom$DP

ASoC_w_pBinom$neg_REF_ratio <- 0 - log2(ASoC_w_pBinom$REF_ratio)
ASoC_w_pBinom$log_DP <- log2(ASoC_w_pBinom$DP)


# sample_n(ASoC_w_pBinom, size = 10)
ASoC_w_pBinom_plot <- sample_n(ASoC_w_pBinom, size = 10000, replace = F)

# ggplot(ASoC_w_pBinom_plot, aes(x = neg_REF_ratio, y = log_DP)) +
#   
#   geom_point(color = ASoC_w_pBinom_plot$type, cex = 0.5, alpha = 1) +
#   stat_density_2d(aes(fill = 0-stat(density)), 
#                   geom = "raster", alpha = 0.5, 
#                   contour = F, 
#                   n = c(1000,1000)) +
#   scale_fill_gradientn(colours = c("red", "yellow", "green", "cyan", "white"), 
#                        values = c(0, 0.25, 0.5,0.75, 1)) +
#   
#   theme_classic() +
#   ylim(0, 20)


ggplot(ASoC_w_pBinom_plot, aes(x = neg_REF_ratio, y = DP)) +
  
  geom_point(color = ASoC_w_pBinom_plot$type, cex = 0.2, alpha = 1) +
  # theme_classic() +
  # ylim(0, 20) +
  labs(x = "-log2 value of the reference allele ratio", 
       y = "Read depth") +
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2') +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8)) +
  xlim(0, 2) +
  ylim(0, 2000) +
  theme_classic()

ggplot(ASoC_w_pBinom_plot, aes(x = `-log10P`, y = DP)) +
  
  geom_point(color = ASoC_w_pBinom_plot$type, 
             cex = 0.2, alpha = 1) +
  # theme_classic() +
  # ylim(0, 20) +
  labs(x = "-log2 value of the reference allele ratio", 
       y = "Read depth") +
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2') +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8)) +
  # xlim(0, 2) +
  # ylim(0, 2000) +
  theme_classic()


####

ASoC_w_pBinom$read_51bp <- Glut_20_100K_SNP_51bp_output$CN20_merged_sorted.bam
# sample_n(ASoC_w_pBinom, size = 10)
ASoC_w_pBinom_plot <- sample_n(ASoC_w_pBinom, size = 10000, replace = F)
ASoC_w_pBinom_plot$minRatio <- 1
ASoC_w_pBinom_plot$minRatio  <- sapply(1:nrow(ASoC_w_pBinom_plot), 
                                       function(x) min(ASoC_w_pBinom_plot$REF_ratio[x], 
                                                       (1 - ASoC_w_pBinom_plot$REF_ratio[x])))
ASoC_w_pBinom_plot$FDR_sig <- "Not significant"
ASoC_w_pBinom_plot$FDR_sig[ASoC_w_pBinom_plot$FDR < 0.05] <- "Significant"
ASoC_w_pBinom_plot$FDR_sig <- as.factor(ASoC_w_pBinom_plot$FDR_sig)

ggplot(ASoC_w_pBinom_plot, aes(x = minRatio, y = DP)) +
  geom_point(aes(colour = `-log10P`, 
                 shape = FDR_sig),
             # shape = ASoC_w_pBinom_plot$FDR_sig,
             cex = 0.5, alpha = 1) +
  # theme_classic() +
  # ylim(0, 20) +
  scale_shape_manual(values = c(19,23)) +
  # scale_fill_manual(values = c(""))
  labs(x = "Minor allele ratio", 
       y = "Read depth") +
  # scale_color_gradientn(colours = c("lightblue", 
  #                                   "red", "red", "red", 
  #                                   "red", "red", "red", "darkred")) +
  scale_color_gradientn(colours = c("lightskyblue", "orange", "red", "red"),
                        # values = c(0, 1, 5, 10, 20, 50), 
                        limits = c(0, 5),
                        na.value = "red") +
  # scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11, "matlab.like")))(100), 
  #                        limits = c(0, 10)) +
  # scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8)) +
  # xlim(0, 2) +
  # ylim(0, 2000) +
  # ggtitle(label = "Use read depth at SNP site, Y axis as log scale") +
  geom_vline(xintercept = 1/3, linetype = 2, colour = "darkblue") +
  theme_classic()


ASoC_w_pBinom$minRatio <- 1
ASoC_w_pBinom$minRatio <- sapply(1:nrow(ASoC_w_pBinom), 
                                  function(x) min(ASoC_w_pBinom$REF_ratio[x], 
                                                  (1 - ASoC_w_pBinom$REF_ratio[x])))


ggplot(ASoC_w_pBinom, aes(x = minRatio, y = DP)) +
  geom_point(aes(colour = as.factor(type)),
             cex = 0.5, alpha = 0.8,
             shape = 19) +
  # scale_shape_manual(values = c(19,23)) +
  labs(x = "Minor allele ratio", 
       y = "Read depth") +
  scale_colour_manual(values = c("black", "red")) +
  scale_y_continuous(trans = 'log2') +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8)) +
  geom_vline(xintercept = 1/3, linetype = 2, colour = "darkblue") +
  theme_classic()

### bar plot of the FDR-significant SNP in 50 increments, max DP = 5711
#####
DP_hist <- hist(ASoC_w_pBinom$DP, breaks = 40, xlim = c(20, 4020))
FDR_sig_hist <- hist(ASoC_w_pBinom$DP[ASoC_w_pBinom$FDR < 0.05], breaks = 40, xlim = c(20, 4020))

df_bar_plot <- data.frame(DP = DP_hist$counts, 
                          FDR_sig = FDR_sig_hist$counts)
df_bar_plot$non_sig <- df_bar_plot$DP - df_bar_plot$FDR_sig
df_bar_plot$index <- DP_hist$mids
df_bar_plot <- df_bar_plot[1:20, ]

df_bar_ggplot <- as.data.frame(rbind(data.frame(SNP_count = df_bar_plot$FDR_sig,
                                                SNP_index = df_bar_plot$index,
                                                SNP_colour = "sig"),
                                     data.frame(SNP_count = df_bar_plot$non_sig,
                                                SNP_index = df_bar_plot$index,
                                                SNP_colour = "non_sig")))

ggplot(df_bar_ggplot, aes(x = SNP_index, y = SNP_count,
                          fill = SNP_colour)) +
  geom_bar(stat = "identity", 
           position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("Read depth at SNP site") +
  ylab(label = "Categorial percentages of all SNP types in the bin") +
  # scale_y_log10() +
  theme_classic()

ggplot(df_bar_ggplot, aes(x = SNP_index, y = SNP_count,
                          fill = SNP_colour)) +
  geom_bar(stat = "identity", 
           position = "stack") +
  # scale_y_continuous(labels = scales::percent) +
  xlab("Read depth at SNP site") +
  ylab(label = "Categorial counts of all SNP types\nin the bin (log scale)") +
  # scale_y_log10() +
  theme_classic()
#####

ASoC_w_pBinom$sig <- 0
ASoC_w_pBinom$sig[ASoC_w_pBinom$FDR < 0.05] <- 1

ggplot(ASoC_w_pBinom, aes(x = DP, 
                          fill = sig)) +
  # geom_bar(stat = "bin") +
  stat_bin(binwidth = 50) +
  theme_classic()

###
df_sig_DP_count <- data.frame(DP = 1:250, 
                              sig_ASoC = 0, 
                              total_SNP = 0)
i <- 1
for (i in 1:250) {
  print(i)
  df_sig_DP_count$DP[i] <- i * 20
  df_sig_DP_count$sig_ASoC[i] <- sum(ASoC_w_pBinom$sig[ASoC_w_pBinom$DP >= (i * 20)])
  df_sig_DP_count$total_SNP[i] <- sum(ASoC_w_pBinom$DP >= (i * 20))
}
df_sig_DP_count$ASoC_perc <- df_sig_DP_count$sig_ASoC / df_sig_DP_count$total_SNP

ggplot(df_sig_DP_count, aes(x = DP, y = ASoC_perc)) +
  geom_point() +
  # xlim(0, 1000) +
  # ylim(0, 0.2) +
  # geom_vline(xintercept = c(100, 200), 
  #            colour = "red",
  #            show.legend = F) +
  scale_y_continuous(labels = scales::percent, 
                     limits = c(0, 0.2)) +
  labs(x = "Minimum OCR read depth at a SNP site", 
       y = "Detectable ASoC SNPs (FDR < 0.05) \nout of all heterozygous SNPs") +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000), 
                     limits = c(0, 1000)) +
  theme_classic() +
  # theme(axis.ticks.length.x = eleme)
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14))



#######

sig_ASoC_w_pBinom <- ASoC_w_pBinom[ASoC_w_pBinom$FDR < 0.05, ]
sig_ASoC_w_pBinom$minRatio <- 1
sig_ASoC_w_pBinom$minRatio <- sapply(1:nrow(sig_ASoC_w_pBinom), 
                                     function(x) min(sig_ASoC_w_pBinom$REF_ratio[x], 
                                                     (1 - sig_ASoC_w_pBinom$REF_ratio[x])))
sig_ASoC_w_pBinom$DP_threshold <- "20-100"
sig_ASoC_w_pBinom$DP_threshold[sig_ASoC_w_pBinom$DP > 99] <- "101-200"
sig_ASoC_w_pBinom$DP_threshold[sig_ASoC_w_pBinom$DP > 199] <- "> 201"
sig_ASoC_w_pBinom$DP_threshold <- as.factor(sig_ASoC_w_pBinom$DP_threshold)
mu <- plyr::ddply(sig_ASoC_w_pBinom, "DP_threshold", summarise, grp.mean = mean(minRatio))

ggplot(sig_ASoC_w_pBinom, aes(x = DP, y = `-log10P`, 
                              colour = minRatio)) +
  geom_point() +
  xlab("log10(DP)") +
  scale_color_gradientn(colours = c("red", "orange", "lightskyblue", "blue", "darkblue")
                        # values = c(0, 1, 5, 10, 20, 50), 
                        ) +
  scale_x_log10() +
  theme_classic()

ggplot(sig_ASoC_w_pBinom, aes(x = minRatio, 
                              colour = DP_threshold, 
                              fill = DP_threshold)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  geom_vline(data = mu, aes(xintercept = grp.mean, colour = DP_threshold), 
             linetype = "dashed") +
  labs(y = "ASoC SNP count") +
  theme_classic()

#
ggplot() +
  geom_histogram(data = sig_ASoC_w_pBinom,
                 aes(x = minRatio),
                 colour = "blue",
                 fill = "blue",
                 alpha = 0.1) +
  geom_histogram(data = sig_ASoC_w_pBinom[sig_ASoC_w_pBinom$DP > 99, ],
                 aes(x = minRatio),
                 colour = "green",
                 fill = "green",
                 alpha = 0.1) +
  geom_histogram(data = sig_ASoC_w_pBinom[sig_ASoC_w_pBinom$DP > 199, ],
                 aes(x = minRatio),
                 colour = "red",
                 fill = "red",
                 alpha = 0.1) +
  # geom_histogram(position = "identity", alpha = 0.5) +
  # geom_vline(data = mu, aes(xintercept = grp.mean, colour = DP_threshold), 
  #            linetype = "dashed") +
  labs(y = "ASoC SNP count") +
  ggtitle(label = "blue: all; green: DP>100; red: DP>200") +
  theme_classic()

df_minRatio <- as.data.frame(rbind(data.frame(minRatio = sig_ASoC_w_pBinom$minRatio, 
                                              index = "> 20"), 
                                   data.frame(minRatio = sig_ASoC_w_pBinom$minRatio[sig_ASoC_w_pBinom$DP > 99], 
                                              index = "> 100"), 
                                   data.frame(minRatio = sig_ASoC_w_pBinom$minRatio[sig_ASoC_w_pBinom$DP > 199], 
                                              index = "> 200")))

ggplot(df_minRatio, aes(x = minRatio, 
                        colour = index, 
                        fill = index)) +
  geom_histogram(position = "dodge", alpha = 0.8, bins = 10) +
  xlim(c(0, 0.5)) +
  scale_colour_manual(values = c("red", "darkgreen", "darkblue")) +
  scale_fill_manual(values = c("red", "darkgreen", "darkblue")) +
  labs(y = "ASoC SNP count (FDR < 0.05)") +
  theme_classic()

# save plot information 
h_plot <- ggplot(df_minRatio, aes(x = minRatio, 
                                  colour = index, 
                                  fill = index)) +
  geom_histogram(position = "dodge", alpha = 0.8, bins = 11) +
  xlim(c(0, 0.5)) +
  scale_colour_manual(values = c("red", "darkgreen", "darkblue")) +
  scale_fill_manual(values = c("red", "darkgreen", "darkblue")) +
  labs(y = "ASoC SNP count (FDR < 0.05)") +
  theme_classic()

## get plotting information as data.frame
h_plotdata <- ggplot_build(h_plot)$data[[1]]
h_plotdata$group <- as.factor(h_plotdata$group)
levels(h_plotdata$group) <- c("> 20", "> 100", "> 200")
h_plotdata$index_x <- 1
h_plotdata$index_x <- h_plotdata$x[12:22]
h_plotdata$y[1:11] <- h_plotdata$y[1:11]/5612
h_plotdata$y[12:22] <- h_plotdata$y[12:22]/3667
h_plotdata$y[23:33] <- h_plotdata$y[23:33]/2249


## plot with geom_bar
ggplot(h_plotdata, aes(x = factor(index_x), 
                       y = y, fill = group)) +
  geom_bar(stat = "identity", 
           position = "dodge") +
  labs(x = "minRatio", y = "ASoC SNP count (FDR < 0.05)") +
  scale_colour_manual(values = c("red", "darkgreen", "darkblue")) +
  scale_fill_manual(values = c("red", "darkgreen", "darkblue")) +
  scale_y_continuous(labels = scales::percent) +
  # scale_x_continuous(labels = scales::comma, 
  #                    breaks = seq(0, 0.5, by = 0.1)) +
  theme_classic() #+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5))


ggplot(df_bar_ggplot, aes(SNP_index, colour = SNP_colour)) +
  geom_bar(position = "identity", )
#####

df_more_than_200 <- data.frame(x = (1:5) * 0.1, 
                               y = 0)
i <- 1
for (i in 1:5) {
  print(i)
  df_more_than_200$y[i] <- sum(sig_ASoC_temp$minRatio <= (i * 0.1)) / nrow(sig_ASoC_temp)
}

ggplot(df_more_than_200, aes(x = x, y = y)) +
  geom_bar(stat = "identity", 
           colour = "darkblue", 
           fill = "darkblue") +
  scale_x_reverse() +
  scale_y_continuous(labels = scales::percent, 
                     breaks = seq(0, 1, 0.1)) +
  labs(x = "Minimum SNP allelic ratio of ATAC-seq reads", 
       y = "Percentage of total ASoC SNPs (FDR < 0.05)") +
  theme_classic() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14))

#####
i <- 1
read_depth <- 200
ASoC_temp <- ASoC_w_pBinom[ASoC_w_pBinom$DP >= read_depth, ]
sig_ASoC_temp <- sig_ASoC_w_pBinom[sig_ASoC_w_pBinom$DP >= read_depth, ]
# ASoC_percentage_at_ratio <- 1:25
df_temp <- data.frame(x = 1:50, 
                      y = 1:50, 
                      colour = read_depth)

for (i in 1:50) {
  
  df_temp$x[i] <- 0.01 * i
  df_temp$y[i] <- sum(sig_ASoC_temp$minRatio < (0.01 * i)) /
    sum(ASoC_temp$minRatio < (0.01 * i))
  print(paste(i, df_temp$y[i], sep = ", "))
}
# df_temp <- data.frame(x = 0.02 * ASoC_percentage_at_ratio, 
#                       y = ASoC_percentage_at_ratio, 
#                       colour = read_depth)
# df_plot_ratio_vs_sig_perc <- df_temp
df_plot_ratio_vs_sig_perc <- as.data.frame(rbind(df_plot_ratio_vs_sig_perc,
                                                 df_temp))

# rm(df_plot_ratio_vs_sig_perc)

ggplot(df_plot_ratio_vs_sig_perc, aes(x = x, y = y, 
                                      colour = factor(colour))) +
  geom_point(alpha = 0.8) +
  scale_y_continuous(labels = scales::percent, 
                     breaks = seq(0, 1, 0.1)) +
  labs(x = "Minimum SNP allelic ratio of ATAC-seq reads", 
       y = "Detectable  ASoC SNPs (FDR < 0.05)  \nout of all heterozygous SNPs", 
       colour = "Read Depth") +
  scale_colour_manual(values = c("red", "darkgreen", "darkblue")) +
  scale_x_reverse() +
  theme_classic() +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14))




#######
ggplot(ASoC_w_pBinom_plot, aes(x = REF_ratio, y = DP)) +
  geom_point(aes(colour = `-log10P`), 
             cex = 0.2, alpha = 1) +
  # theme_classic() +
  # ylim(0, 20) +
  labs(x = "reference allele ratio", 
       y = "Read depth") +
  # scale_color_gradientn(colours = c("lightblue", 
  #                                   "red", "red", "red", 
  #                                   "red", "red", "red", "darkred")) +
  scale_color_gradientn(colours = c("lightskyblue", "orange", "red", "red"),
                        # values = c(0, 1, 5, 10, 20, 50), 
                        limits = c(0, 20),
                        na.value = "red") +
  # scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11, "matlab.like")))(100), 
  #                        limits = c(0, 10)) +
  # scale_x_continuous(trans='log2') +
  # scale_y_continuous(trans='log2') +
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8)) +
  # xlim(0, 2) +
  # ylim(0, 2000) +
  ggtitle(label = "Use read depth at SNP site +/- 25 bp") +
  theme_classic()
