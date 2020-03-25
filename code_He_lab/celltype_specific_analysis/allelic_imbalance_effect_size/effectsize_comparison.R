## Correlation between the allelic ratios at ASoC SNPs 
## in the ascertained neuronal cell type vs in iPSC
## Fig S18 right panels


library(dplyr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
 
get_cell_spec_FC = function(A,A.formal_name,REF,REF.formal_name){
  REF.full = readRDS(paste0('../SNP_accessibilty/',REF,'_all_SNPs.rds'))
  REF.full = mutate(REF.full,log_ratio=log(REF_C/ALT_C))
  A.full = readRDS(paste0('../SNP_accessibilty/',A,'_all_SNPs.rds'))
  A.full = mutate(A.full,log_ratio=log(REF_C/ALT_C))
  A.sig_in_REF = inner_join(filter(A.full, FDR<0.05), REF.full, by='ID')
  ## Exclude SNPs that are also significant in REF cell type:
  A.spec_in_REF = filter(A.sig_in_REF,pBinom.y>0.05)
  
  ## Linear regression
  lm_sum = summary(lm(log_ratio.y ~ log_ratio.x, data = A.spec_in_REF))
  r2 = lm_sum$r.squared
  pval = lm_sum$coefficients[2,4]
  slope = lm_sum$coefficients[2,1]
  intercept = lm_sum$coefficients[1,1]
  
  plot_out = ggplot(data = A.spec_in_REF, aes(x = log_ratio.x, y = log_ratio.y)) +
    geom_point(aes(color = log2(DP.y)), size = 1) +
    scale_color_gradientn(colours = c("lightgrey", "blue", "darkblue"), 
                          name = bquote(~log[2]~'(read depth) in '~.(REF.formal_name))) +
                          # name = paste('log2 (read depth) in',REF.formal_name)) +
    geom_abline(slope = slope, intercept = intercept, 
                color = "red", linetype = 'dashed', size = 1.2) +
    geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    xlim(c(-2.3, 2.3)) + ylim(c(-1,1)) +
    xlab(paste("log(REF/ALT) of ASoC SNPs in",A.formal_name)) +
    ylab(paste("log(REF/ALT) in",REF.formal_name)) + 
    ggtitle(label = paste(nrow(A.spec_in_REF), A.formal_name,'specific ASoC SNPs vs',REF.formal_name)) +
    theme_minimal() + theme(plot.title = element_text(size = 14), 
          axis.title = element_text(size = 16), axis.text = element_text(size = 16), 
          legend.title = element_text(size = 12), legend.text = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(list(plot = plot_out, df = A.spec_in_REF, r2 = r2, pval = pval))
}
 
ascertain_FC = function(A,REF){
  REF.full = readRDS(paste0('../SNP_accessibilty/',REF,'_all_SNPs.rds'))
  REF.full = mutate(REF.full,log_ratio=log(REF_C/ALT_C))
  A.full = readRDS(paste0('../SNP_accessibilty/',A,'_all_SNPs.rds'))
  A.sig = filter(A.full, FDR<0.05)$ID
  A.sig_in_REF = A.sig[A.sig %in% REF.full$ID]
  A.sig_fullinfo = A.full[A.full$ID %in% A.sig_in_REF, ] %>% mutate(log_ratio=log(REF_C/ALT_C))
  sig_A_vs_ips = merge(A.sig_fullinfo, REF.full, by = "ID")
  
  ## color with lfdr level
  sig_A_vs_ips$lfdr = qvalue::lfdr(sig_A_vs_ips$pBinom.y)
  sig_A_vs_ips$readdepth_100 = ifelse(sig_A_vs_ips$DP.y > 100, '>100','<=100')
  
  ## Linear regression
    lm_sum.whole = summary(lm(log_ratio.y ~ log_ratio.x, data = sig_A_vs_ips))
    r2_whole = lm_sum.whole$r.squared
    pval_whole = lm_sum.whole$coefficients[2,4]
    slope_whole = lm_sum.whole$coefficients[2,1]
    intercept_whole = lm_sum.whole$coefficients[1,1]
  
  lm_sum.signif = summary(lm(log_ratio.y ~ log_ratio.x, data = filter(sig_A_vs_ips, lfdr<=0.1)))
  r2_signif = lm_sum.signif$r.squared
  pval_signif = lm_sum.signif$coefficients[2,4]
  slope_signif = lm_sum.signif$coefficients[2,1]
  intercept_signif = lm_sum.signif$coefficients[1,1]
  
  return(list(df = sig_A_vs_ips, slope_whole=slope_whole, intercept_whole=intercept_whole,
              slope_signif=slope_signif, intercept_signif=intercept_signif,
              r2_whole = r2_whole, pval_whole = pval_whole,
              r2_signif = r2_signif, pval_signif = pval_signif))
}

scatter_plot = function(A.res,A.formal_name,REF.formal_name){
  sig_A_vs_ips = A.res$df
  xmax = ceiling(max(sig_A_vs_ips$log_ratio.x)*10)/10
  xmin = floor(min(sig_A_vs_ips$log_ratio.x)*10)/10
  ymax = ceiling(max(sig_A_vs_ips$log_ratio.y)*10)/10
  ymin = floor(min(sig_A_vs_ips$log_ratio.y)*10)/10
  
  plot_out = ggplot(data = sig_A_vs_ips, aes(x = log_ratio.x, y = log_ratio.y)) +
    geom_point(aes(fill = -log10(lfdr), shape = factor(readdepth_100)), size = 1, col = "#888888") +
    scale_fill_gradientn(colours = c("white", "blue", "darkblue"), name = paste('-log10 (lfdr) in',REF.formal_name)) +
    scale_shape_manual(values=c(21, 24), name = paste('Read depth in',REF.formal_name))+
    theme_minimal() +
    geom_abline(slope = A.res$slope_signif, intercept = A.res$intercept_signif, color = "magenta", linetype = 'dashed', size = 1.2) +
    geom_abline(slope = A.res$slope_whole, intercept = A.res$intercept_whole, color = "royalblue", size = 1.2) +
    xlab(paste("log(REF/ALT) of ASoC SNPs in",A.formal_name)) +
    ylab(paste("log(REF/ALT) in",REF.formal_name)) + 
    xlim(c(xmin, xmax)) + ylim(c(ymin, ymax)) +
    geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    ggtitle(label = paste(nrow(sig_A_vs_ips), A.formal_name,'ASoC SNPs w/ DP >= 20 in',REF.formal_name)) +
    theme(plot.title = element_text(size = 14), 
          axis.title = element_text(size = 16), axis.text = element_text(size = 16), 
          legend.title = element_text(size = 14), legend.text = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(plot_out)
}

# 1. Plot effect sizes for all ASoCs ascertained in the given neuronal cell type (vs iPS) ####

## CN20 vs iPS
A.res = ascertain_FC('CN20','iPS_8')
plot_out = scatter_plot(A.res,'Glut-20','iPSC-8')
plot_out = plot_out + 
  geom_text(x = -1.5, y = 2, colour = "royalblue", size = 5,
            label = as.expression(italic(r)^2~"="~0.32*"," ~ italic(P)~"="~1.2 %.% 10^-184), 
            parse = TRUE) +
  geom_text(x = -1.5, y = 1.8, colour = "magenta", size = 5,
            label = as.expression(italic(r)^2~"="~0.79*"," ~ italic(P)~"="~4.1 %.% 10^-74), 
            parse = TRUE)
ggsave(paste0('Figs/effectsize_all_ASoCs_CN20_vs_iPS.pdf'), plot_out, width = 8, height = 6)

## CN-8 vs iPS
A.res = ascertain_FC(A = 'CN_8', REF = 'iPS_8')
plot_out = scatter_plot(A.res,'Glut-8','iPSC-8')
plot_out = plot_out + 
  geom_text(x = -1.5, y = 2, colour = "royalblue", size = 5,
            label = as.expression(italic(r)^2~"="~0.31*"," ~ italic(P)~"="~9.8 %.% 10^-89), 
            parse = TRUE) +
  geom_text(x = -1.5, y = 1.8, colour = "magenta", size = 5,
            label = as.expression(italic(r)^2~"="~0.77*"," ~ italic(P)~"="~5.3 %.% 10^-50), 
            parse = TRUE)
ggsave(paste0('Figs/effectsize_all_ASoCs_CN8_vs_iPS.pdf'), plot_out, width = 8, height = 6)

## GA-8 vs iPS
A.res = ascertain_FC(A = 'GA_8', REF = 'iPS_8')
plot_out = scatter_plot(A.res,'GA-8','iPSC-8')
plot_out = plot_out + 
  geom_text(x = -1.4, y = 2, colour = "royalblue", size = 5,
            label = as.expression(italic(r)^2~"="~0.47*"," ~ italic(P)~"="~8.5 %.% 10^-106), 
            parse = TRUE) +
  geom_text(x = -1.4, y = 1.8, colour = "magenta", size = 5,
            label = as.expression(italic(r)^2~"="~0.85*"," ~ italic(P)~"="~2.0 %.% 10^-67), 
            parse = TRUE)
ggsave(paste0('Figs/effectsize_all_ASoCs_GA8_vs_iPS.pdf'), plot_out, width = 8, height = 6)

## DN-8 vs iPS
A.res = ascertain_FC(A = 'DN_8', REF = 'iPS_8')
plot_out = scatter_plot(A.res,'DN-8','iPSC-8')
plot_out = plot_out + 
  geom_text(x = -1.5, y = 2, colour = "royalblue", size = 5,
            label = as.expression(italic(r)^2~"="~0.39*"," ~ italic(P)~"="~2.1 %.% 10^-54), 
            parse = TRUE) +
  geom_text(x = -1.5, y = 1.8, colour = "magenta", size = 5,
            label = as.expression(italic(r)^2~"="~0.79*"," ~ italic(P)~"="~8.8 %.% 10^-39), 
            parse = TRUE)
ggsave(paste0('Figs/effectsize_all_ASoCs_DN8_vs_iPS.pdf'), plot_out, width = 8, height = 6)

## NSC-8 vs iPS
A.res = ascertain_FC(A = 'NSC_8',REF = 'iPS_8')
plot_out = scatter_plot(A.res,'NPC-8','iPSC-8')
plot_out = plot_out + 
  geom_text(x = -1.4, y = 2, colour = "royalblue", size = 5,
            label = as.expression(italic(r)^2~"="~0.43*"," ~ italic(P)~"="~8.8 %.% 10^-48), 
            parse = TRUE) +
  geom_text(x = -1.4, y = 1.8, colour = "magenta", size = 5,
            label = as.expression(italic(r)^2~"="~0.89*"," ~ italic(P)~"="~1.6 %.% 10^-35), 
            parse = TRUE)
ggsave(paste0('Figs/effectsize_all_ASoCs_NSC8_vs_iPS.pdf'), plot_out, width = 8, height = 6)

# 2. Plot effect sizes only for cell-type-specific ASoCs (w/ pval>0.05 in the other) ####
## CN-8 vs iPS
CN8_vs_iPS = get_cell_spec_FC(A = 'CN_8',A.formal_name = 'Glut-8',REF = 'iPS_8',REF.formal_name = 'iPSC-8')
final_plot =  CN8_vs_iPS$plot
final_plot = final_plot + geom_abline(slope = 1,intercept = 0,color='grey46') +
  geom_text(x = -2, y = 0.9, colour = "red", size = 5, label = as.expression(italic(r)^2~"="~0.16),parse = T) +
  geom_text(x = 1.3, y = 0.9, colour = "grey46", size = 5, label = '1:1 line') 
ggsave(paste0('Figs/effectsize_CN8_spec_vs_iPS.pdf'), final_plot, width = 9, height = 4)

## Plot GA-8 vs NSC-8
GA8_vs_NSC = get_cell_spec_FC(A = 'GA_8',A.formal_name = 'GA-8',REF = 'NSC_8',REF.formal_name = 'NPC-8')
final_plot =  GA8_vs_NSC$plot 
final_plot = final_plot +  geom_abline(slope = 1,intercept = 0,color='grey46') +
  geom_text(x = -2, y = 0.9, colour = "red", size = 5,label = as.expression(italic(r)^2~"="~0.38),parse = T) +
  geom_text(x = 1.3, y = 0.9, colour = "grey46", size = 5, label = '1:1 line')
ggsave(paste0('Figs/effectsize_GA8_spec_vs_NSC8.pdf'), final_plot, width = 9, height = 4)

# 3. Count # of data points in quadrants II and IV of cell-type-specific scatter plots ####
count_quadrant24 = function(df){
  tmp2.1 = nrow(filter(df,log_ratio.x<0) %>% filter(log_ratio.y>0))
  tmp2.2 = nrow(filter(df,log_ratio.x<0) %>% filter(log_ratio.y>=0))
  num_in_Q2 = (tmp2.1+tmp2.2)/2
  tmp4.1 = nrow(filter(df,log_ratio.x>0) %>% filter(log_ratio.y<0))
  tmp4.2 = nrow(filter(df,log_ratio.x>0) %>% filter(log_ratio.y<=0))
  num_in_Q4 = (tmp4.1+tmp4.2)/2
  precent_24 = (num_in_Q2+num_in_Q4)/nrow(df)
  return(precent_24)
}

count_quadrant24(CN8_vs_iPS$df) # 33.5%
count_quadrant24(GA8_vs_NSC$df) # 21.2%

# 4. Distribution of slopes of data points in cell-type-specific scatter plots ####
## CN-8 vs iPS
df = CN8_vs_iPS$df
sum(df$log_ratio.y/df$log_ratio.x < 0.5)/nrow(df) # 84.6%
ggplot(df,aes(x=log_ratio.y/log_ratio.x)) + geom_histogram(bins=50) +
  geom_vline(xintercept = 0.5, color='red') + 
  labs(x='Y.log-ratio/X.log-ratio (iN-Glut>iPSC)', y='counts') +
  scale_x_continuous(breaks=seq(-2,2,0.5)) +
  theme_classic() + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16)) +
  geom_text(x = -0.8, y = 55, colour = "red", size = 5,
            label = '84.6% < 0.5')
ggsave('Figs/slope_hist_CN8_vs_iPS.pdf', last_plot(), width = 6,height = 4)

## GA-8 vs NSC-8
df = GA8_vs_NSC$df
sum(df$log_ratio.y/df$log_ratio.x < 0.5)/nrow(df) # 68.3%
ggplot(df,aes(x=log_ratio.y/log_ratio.x)) + geom_histogram(bins=50) +
  geom_vline(xintercept = 0.5,color='red') + 
  labs(x='Y.log-ratio/X.log-ratio (iN-GA>NPC)',y='counts') +
  scale_x_continuous(breaks=seq(-2,2,0.5)) +
  theme_classic() + 
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16)) +
  geom_text(x = -0.7, y = 30, colour = "red", size = 5,
            label = '68.3% < 0.5')
ggsave('Figs/slope_hist_GA8_vs_NSC8.pdf', last_plot(), width = 6,height = 4)
