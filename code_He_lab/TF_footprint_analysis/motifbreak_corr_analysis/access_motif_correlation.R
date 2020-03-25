# Following is the code to correlate the previously obtained
# SNP disruption effects on TF motifs ('SNP_motif_disruption.R')
# with the allelic imbalance of accessibilty in given cell type,
# and to generate Figures 2E,2F,S15 --
# "Correlation analysis of predicted TF motif disruption scores with the allelic imbalance of chromatin accessibility"


library(data.table)
library(dplyr)
library(motifbreakR)
library(ggplot2)
library(gridExtra)

options(stringsAsFactors = FALSE)
SNP_dir = '~/Downloads/ASoC/data/ASoC_per_celltype/'
footprint_dir = '~/Downloads/ASoC/TF_motifs/footprint_jaspar_hsap_redundant/'

# 1. Correlation analysis ####
merge_result = function(SNP_dir,celltype,mb_res_suffix){
  ## Merge SNP motif disruption result with allelic imbalance stats
  SNP.stats = data.frame(fread(file = paste0(SNP_dir,celltype,'_w_FDR.txt'),
                               sep = '\t',header = T,na.strings = "NA"))
  names(SNP.stats)[1:2] = c('chr','end')
  SNP.stats = mutate(SNP.stats,start = end-1) %>% mutate(zscore = (ALT-DP/2)/sqrt(DP*0.25)) %>% 
    mutate(log_alt_ref = log2(ALT/REF))
  mb.result = readRDS(paste0('data/',celltype,mb_res_suffix))
  mb.result.df = data.frame(mb.result@ranges)
  tmp.df = data.frame(mcols(mb.result)[,-c(6,7,14,15)])
  mb.result.df = cbind(mb.result.df,tmp.df)
  mb.result.df = mutate(mb.result.df, Dpct = pctAlt-pctRef) %>% mutate(Dscore = scoreAlt-scoreRef)
  mb.result.df$rsID = sapply(stringr::str_split(mb.result.df$names,':'), function(x){x[1]})
  
  mb.result.df = left_join(mb.result.df, select(SNP.stats,c('rsID','zscore','log_alt_ref')), by='rsID')
  mb.result.df = filter(mb.result.df,effect=='strong')
  return(mb.result.df)
}

regress_summary = function(mb.result.df, motif_cluster, num_thres = 5){
  ## Linear regression of the two stats across SNPs in the same TF motif
  mb.result.df = group_by(mb.result.df, providerId) %>% 
    filter(n()>num_thres) %>% ungroup() # Filter motif groups w/ too few SNPs
  uniq_motifs = unique(mb.result.df$providerId)
  len = length(uniq_motifs)
  motif.stats = data.frame(motif = uniq_motifs, beta1 = rep(NA,len), 
                           r2 = rep(NA,len), pval = rep(NA,len))
  for (i in 1:len){
    mb.result.tmp = filter(mb.result.df, providerId == uniq_motifs[i])
    lm.summary = summary(lm(log_alt_ref ~ Dpct - 1, data=mb.result.tmp))
    if (lm.summary$r.squared==0){
      next
    } else {
      motif.stats$beta1[i] = lm.summary$coefficients[1]
      motif.stats$r2[i] = lm.summary$r.squared
      motif.stats$pval[i] = lm.summary$coefficients[4]
    }
  }
  motif.stats = na.omit(motif.stats)
  motif.stats$fdr = p.adjust(motif.stats$pval, method = 'BH')
  motif.stats = inner_join(motif.stats, motif_cluster, by = 'motif')
  signif_motif.stats = filter(motif.stats, fdr<0.05) %>% group_by(cluster_id) %>% arrange(fdr,.by_group = T)
  return(list(all = motif.stats, signif = signif_motif.stats))
}

# Load motif cluster info (mainly based on the clustering trees in the JASPAR database)
motif_cluster = read.delim('~/Downloads/ASoC/TF_motifs/jaspar2018_hsapiens_motif_cluster_imputed.txt',
                           sep = '\t', header = T, na.strings = 'NA')
motif_cluster = na.omit(motif_cluster)

## Correlation results:

for (celltype in c('CN20','NSC20')) {
  mb.result = merge_result(SNP_dir,celltype,'_in_called_fp.motifBreak.rds')
  motif.stats = regress_summary(mb.result, motif_cluster)
  saveRDS(motif.stats,paste0('data/',celltype,'_in_called_fp.corr_summary.rds'))
}

# 2. Heatmap comparison of correlation results across cell types (Fig S15A, Fig 2E) ####
## Plot correlation heatmap for all significant TF motifs in either cell types (Fig S15A)
joint_heatmap = function(corr_res_suffix, fdr_thres = 0.05){
  celltype = c('CN20','NSC20')
  motif.stats = list()
  for (i in 1:length(celltype)){
    motif.stats.tmp = readRDS(paste0('data/',celltype[i],corr_res_suffix))
    motif.stats.tmp = motif.stats.tmp$all
    motif.stats.tmp = mutate(motif.stats.tmp, zscore = ifelse(beta1>0, abs(qnorm(pval)), -abs(qnorm(pval))) )
    motif.stats[[i]] = motif.stats.tmp
  }
  
  signif_motif = unique(filter(plyr::ldply(motif.stats), fdr<fdr_thres)$motif) # 48 unique motifs (inFP log(alt/ref))
  template.df = plyr::join_all(motif.stats, by = c("motif","TF","cluster_id","sequence"), type = 'full') %>% 
    filter(motif %in% signif_motif) %>% select(c("motif","TF","cluster_id","sequence"))
  
  joined.lst = list()
  i = 1
  tmp = right_join(motif.stats[[i]],template.df,by=c("motif","TF","cluster_id","sequence"))
  tmp$celltype = rep(celltype[i],nrow(template.df))
  tmp = group_by(tmp, cluster_id) %>% arrange(zscore, .by_group = T)
  joined.lst[[i]] = tmp
  motif_ordered = tmp$motif
  
  for (i in 2:length(celltype)){
    tmp = right_join(motif.stats[[i]],template.df,by=c("motif","TF","cluster_id","sequence"))
    tmp$celltype = rep(celltype[i],nrow(template.df))
    tmp = tmp[match(motif_ordered,tmp$motif),]
    joined.lst[[i]] = tmp
  }
  
  joined.df = plyr::ldply(joined.lst) %>% select(c('celltype','zscore','pval','cluster_id','motif','TF'))
  joined.df$cluster_id = factor(joined.df$cluster_id, levels = unique(joined.df$cluster_id))
  rep_times = as.numeric(table(joined.df$cluster_id[1:length(signif_motif)]))
  
  joined.df$bracket = rep(c(1:length(rep_times)), times = rep_times)
  str.tmp = paste(joined.df$bracket,collapse = '],[')
  str.tmp = paste0('[',str.tmp,']')
  str.tmp = stringr::str_split(str.tmp,',')[[1]]
  joined.df$bracket = str.tmp
  
  joined.df = mutate(joined.df, cluster_TF = paste(motif,TF,bracket,sep = ' ')) 
  joined.df$cluster_TF = factor(joined.df$cluster_TF, levels = joined.df$cluster_TF[1:length(signif_motif)])
  joined.df$celltype = factor(joined.df$celltype,levels = celltype)
  joined.df$celltype = plyr::revalue(joined.df$celltype, c("CN20"="iN-Glut 20", "NSC20"="NPC 20"))
  joined.df = mutate(joined.df, polar_pval = -log10(pval)*sign(zscore))
  
  plot_out = ggplot(data = joined.df, aes(x=celltype, y=cluster_TF)) + geom_tile(aes(fill=polar_pval)) + 
    scale_fill_gradient2(high = 'darkred', mid = 'white', low = "darkblue",
                         na.value = 'gray', name = expression('polarized '*-log[10]*' P value')) + 
    theme_minimal() + xlab('Cell Type') + ylab('Motif') + 
    theme(axis.title = element_text(size = 14,face="bold"), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 11),
          legend.title = element_text(size = 13,face="bold"), legend.text = element_text(size = 11))
  return(list(data=joined.df,plot=plot_out))
}

joint_res = joint_heatmap('_in_called_fp.corr_summary.rds')
joined.df = joint_res$data
saveRDS(joined.df,'data/Glut_NPC_joined_heatmap_table.rds')
ggsave('Figs/heatmap_CN20_vs_NSC20_inFP.pdf',joint_res$plot,width = 8,height = 16)

## Plot heatmap of selected TFs (Fig 2E)
TF_string = 'FOXH1, FOXF2, SRY, CTCF, NHLH1, TCF3, TCF4, PLAG1, MAC, CLOCK, RORA, NR2F1, SP1, KLF16, EGR4, MEIS2, SOX9, SOX15, FOSB:JUN, ELK1, TBX5, BACH2'
TF_split = stringr::str_split(TF_string,', ')[[1]]
selected_joined.df = filter(joined.df, TF %in% TF_split)
selected_joined.df = selected_joined.df[-c(14,31),] # Remove duplicated TFs

ggplot(data = selected_joined.df, aes(x=celltype, y=cluster_TF)) + geom_tile(aes(fill = polar_pval)) + 
  scale_fill_gradient2(high = 'darkred', mid = 'white', low = "darkblue",
                       na.value = 'gray', name = expression('polarized '*-log[10]*' P value')) + 
  theme_minimal() + xlab('Cell Type') + ylab('Motif') + 
  theme(axis.title = element_text(size = 14,face="bold"), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 13,face="bold"), legend.text = element_text(size = 11))
ggsave('Figs/heatmap_CN20_vs_NSC20_inFP.selected.pdf',last_plot(),width = 8,height = 10)

# 3. Plot individual TF motif cases in Glut-20 (Fig 2F) ####
plot_motif = function(plot.df, title_text, x.coords = c(-0.17,-0.155), y.coords = c(3,2.5)){
  lm.summary = summary(lm(log_alt_ref ~ Dpct - 1, data = plot.df))
  r2 = lm.summary$r.squared
  pval = lm.summary$coefficients[4]
  plot_out = ggplot(plot.df, aes(x=Dpct, y=log_alt_ref)) + geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype="dashed") +
    geom_abline(slope = lm.summary$coefficients[1],intercept = 0, 
                color='red', size=0.8) +
    scale_x_continuous(name = 'motif score change', limits = c(-0.205,0.205), breaks=seq(-0.2,0.2,0.1)) +
    scale_y_continuous(name = expression('log'[2]*' (alt count / ref count)'), limits = c(-3,3), breaks=seq(-3,3,1)) +
    labs(title = title_text) +
    theme_classic() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=18),
          legend.title = element_text(size=15),
          legend.text = element_text(size=14)) +
    annotate('text', x = x.coords, y = y.coords, 
             label=c(paste0('R^2 == ',signif(r2,2)), paste0('P == ',signif(pval,2))),
             parse = T,size = 5)
  return(plot_out)
}

mb.result = merge_result(SNP_dir,'CN20','_in_called_fp.motifBreak.rds') # result from CN20 ASoC SNPs in FP
motif.stats = regress_summary(mb.result, motif_cluster)
## Scatter plot of correlation for SNPs in SOX9 motifs (Fig 2F)
mb.result.SOX9 = filter(mb.result, providerId=='MA0077.1')
plot_out = plot_motif(mb.result.SOX9,'SOX9 Correlation - Glut20 ASoC SNPs in Footprints')
ggsave('Figs/SOX9_scatterplot_CN20_ASoC_inFP.pdf',plot_out,width = 5,height = 6)

## Merged scatter plot for selected TF motifs
selected_motifs = data.frame(motif=c('MA0077.1','MA0017.2','MA0048.2','MA0479.1','MA0830.1'),
                             TF_name=c('SOX9','NR2F1','NHLH1','FOXH1','TCF4'))
selected_motifs$pval = motif.stats$all$pval[match(selected_motifs$motif,motif.stats$all$motif)]
selected_motifs$pval = formatC(selected_motifs$pval, format = "e", digits = 1) # Scientific Notation
selected_motifs = mutate(selected_motifs, TF=paste0(TF_name,' (',pval,')'))

selected_mb.result = filter(mb.result, providerId %in% selected_motifs$motif) %>% 
  select(c('rsID','Dpct','log_alt_ref','providerId'))
selected_mb.result$TF = selected_motifs$TF[match(selected_mb.result$providerId,selected_motifs$motif)]
selected_mb.result$providerId = factor(selected_mb.result$providerId, levels = selected_motifs$motif)
selected_mb.result$TF = factor(selected_mb.result$TF, levels = selected_motifs$TF)

ggplot(selected_mb.result,aes(x=Dpct, y=log_alt_ref, color=TF, shape=TF)) + 
  geom_point(size = 2) + 
  scale_color_manual(values=c("#E64B35B2","#3C5488B2","#7CAE00","#4DBBD5B2","#F39B7FB2")) +
  geom_smooth(method = 'lm', fill = NA,formula = 'y~x-1', fullrange = T) +
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(name = 'motif score change', limits = c(-0.25,0.25), expand = c(0,0))+ #, breaks=seq(-0.2,0.2,0.1)) +
  scale_y_continuous(name = expression('log'[2]*' (alt count / ref count)'), limits = c(-3,3), expand = c(0,0)) + #, breaks=seq(-3,3,1)) +
  theme_bw() + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18),
        legend.title = element_blank(),legend.text = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave('Figs/Merged_scatterplot_CN20_ASoC_inFP.pdf', last_plot(), width = 8, height = 5)

# 4. Loop through ALL significant TF motifs in CN20 and generate scatter plots (Fig S15B)  ####
# Obtain CN20 'mb.result', 'motif.stats' the same way as in Part 3.
joined.df = readRDS('data/Glut_NPC_joined_heatmap_table.rds')
signif_motif.df = inner_join(joined.df[1:48,c("motif","TF","cluster_TF","polar_pval")],
                             motif.stats$signif, by=c("motif","TF"))

signif_motif.df$cluster_TF = as.character(signif_motif.df$cluster_TF)
signif_motif.df$cluster_TF[c(2,3,4,14)] = c('MA1139.1 FOSL2::JUNB [4]','MA1131.1 FOSL2::JUN [4]',
                                            'MA1136.1 FOSB::JUNB [4]','MA0828.1 SREBF2 [11]')
signif_mb.result.df = filter(mb.result, providerId %in% signif_motif.df$motif)

plot.lst = list()
for (i in 1:length(signif_motif.df$motif)){
  mb.result.tmp = filter(mb.result, providerId == signif_motif.df$motif[i])
  lm.summary = summary(lm(log_alt_ref ~ Dpct - 1, data=mb.result.tmp))
  beta1 = lm.summary$coefficients[1]
  r2 = lm.summary$r.squared
  pval = lm.summary$coefficients[4]
  ymax = max(mb.result.tmp$log_alt_ref)
  ymin = min(mb.result.tmp$log_alt_ref)
  gap = (ymax-ymin)/20
  outplot = ggplot(mb.result.tmp, aes(x=Dpct, y=log_alt_ref)) + geom_point(size=1) +
    geom_abline(slope = beta1, intercept = 0, linetype="dashed",
                color = ifelse(beta1>0,'red','blue')) +
    scale_x_continuous(breaks = seq(-0.2,0.2,0.1)) +
    # labs(x = 'motif score change', y = expression('log'[2]*' (alt count / ref count)'), 
    labs(title = signif_motif.df$cluster_TF[i]) +
    theme_bw() + theme(axis.text=element_text(size=9),axis.title=element_blank(),
                       plot.title = element_text(size = 9,face = 'bold',hjust = 0.5),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    annotate('text', x = c(0,0), y = c(ymax-gap,ymax-3*gap), # y = c(0.22,-0.22), 
             label=c(paste0('R^2 == ',signif(r2,3)), paste0('pval == ',signif(pval,3))),
             color = ifelse(beta1>0,'red','blue'), parse = T,size = 3)
  plot.lst[[i]] = outplot
}
args = c(plot.lst, list(nrow = 8, ncol = 4, bottom = "motif score change",
                        left = 'allelic imbalance of chromatin accessibility'))
pdf('Figs/signif_TF_regression_grid_plots.pdf',width = 8,height = 16)
do.call(grid.arrange,args)
dev.off()
