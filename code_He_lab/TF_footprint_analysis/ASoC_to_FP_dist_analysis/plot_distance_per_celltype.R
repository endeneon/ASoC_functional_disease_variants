# Following is the code to analyze the distance of ASoC SNPs 
# from identitified TF footprints per cell type,
# and to generate Figures S14B-G --
# "TF footprint analysis using ATAC-seq data in each cell type"


library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = F)
snp_dir = '~/Downloads/ASoC/data/ASoC_per_celltype/'
footprint_dir = '~/Downloads/ASoC/TF_motifs/footprint_jaspar_hsap_redundant/'

## Distribution of SNP distance from its nearest TF footprint (5 core 8 cell types) ####
# Figures S14C-G
make.ranges = function(df, metadata){
  dfRanges = GRanges(seqnames = df$chr, 
                      ranges = IRanges(start = df$start,
                                       end = df$end))
  for(m in metadata){
    dfRanges = plyranges::mutate(dfRanges, !!m := unname(df[,m]))
  }
  return(dfRanges)
}
dist_plot = function(stats.df, title_text){
  plot1 = ggplot(stats.df, aes(log10(dist+1))) + 
    geom_histogram(bins = 30) + 
    labs(x = expression('log'[10]*'(distance of SNP from nearest TF footprint + 1)'),
         y = expression('count'),title = title_text) + 
    theme_classic() +
    theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
          plot.title = element_text(size=16,face = 'bold'))
  print(plot1)
}

covert_tb = data.frame(orig=c('NSC','CN','DN','GA','IPS'),
                       new=c('NPC Core 8','iN-Glut Core 8','iN-DN Core 8','iN-GA Core 8','iPS Core 8'))
for (celltype in covert_tb$orig){
  title_text = covert_tb$new[covert_tb$orig==celltype]
  ASoC_snp = as_tibble(fread(file = paste0(snp_dir,celltype,'_core8_w_FDR0.05.bed'), 
                             sep = '\t',header = F,na.strings = "NA",
                             col.names = c('chr','start','end','rsID','FDR','dp','ref','alt')))
  snpRanges = make.ranges(ASoC_snp, c('rsID','FDR'))
  
  merged_footprint = as_tibble(fread(file = paste0(footprint_dir,celltype,'_merged_footprint.jaspar2018_hsapiens.dnase2tf_fdr0.05.bed'), 
                                     sep = '\t', header = F,
                                     col.names = c('chr','start','end','TF','score','strand')))
  tfRanges = make.ranges(merged_footprint, c('TF'))
  
  dist = distanceToNearest(snpRanges,tfRanges)
  ASoC_snp$dist = dist@elementMetadata$distance
  
  plot_out = dist_plot(ASoC_snp,title_text)
  ggsave(paste0('Figs/dist_',celltype,'8_ASoC_fp.pdf'),plot = plot_out, width = 6,height = 5)
}

## Enrichment in footprint (ASoC vs non-ASoC) in iN-Glut 20 ####
# Figure S14B
fisher_enrich = function(snp.df, dist_cutoff = 0){
  ASoC_snp = filter(snp.df,asoc==1)
  nonASoC_snp = filter(snp.df,asoc==0)
  
  f1 = sum(ASoC_snp$dist <= dist_cutoff)
  f3 = dim(ASoC_snp)[1]-f1
  f2 = sum(nonASoC_snp$dist <= dist_cutoff)
  f4 = dim(nonASoC_snp)[1]-f2
  fisher_tb = matrix(c(f1,f2,f3,f4),nrow = 2)
  res = fisher.test(fisher_tb)
  
  return(list(tb = fisher_tb,log2_OR = log2(res$estimate), log10_pval = -log10(res$p.value)))
}

full_snp = as_tibble(fread(file = paste0(snp_dir,'CN20_w_FDR.bed'), 
                           sep = '\t',header = F,na.strings = "NA",
                           col.names = c('chr','start','end','rsID','FDR','strand')))
full_snp = mutate(full_snp, asoc = ifelse(FDR<0.05,1,0))
full_snpRanges = make.ranges(full_snp, c('rsID','FDR'))

merged_footprint = as_tibble(fread(file = paste0(footprint_dir,'CN20_merged_footprint.jaspar2018_hsapiens.dnase2tf_fdr0.05.bed'), 
                                   sep = '\t', header = F,
                                   col.names = c('chr','start','end','TF','score','strand')))
tfRanges = make.ranges(merged_footprint, c('TF'))

dist = distanceToNearest(full_snpRanges,tfRanges)
full_snp$dist = dist@elementMetadata$distance

fisher_stats = fisher_enrich(full_snp,dist_cutoff = 0)

full_snp$asoc = factor(full_snp$asoc,levels = c('1','0'))
full_snp$asoc = plyr::revalue(full_snp$asoc, c('1'='ASoC SNP','0'='non-ASoC SNP'))
ggplot(full_snp, aes(log10(dist+1), stat(width*density), fill=asoc)) + 
  geom_histogram(bins=20, position = 'dodge') +
  labs(x=expression('log'[10]*'(distance of SNP from nearest TF footprint + 1)'),
       y=expression('frequency'), title = 'iN-Glut 20', fill='SNP type') + 
  theme_classic() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16),
                          legend.title = element_text(size=14), legend.text = element_text(size=12),
                          plot.title = element_text(size=16,face = 'bold'))
ggsave('fp_dist_CN20_paired.pdf',plot = last_plot(), width = 6, height = 5)
