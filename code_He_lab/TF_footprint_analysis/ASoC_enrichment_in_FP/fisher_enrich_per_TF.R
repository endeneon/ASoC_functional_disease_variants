# Following is the code to analyze the enrichment of ASoC SNPs vs non-ASoC SNPs 
# in identitified TF footprints per cell type,
# and to generate Figures 2G and S18 -- 
# "Enrichment of TF-binding footprints in ASoC SNPs in each cell type"


library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = F)
snp_dir = '~/Downloads/ASoC/data/ASoC_per_celltype/'
footprint_dir = '~/Downloads/ASoC/TF_motifs/footprint_jaspar_hsap_redundant/'

make.ranges = function(df, metadata){
  dfRanges = GRanges(seqnames = df$chr, 
                      ranges = IRanges(start = df$start,
                                       end = df$end))
  for(m in metadata){
    dfRanges = plyranges::mutate(dfRanges, !!m := unname(df[,m]))
  }
  return(dfRanges)
}

fisher_enrich = function(snp.df, dist, dist_cutoff = 0){
  snp.df$dist = dist
  ASoC_snp = filter(snp.df,asoc==1)
  nonASoC_snp = filter(snp.df,asoc==0)
  
  f1 = sum(ASoC_snp$dist <= dist_cutoff)
  f3 = dim(ASoC_snp)[1]-f1
  f2 = sum(nonASoC_snp$dist <= dist_cutoff)
  f4 = dim(nonASoC_snp)[1]-f2
  fisher_tb = matrix(c(f1,f2,f3,f4),nrow = 2)
  res = fisher.test(fisher_tb)

  return(list(log2_OR = log2(res$estimate), log10_pval = -log10(res$p.value)))
}

enrich_per_celltype = list()
for (celltype in c('CN','NSC','DN','GA','IPS','CN20','NSC20')){
  if (celltype %in% c('CN20','NSC20')){
    ocr_snp = as_tibble(fread(file = paste0(snp_dir,celltype,'_w_FDR.txt'),
                              sep = '\t',header = T,na.strings = "NA"))
  } else {
    ocr_snp = as_tibble(fread(file = paste0(snp_dir,celltype,'_core8_w_FDR.txt'),
                              sep = '\t',header = T,na.strings = "NA"))
  }
  
  names(ocr_snp)[c(1,2)] = c('chr','end')
  ocr_snp = mutate(ocr_snp, start = end-1) %>% mutate(asoc = (FDR<0.05)*1 )
  snpRanges = make.ranges(ocr_snp, c('rsID','asoc'))
  
  merged_footprint = as_tibble(fread(file = paste0(footprint_dir,celltype,'_merged_footprint.jaspar2018_hsapiens.dnase2tf_fdr0.05.bed'),
                                     sep = '\t',header = F,col.names = c('chr','start','end','TF','score','strand')))
  motif = unique(merged_footprint$TF)
  enrich_by_TF = data.frame(matrix(nrow = length(motif), ncol = 2), row.names = motif)
  names(enrich_by_TF) = c('log2_OR','log10_pval')
  for (m in motif){
    single_footprint = filter(merged_footprint, TF==m)
    tfRanges = make.ranges(single_footprint, c('TF'))
    dist = distanceToNearest(snpRanges,tfRanges)
    dist = dist@elementMetadata$distance
    
    fisher_stats = fisher_enrich(ocr_snp,dist,dist_cutoff = 0)
    enrich_by_TF[m,] = c(fisher_stats$log2_OR,fisher_stats$log10_pval)
  }
  enrich_by_TF$motif = row.names(enrich_by_TF)
  enrich_by_TF$TF = sapply(enrich_by_TF$motif, function(x) {strsplit(x,split = '[.]')[[1]][3]})
  
  tmp.df = enrich_by_TF %>% group_by(TF) %>% summarise(motif = max(motif))
  enrich_by_TF.simp = inner_join(enrich_by_TF, tmp.df, by = c('TF','motif'))
  
  enrich_per_celltype[[celltype]] = enrich_by_TF.simp
  
  plot_out = ggplot(enrich_by_TF.simp, aes(x=log10_pval, y=log2_OR, label = ifelse(log10_pval>3,TF,''))) +
    geom_point(color = ifelse(enrich_by_TF.simp$log10_pval>3,'red','black')) + geom_text_repel(size = 3) +
    geom_hline(yintercept = 0, color='blue') + geom_vline(xintercept = 3, color='brown1', linetype="dashed") +
    labs(x = '-log10(p value)', y = 'log2(odds ratio)', title = paste(celltype, 'ASoC enrichment in TF footprints')) +
    theme_classic()
  ggsave(paste0('Figs/',celltype,'_ASoC_enrich_per_TFfp.pdf'), plot = plot_out, width = 8, height = 6)
}
