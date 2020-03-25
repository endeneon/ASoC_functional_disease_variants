# Scatter plots of chromatin accessibility in neuronal cell types vs in iPSC
# of peaks that contain neuron-specific ASoC SNPs (w/ pval>0.05 or DP<20 in iPS)
# Fig 1E, Fig S10


library(data.table)
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
options(stringsAsFactors = F)

# Obtain the chromatin accessibility of neuron-specific ATAC-seq peaks in various cell types ####
make.ranges = function(df, metadata){
  dfRanges <- GRanges(seqnames = df$chr, 
                      ranges = IRanges(start = df$start,
                                       end = df$end))
  for(m in metadata){
    dfRanges <- plyranges::mutate(dfRanges, !!m := unname(df[,m]))
  }
  return(dfRanges)
}

peak_dir = '~/Google Drive/Research/DNA_footprint_and_ATAC-seq/ATAC-seq_PeakData/new/'
peak_counts.rpkm_log2= readRDS(paste0(peak_dir,'peak_counts.log2_RPKM.rds'))
peak_info = data.frame(fread(paste0(peak_dir,'all_peaks.gcContent.txt'),sep='\t',header = T))
peakRanges = make.ranges(peak_info, 'peak_name')

ASoC_spec.access = list()
for (A in c('NSC','CN','DN','GA')){
  print(paste('Cell type:',A))
  infile <- paste0('../celltype_spec_ASoC/ASoC_',A,'_specific.IPS.bed')
  ASoC_spec.not_ips <- read.table(infile,header = T)
  names(ASoC_spec.not_ips)[1:3] = c('chr','start','end')
  snpRanges = make.ranges(ASoC_spec.not_ips, 'ID')
  
  dist = distanceToNearest(snpRanges,peakRanges)
  ASoC_spec.not_ips$dist = dist@elementMetadata$distance
  ASoC_spec.not_ips$peak_name = peak_info$peak_name[dist@to]
  ASoC_spec.not_ips = filter(ASoC_spec.not_ips,dist==0)
  peak_count.A_ips = cbind(select(peak_counts.rpkm_log2, starts_with(A)), 
                           select(peak_counts.rpkm_log2, starts_with('ips')))
  ASoC_spec.not_ips[[A]] = rep(NA,nrow(ASoC_spec.not_ips))
  ASoC_spec.not_ips$ips = rep(NA,nrow(ASoC_spec.not_ips))
  
  for (i in 1:nrow(ASoC_spec.not_ips)){
    pk_name = ASoC_spec.not_ips$peak_name[i]
    count.tmp = as.numeric(peak_count.A_ips[pk_name,])
    ASoC_spec.not_ips[[A]][i] = mean(count.tmp[1:8])
    ASoC_spec.not_ips$ips[i] = mean(count.tmp[9:16])
  }
  ASoC_spec.access[[A]] = ASoC_spec.not_ips
}

saveRDS(ASoC_spec.access,'data/neuron_spec_ASoC_peak_access_rpkm.rds')

# Plotting (Fig 1E, Fig S10) ####
ASoC_spec.access = readRDS('data/neuron_spec_ASoC_peak_access_rpkm.rds')
convert.tb = data.frame(fname = names(ASoC_spec.access), 
                        formal = c('NPC','iN-Glut','iN-DN','iN-GA'))

for (i in 1:nrow(convert.tb)){
  A = convert.tb$fname[i]
  A_formal = convert.tb$formal[i]
  ASoC_spec.not_ips = ASoC_spec.access[[A]]
  # Remove duplicated peaks:
  ASoC_spec.not_ips = ASoC_spec.not_ips[,8:10]
  ASoC_spec.not_ips = unique(ASoC_spec.not_ips)
  # Calculate % of SNPs in the shaded area:
  access_diff = ASoC_spec.not_ips$ips - ASoC_spec.not_ips[[A]]
  percentage = sum(access_diff>-1 & access_diff<1) / nrow(ASoC_spec.not_ips)
  # NSC: 44.3%, CN: 47.6%, DN: 49.9%, GA: 56.9%
  
  ## Overlap with differential stats obtained from 'differential_peak_analysis.R':
  # peak_diff = data.frame(fread(paste0('peak_DE_result/',A,'_vs_iPS_peak_edgeR_res.txt'),
  #                              sep = '\t',header = T),row.names = 1)
  # peak_diff$peak_name = rownames(peak_diff)
  # ASoC_spec.not_ips = inner_join(ASoC_spec.not_ips,peak_diff,by='peak_name')
  
  xmax = max(ASoC_spec.not_ips[[A]])
  ymax = max(ASoC_spec.not_ips$ips)
  
  scatterplot = ggplot(ASoC_spec.not_ips,aes_string(x=A,y='ips')) + geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'blue') +
    labs(title = paste(A_formal,'specific:',nrow(ASoC_spec.not_ips))) +
    scale_x_continuous(name = bquote('Accessibility in'~.(A_formal)~'('~log[2]~'ATAC-seq peak RPKM'~')'),
                       expand = c(0, 0), breaks = seq(0,6,1), limits = c(0,7+0.1)) +
    scale_y_continuous(name = bquote('Accessibility in iPS ('~log[2]~'ATAC-seq peak RPKM'~')'),
                       expand = c(0, 0), limits = c(0,ymax+0.1)) +
    theme_classic() + theme(axis.title = element_text(size = 16),
                            axis.text = element_text(size = 14)) +
    annotate('polygon',x=c(0,0,1,ymax+1.1,ymax-0.9), y=c(1,0,0,ymax+0.1,ymax+0.1),
             fill="blue", alpha=0.2) +
    geom_text(x=0.5,y=0.5,label=paste0(signif(percentage*100,digits = 3),'%'),size=5)
  
  ggsave(paste0('Figs/',A,'_spec_ASoC_peak_access.rpkm_vs_iPS.pdf'),
         scatterplot,width = 8,height = 7)
}
