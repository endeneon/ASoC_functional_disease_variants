# Following is the code to obtain the disruption effect of 
# ASoC SNPs of a given cell type on JASPAR TF motifs 
# using the R package 'motifbreakR',
# (and to generate the motif logo figures at selected SNPs.)


library(data.table)
library(dplyr)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP149.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)

options(stringsAsFactors = FALSE)
SNP_dir = '~/Downloads/ASoC/data/ASoC_per_celltype/'
footprint_dir = '~/Downloads/ASoC/TF_motifs/footprint_jaspar_hsap_redundant/'

JASPAR2018.indx = mcols(MotifDb)$dataSource=='jaspar2018'
JASPAR2018.motif = MotifDb[JASPAR2018.indx]
JASPAR2018.Hsapiens.motif = JASPAR2018.motif[mcols(JASPAR2018.motif)$organism=='Hsapiens']
rm(JASPAR2018.motif,JASPAR2018.indx)

selected.snps = snps.from.rsid(rsid = c('rs2027349','rs12895055','rs11624408'),
                               dbSNP = SNPlocs.Hsapiens.dbSNP149.GRCh38,
                               search.genome = BSgenome.Hsapiens.UCSC.hg38)
selected.mb.result = motifbreakR(snpList = selected.snps, pwmList = JASPAR2018.Hsapiens.motif,
                                 method = 'ic',threshold = 2.5e-4, filterp = TRUE, show.neutral = TRUE)
# Motif disruption at rs2027349 (Fig 4E)
plotMB(results = selected.mb.result, rsid = "rs2027349", effect = c("strong","weak"))
# Motif disruption at rs12895055 (Fig S25C)
plotMB(results = selected.mb.result, rsid = "rs12895055", effect = "strong")
# Motif disruption at rs11624408 (Fig S25D)
plotMB(results = selected.mb.result, rsid = "rs11624408", effect = "strong")

## Run motifbreakR() for all SNPs of interest in a given cell type and save the results:

make.ranges = function(df, metadata){
  dfRanges = GRanges(seqnames = df$chr, 
                     ranges = IRanges(start = df$start,
                                      end = df$end))
  for(m in metadata){
    dfRanges = plyranges::mutate(dfRanges, !!m := unname(df[,m]))
  }
  return(dfRanges)
}

for (celltype in c('CN20','NSC20')){
  ## Load SNPs of interest
  SNP.stats = data.frame(fread(file = paste0(SNP_dir,celltype,'_w_FDR.txt'),
                               sep = '\t',header = T,na.strings = "NA"))
  names(SNP.stats)[1:2] = c('chr','end')
  SNP.stats = mutate(SNP.stats,start = end-1)
  snpRanges = make.ranges(SNP.stats, 'rsID')
  
  ## Load imputed footprints from all TFs
  merged_footprint = data.frame(fread(file = paste0(footprint_dir,celltype,'_merged_footprint.jaspar2018_hsapiens.dnase2tf_fdr0.05.bed'),
                                      sep = '\t',header = F,col.names = c('chr','start','end','TF','score','strand')))
  fp.Ranges = make.ranges(merged_footprint, 'TF')
  
  ## Measure SNP distance to nearest TF footprints
  dist = distanceToNearest(snpRanges,fp.Ranges)
  SNP.stats$dist_fp = dist@elementMetadata$distance
  rm(snpRanges,merged_footprint,fp.Ranges,dist)
  
  ## SNP filtering
  # (Number of all SNPs) CN20: 111718; NSC20: 99402
  SNP_inFP.stats = filter(SNP.stats, dist_fp==0) # CN20: 25660; NSC20: 22071
  ASoC.stats = filter(SNP.stats, FDR<0.05) # CN20: 5611; NSC20: 3547
  
  # ASoC SNPs in TF footprints:
  ASoC_inFP.stats = filter(ASoC.stats, dist_fp==0) # CN20: 1802; NSC20: 1201
  snps.obj = snps.from.rsid(rsid = ASoC_inFP.stats$rsID,
                            dbSNP = SNPlocs.Hsapiens.dbSNP149.GRCh38,
                            search.genome = BSgenome.Hsapiens.UCSC.hg38)

  ## Motif disruption analysis
  mb.result = motifbreakR(snpList = snps.obj, pwmList = JASPAR2018.Hsapiens.motif,
                          method = 'ic',threshold = 2.5e-4, filterp = TRUE, show.neutral = TRUE, 
                          BPPARAM = BiocParallel::MulticoreParam(workers = 16))
  
  saveRDS(mb.result, paste0('data/',celltype,'_in_called_fp.motifBreak.rds'))
}
