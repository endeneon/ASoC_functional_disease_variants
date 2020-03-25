# 23 Dec 2018
# Siwei read de novo assembled and called GAs

library(DiffBind)
library(gplots)
library(factoextra)

#
#######
peak_new_GA_counted <- dba()
GA_read_matrix_23Dec2018 <- dba()
# colnames(GA_read_matrix_23Dec2018) <- c("SampleID",	"Tissue",	"Peaks",	"bamReads",	
#                                         "PeakCaller",	"PeakFormat",	"ScoreCol",	"LowerBetter")

peaks_new_GA <- dba(sampleSheet = GA_read_matrix_23Dec2018,
                    config = data.frame(RunParallel = TRUE),
                    bRemoveM = TRUE)
#

####### Occupancy test

peaks_new_GA_minOverlap_4_peakset <- dba.peakset(peaks_new_GA, minOverlap = 0.5,
                                                 consensus = c(DBA_TISSUE),
                                                 bRetrieve = F)


peaks_new_GA_minOverlap_4_consensus_tissue <- dba(peaks_new_GA_minOverlap_4_peakset,
                                                  mask = peaks_new_GA_minOverlap_4_peakset$masks$Consensus,
                                                  minOverlap = 1)
