# 21 Dec 2018
# extract the count matrix from 5x8 table

library(edgeR)
library(DiffBind)
peak_5x8_raw_counts <- data.frame(1:241911)
peak_5x8_index <- data.frame(1:241911)

peak_occupancy_consensus_edgeR <- dba()
i <- 1

# peak_5x8_raw_counts[[1]] <- peak_occupancy_consensus_edgeR$peaks[[1]]$Reads
for(i in 1:40) {
  peak_5x8_raw_counts[[i]] <- peak_occupancy_consensus_edgeR$peaks[[i]]$Reads
}
peak_5x8_raw_counts$Geneid_hg38 <- paste(peak_occupancy_consensus_edgeR$peaks[[1]]$Chr,
                        peak_occupancy_consensus_edgeR$peaks[[1]]$Start,
                        peak_occupancy_consensus_edgeR$peaks[[1]]$End,
                        sep = "_")

save(peak_5x8_raw_counts, file = "peak_5x8_counts.RData")

# for(i in 1:40) {
#   peak_5x8_raw_counts$Geneid_hg38 <-
#     paste(peak_5x8_raw_counts)
#
# }


colnames(peak_5x8_raw_counts) <- DiffBind_read_matrix_29Nov2018$SampleID
save(peak_5x8_raw_counts, file = "peak_5x8_counts.RData")

paste("abc", "ABC", sep = "_", collapse = "_")
