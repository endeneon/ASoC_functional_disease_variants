# 01 Dec 2019
# match exclusive peak sets to their corresponding coordination

CN_summits_exclusive <- merge(x = CN_summits_extracted_max, 
                              y = CN_15_exclusive, 
                              by.x = CN_summits_extracted_max$X4, 
                              by.y = CN_15_exclusive$X1)


CN_summits_exclusive <- CN_summits_extracted_max[CN_summits_extracted_max$X4 %in% 
                                                   CN_15_exclusive$X1, ]
CN_summits_exclusive <- CN_summits_exclusive[!duplicated(CN_summits_exclusive$X4), ]

DN_summits_exclusive <- DN_summits_extracted_max[DN_summits_extracted_max$X4 %in% DN_exclusive$X1, ]
DN_summits_exclusive <- DN_summits_exclusive[!duplicated(DN_summits_exclusive$X4), ]

iPS_summits_exclusive <- iPS_summits_extracted_max[iPS_summits_extracted_max$X4 %in% iPS_exclusive$X1, ]
iPS_summits_exclusive <- iPS_summits_exclusive[!duplicated(iPS_summits_exclusive$X4), ]

NPC_summits_exclusive <- NPC_summits_extracted_max[NPC_summits_extracted_max$X4 %in% NPC_exclusive$X1, ]
NPC_summits_exclusive <- NPC_summits_exclusive[!duplicated(NPC_summits_exclusive$X4), ]

GA_summits_exclusive <- GA_summits_extracted_max[GA_summits_extracted_max$X4 %in% GA_exclusive$X1, ]
GA_summits_exclusive <- GA_summits_exclusive[!duplicated(GA_summits_exclusive$X4), ]

write.table(CN_summits_exclusive, 
            file = "CN_summits_exclusive.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(DN_summits_exclusive, 
            file = "DN_summits_exclusive.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(iPS_summits_exclusive, 
            file = "iPS_summits_exclusive.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(NPC_summits_exclusive, 
            file = "NPC_summits_exclusive.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(GA_summits_exclusive, 
            file = "GA_summits_exclusive.bed", 
            quote = F, row.names = F, col.names = F, sep = "\t")
