# 29 Jan 2019 Siwei
# make DE gene list for each group

library(edgeR)
library(gplots)
library(stats)
library(corrplot)
library(ggplot2)
library(factoextra)

###

group <- factor(c(rep("CN", 8), rep("DN", 8), rep("iPS", 8), rep("NPC", 8), rep("GA", 8)))
raw_DGE <- DGEList(counts = as.matrix(raw_reads), remove.zeros = T)


DGE_filtered <- raw_DGE[rowSums(cpm(raw_DGE) > 1) >= 4, , keep.lib.sizes = F]
filtered_gene_names <- as.data.frame(rownames(DGE_filtered$counts))

DGE_filtered <- calcNormFactors(DGE_filtered)
DGE_filtered <- estimateDisp(DGE_filtered)
design <- model.matrix(~ 0 + group)
DGE_fitted <- glmFit(DGE_filtered, design = design) # Use likelihood ratio test glmFit() and glmLRT()
# > design
# groupCN groupDN groupGA groupiPS groupNPC
#####

NPC_vs_iPS <- glmLRT(DGE_fitted, contrast = c(0, 0, 0, 1, -1))
decideTests(NPC_vs_iPS)
summary(decideTests(NPC_vs_iPS))
# 1*groupiPS -1*groupNPC
# Down                     4484
# NotSig                   9617
# Up                       4117

NPC_vs_iPS_result <- as.data.frame(decideTests(NPC_vs_iPS))
NPC_vs_iPS_result$Geneid <- rownames(NPC_vs_iPS_result)
NPC_vs_iPS_result_sig <- NPC_vs_iPS_result[NPC_vs_iPS_result$`1*groupiPS -1*groupNPC` != 0, ]

NPC_vs_iPS_FC <- NPC_vs_iPS$table[rownames(NPC_vs_iPS$table) %in% NPC_vs_iPS_result_sig$Geneid, ]
NPC_vs_iPS_FC <- NPC_vs_iPS$table[rownames(NPC_vs_iPS$table) %in% NPC_vs_iPS_TSS_1K_table_sig$X1, ]
write.table(NPC_vs_iPS_FC, file = "NPC_vs_iPS_FC.txt", sep = "\t", quote = F)

NPC_vs_iPS_TSS_1K_table_sig <- gencode_v28_TSS_1K_import[gencode_v28_TSS_1K_import$X1 %in% 
                                                           rownames(NPC_vs_iPS_FC), ]
NPC_vs_iPS_TSS_1K_table_sig <- NPC_vs_iPS_TSS_1K_table_sig[!duplicated(NPC_vs_iPS_TSS_1K_table_sig$X1), ]
write.table(NPC_vs_iPS_TSS_1K_table_sig, file = "NPC_vs_iPS_TSS_1K_table_sig.saf", sep = "\t", 
            row.names = F, col.names = F, quote = F)
#####

Glut_vs_NPC <- glmLRT(DGE_fitted, contrast = c(1, 0, 0, 0, -1))
decideTests(Glut_vs_NPC)
summary(decideTests(Glut_vs_NPC))

# 1*groupCN -1*groupNPC
# Down                     111
# NotSig                 17461
# Up                       646
#####

GA_vs_NPC <- glmLRT(DGE_fitted, contrast = c(0, 0, 1, 0, -1))
decideTests(GA_vs_NPC)
summary(decideTests(GA_vs_NPC))

# 1*groupGA -1*groupNPC
# Down                    4054
# NotSig                  9362
# Up                      4802
#####

DN_vs_NPC <- glmLRT(DGE_fitted, contrast = c(0, 1, 0, 0, -1))
decideTests(DN_vs_NPC)
summary(decideTests(DN_vs_NPC))

# 1*groupDN -1*groupNPC
# Down                    3244
# NotSig                 10477
# Up                      4497

#####################



colnames(filtered_gene_names) <- "Geneid"
filtered_gene_names <- merge(filtered_gene_names, read_length, by = "Geneid")

filtered_RPKM <- rpkm.DGEList(DGE_filtered, gene.length = filtered_gene_names$Length, log = T, prior.count = 1)
filtered_cpm <- cpm.DGEList(DGE_filtered, log = T)
