# Pairwise pi1 analysis across 5 cell types (Fig 1D)


library(dplyr)
library(ggplot2)

get_pi1 <- function(pval,thres=0.5){
  pi1 <- (length(pval) - sum(pval>thres)/(1-thres)) / length(pval)
  return(pi1)
}

pi1_A_in_B <- function(A.full,B.full){
  A.sig <- filter(A.full, FDR<0.05)$ID
  A.sig_in_B <- filter(B.full, ID %in% A.sig)
  return(get_pi1(A.sig_in_B$pBinom) * nrow(A.sig_in_B) / length(A.sig))
  # pi1 = 1-qvalue::pi0est(A.sig_in_B$pBinom)$pi0
  # return(pi1 * nrow(A.sig_in_B) / length(A.sig))
}

cell_types <- c('iPS_8','NSC_8','CN_8','DN_8','GA_8')
cell_names <- c('iPS','NPC','iN-Glut','iN-DN','iN-GA')
pairwise_pi1 <- data.frame(diag(5), row.names = cell_names)
names(pairwise_pi1) <- cell_names

for (i in 1:4){
  for (j in (i+1):5) {
    A <- cell_types[i] # 'iPS'
    B <- cell_types[j] # 'CN'
    print(paste('Pair:',A,B))
    A.full <- readRDS(paste0('../SNP_accessibilty/',A,'_all_SNPs.rds'))
    B.full <- readRDS(paste0('../SNP_accessibilty/',B,'_all_SNPs.rds'))
    pairwise_pi1[i,j] <- pi1_A_in_B(A.full,B.full)
    pairwise_pi1[j,i] <- pi1_A_in_B(B.full,A.full)
  }
}

pairwise_pi1.plot <- expand.grid(dimnames(pairwise_pi1))
pairwise_pi1.plot$value <- unlist(unname(pairwise_pi1))
names(pairwise_pi1.plot) <- c('Leading','Matched','pi1')
ggplot(data = pairwise_pi1.plot, aes(x=Leading, y=Matched)) + geom_tile(aes(fill=pi1)) + 
  scale_fill_gradient(high = 'firebrick3', low = 'white',# low = "deepskyblue4", high = "lightskyblue1",
                       name = expression(pi*'1')) + 
  theme_minimal() + xlab('Leading Cell Type') + ylab('Matched Cell Type') + 
  theme(axis.title = element_text(size = 14,face="bold"), axis.text = element_text(size = 14,face="bold"),
        legend.title = element_text(size = 13,face="bold"), legend.text = element_text(size = 12))
ggsave('Figs/Pairwise_pi1.pdf',last_plot(),width = 7,height = 5)
