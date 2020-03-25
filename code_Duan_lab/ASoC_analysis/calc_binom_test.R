# Siwei 07 Dec 2018
# make new GA set use new seq data

i <- 1
for(i in 1:nrow(new_GA_SNP)) {

  new_GA_SNP$pV[i] <- binom.test(x = new_GA_SNP$ALT_C[i], n = new_GA_SNP$DP[i], p = 0.5, 
                                 alternative = "two.sided")$p.value
  print(paste("i =", i, new_GA_SNP$pV[i], 
              collapse = " ", sep = " "))
}

new_GA_SNP$FDR <- p.adjust(new_GA_SNP$pV, "fdr")
new_GA_SNP$Bonferroni <- p.adjust(new_GA_SNP$pV, "bonferroni")
