This folder contains code for three types of analyses on SNPs using predicted TF footprint data.

1. `ASoC_to_FP_dist_analysis`: 
For a given cell type, measure the distance of ASoC SNPs from their nearest predicted TF footprints.

2. `ASoC_enrichment_in_FP`:
For a given cell type and a given type of TF, evaluate the enrichment level of ASoC vs non-ASoC SNPs based on whether a SNP is in a TF footprint or not using Fisher's exact test.

3. `motifbreak_corr_analysis`:
For a given cell type and a given type of TF, use 'motifbreakR' to measure the disruption effect of ASoC SNPs on the given TF motif, and correlate it with the level of allelic imbalance of accessiblity at the SNPs.