This folder contains code for three types of analyses on SNP accessibility and cell type specificity.

1. `pairwise_ASoC_sharing`: 
Estiamte the proportion of significant allelic imbalance SNPs share bewteen each pair of cell types using Storey's pi1 analysis.

2. `allelic_imbalance_effect_size`:
For a given pair of cell types, evaluate the correlation of allelic imbalance effect size for cell-type-specific ASoCs.

3. `peak_access_analysis`:
For a given neuronal cell type, obtain neuron-specific ASoC SNPs (w/ p value > 0.05 or DP<20 in iPSC). Compare the diffence of accessibility of peaks that contain these neuron-specific ASoC SNPs in the neuronal cell type vs in iPSC.