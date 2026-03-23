This directory contains the scripts and files required to reproduce Figure 4 and Figure S2D-F, which illustrate the association between PE chromatin features and their contact dynamics during the naive-to-primed transition (see Methods in Nocente et al.).

To compute H3K27me3 spreading, the primed H3K27me3 CUT&Tag matrix was generated using the script "Heatmaps_DeepTools_K27me3_Primed_final.sh".
This script produces the files "K27me3_heatmap_matrix.txt" and "K27me3_heatmap_matrix_sorted_regions.txt", which are also available on OSF.

To reproduce these plots, run the script "Figure4_SuppFigS2D-F_final_script_MN_12-02-2026.R" using the following input files:
- H3K27me3_matrix_primed (K27me3_heatmap_matrix.txt), available on OSF.
- sorted_regions_primed (K27me3_heatmap_matrix_sorted_regions.txt), available on OSF.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap), available on OSF.
- peakMatrix (Focused_peakmatrix_PE_TSS_interaction.tsv), available on OSF.
- countsNorm (Counts_DESeq2norm_fragLenNorm_ALL.txt), available on OSF.
- lassoPred (2209_LassoPredictorsMatrix_FragLenNorm_Scaled.tsv), available on OSF.
- Temporalcluster_file (pepm_transition_dtw_final_02-02-2026.txt), available in the "Temporal_clustering_PEs_contacts" directory.