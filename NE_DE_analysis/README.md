This directory contains a script to reproduce Figures 6B-C and Figure S4, which show the dynamics of PE contacts after developmental activation (differentiation in neuroectoderm (NE) or definitive endorderm (DE); see Methods in Nocente et al. for details).

To reproduce these figures, run the script "Figure6_NE_DE_analysis_final_script_v17-02-2026.R", using the following input files:
- H3K4me1, H3K27me3, H3K4me3, and H3K27ac CUT&Tag pileup files generated in the "CUTnTag_transition" directory and available on OSF.
- peakmatrixFull_file (Day5_PRC2totalPooled_NPtransition_NE_DE.txt), available on OSF.
- Focused_peakmatrix_file (Focused_peakmatrix_PE_TSS_interaction.tsv), available on OSF.
- cluster_file (pepm_transition_dtw_final_02-02-2026.txt), available in the "Temporal_clustering_PEs_contacts" directory.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.

This script generates the following output files:
- "Supplementary_Table_H3K27me3.tsv", "Supplementary_Table_H3K27ac.tsv", and "Supplementary_Table_H3K4me1.tsv", which integrate PECHi-C stats, temporal classes and CUT&Tag DESEq2 results (DE vs Day14 or NE vs Day14), and are used to generate Table S4.