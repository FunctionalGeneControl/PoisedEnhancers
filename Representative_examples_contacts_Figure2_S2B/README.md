This directory contains the script required to reproduce Figures 2A (right panel), 2C, and S2B, showing an example of PE-CHiC significant contacts and CUT&Tag signal (Figure 2A), as well as representative examples of PE contacts during the naive-to-primed transition (Figures 2C and S2B).

To reproduce these plots, run the script "Figure2_script_representative_examples_PE_contacts.R" using the flolling input files:
- dtw_clusters_file (pepm_transition_dtw_final_02-02-2026.txt), available in the "Temporal_clustering_PEs_contacts" directory.
- Peakmatrix_raw_score_file (Day5_PRC2totalPooled_NPtransition_NE_DE.txt), available on OSF.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.
- tss_file (OtherEnd_TSS_coordinates.txt), produced by the script "Figure2F_PECHiC_GOTermAnalysis_MN_09-02-2026.Rmd" in the "GeneOntology_OtherEnds_Figure2F" directory.

