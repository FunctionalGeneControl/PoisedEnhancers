This directory contains the script required to reproduce Figure 2E, showing the connectivity of PEs with rewired contacts (see Methods in Nocente et al.).

To reproduce this plot, run the script "NetworkAnalysis_pub_MS.Rmd" using the following input files:
- peakMatrix (Peakmatrix_chicago_score_naive-to-primed_transition_differentiation_allSamples.txt), available on OSF.
- IntClOG (pepm_transition_dtw_final_02-02-2026.txt), available in the "Temporal_clustering_PEs_contacts" directory.
- baitmap (2019.11.15_hg38_PE.baitmap), available on OSF.
- baitmap2 (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.
- rmap (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.rmap), available on OSF.
- TSS (2021.09_hg38_TSSmain.tsv), available in the "GeneOntology_OtherEnds_Figure2F" directory .
- TSSv1 (2021.09_TSS_1toY_bedFormat.tsv), available in the "GeneOntology_OtherEnds_Figure2F" directory.