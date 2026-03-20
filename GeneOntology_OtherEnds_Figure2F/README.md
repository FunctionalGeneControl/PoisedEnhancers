This directory contains the script and files required to perform Gene Ontology enrichment analysis of genes contacted by PEs (other-end TSS) across different contact dynamics during the transition (see Methods in Nocente et al.).

To run the script, the following input files are required:
- dtwClnew_file (pepm_transition_dtw_final_02-02-2026.txt), produiced in the directory "Temporal_clustering_PEs_contacts".
- restriction_fragm_coor_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.rmap), available on OSF
- TSS_file (2021.09_hg38_TSSmain.tsv), available in this directory
- TSS_adjusted_coor_file (2021.09_TSS_1toY_bedFormat.tsv), available in this directory

The script generates the file "OtherEnd_TSS_coordinates.txt".