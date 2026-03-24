This directory contains scripts to reproduce Figures 5D-E, Figure S3F, and Figures 5C and S3C, which present integrative analyses of CUT&Tag and PECHi-C data following PRC2 perturbations, as well as representative examples (see Methods in Nocente et al. for details).

Reproducing Figures 5D–E and Figure S3F:
These figures show changes in H3K27me3 signal (D) or in connectivity (E) at PEs after treatment with PRC2 PROTAC or PRC2 catalytic inhibitor relative to DMSO, stratified by bins of log10-H3K27me3 signal in the DMSO condition. Figure S3F shows changes in the overall connectivity of PEs after treatments.
To reproduce these figures, run the script "Integration_CUTnTag_PRC2data_with_PECHiC_Fig5DE_FigS3F_vf.R" using the input files:
- res_PROTAC_vs_DMSO_merged (res_deseq2PECHiC_PROTAC_DMSO_merged_all.txt), which are the PECHi-C DESeq2 results for PROTAC vs DMSO (generated in the "PECHiC_DESeq2_PRC2_perturbation" directory and available on OSF).
- res_UNC1999_vs_DMSO_merged (res_deseq2PECHiC_UNC_DMSO_merged_all.txt), which are the PECHi-C DESeq2 results for UNC1999 vs DMSO (generated in the "PECHiC_DESeq2_PRC2_perturbation" directory and available on OSF).
- temporal_clusters_file (pepm_transition_dtw_final_02-02-2026.txt), available in the "Temporal_clustering_PEs_contacts" directory.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.
- res_PROTAC_DMSO_CnT_all_dmsoBin_normal (CUTnTag_deseq2_res_PROTAC_DMSO_dmsoMean_RunSeq2.txt), which are the CUTnTag DESeq2 results for PROTAC vs DMSO (generated in the "CUTnTag_PRC2perturbations_analysis" directory and available on OSF).
- res_UNC1999_DMSO_CnT_all_dmsoBin_normal (CUTnTag_deseq2_res_UNC1999_DMSO_dmsoMean_RunSeq2.txt), which are the CUTnTag DESeq2 results for UNC1999 vs DMSO (generated in the "CUTnTag_PRC2perturbations_analysis" directory and available on OSF).

Reproducing Figures 5C and S3C
These figures show representative examples of loss of PECHi-C contacts following treatment with a PRC2 PROTAC or a PRC2 catalytic inhibitor, compared to DMSO, alongside H3K27me3 and H2AK119ub CUT&Tag signals for each condition.
To reproduce these figures, run the script "Figure5C_SuppFigS3C_PRC2perturb_examples_script_final_MN_v16-02-2026.R" with the following input files:
- Peakmatrix_all_file (Day5_PRC2totalPooled_NPtransition_NE_DE.txt), available on OSF.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap), available on OSF.
- CUTnTag bigwig merged files