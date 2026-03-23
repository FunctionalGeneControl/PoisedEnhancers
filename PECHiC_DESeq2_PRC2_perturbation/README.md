This directory contains a script for performing differential analysis of read counts obtained from PECHi-C experiments following DMSO, PRC2 PROTAC, or PRC2 catalytic inhibitor treatment (see Methods in Nocente et al.).

Raw PECHi-C data for the three conditions (three replicates each of DMSO, PROTAC, and UNC1999) were aligned using HiCUP combinations and processed with chicagoTools, as described in Methods in Nocente et al.

For the differential analysis, two replicates for each condition (PROTAC and UNC1999) were merged to balance sequencing depth across samples. Total connectivity in cis per PE bait in each replicate was computed using the custom R script "Read_counts_perBaits_tables_CIS_chinput_Marina.R" with the input list "Input_list_files_Read_counts_PerBaits.txt". This script generates the file "Read_counts_PerBait_CIS_table_Chinput.txt", which is used to assess changes in connectivity upon each PRC2 perturbation using DESeq2 in the script "DESeq2_PECHiC_PRC2-perturbation_read-counts-CIS_analysis_final_13-02-2026.Rmd".
The differential analysis results are available on OSF:
- "res_deseq2PECHiC_PROTAC_DMSO_merged_all.txt" (PROTAC vs DMSO)
- "res_deseq2PECHiC_UNC_DMSO_merged_all.txt" (PRC2 catalytic inhibitor vs DMSO)
