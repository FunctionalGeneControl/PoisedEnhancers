This directory contains scripts and files used for Poised Enhancer (PE) Capture Hi-C data analysis.

The capture Hi-C data were analysed as described in the Methods section of Nocente et al., using HiCUP (Wingett et al. 2015), modified to consider both across- and within-read ligation junctions (https://github.com/StevenWingett/HiCUP/tree/combinations, v0.7.4, (Freire-Pritchett et al. 2021)) and chicagoTools v1.30.0 (Cairns et al. 2016, Freire-Pritchett et al. 2021, https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/).

CHiCAGO results were summarised across samples using "makePeakMatrix.R" from ChicagoTools with default parameters.
Example command line:
qsub -v "prefix=Day5_PRC2totalPooled_NPtransition_NE_DE, listfile=Day5_PRC2_NPtransition_NE_DE_2023_2024_peakMatrix_total_pooled.list" RunMakePeakMatrix.sh

See the following scripts and file:
- RunMakePeakMatrix_perChicagoScores_Marina.sh
- RunMakePeakMatrix_perReadCounts_Marina.sh
- Day5_PRC2_NPtransition_NE_DE_2023_2024_peakMatrix_total_pooled.list

To visualise contacts of candidate PEs (baitIDs 580854, 68054, and 1340) in naive and primed hPSCs, the plotBaits() function from the Chicago R package was used, with a maximum interaction distance of 500 kb. 
Run the script "PlotBaits_PECHiC_CnTag_tracks_Figure1_MN.R" with the following input files:
- CHiCAGO RDS files
- peakMatrix_file (Peakmatrix_chicago_score_cutoff5_NPtransition_differentiation.txt), available on OSF.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.
- CUT&Tag bigWig files

Enrichment of PE-interacting regions for histone modifications was estimated using the peakEnrichment4Features() function in CHiCAGO. 
Locus-specific visualisations combining PECHi-C and CUT&Tag signals were generated using the plotgardener R package (v1.6.4; (Kramer et al. 2022)). 
Run the script "Supp_Figure1_final_script_MN.R" in the directory "FigureS1_annotations_enrichment".

Run this script "InteractionClusters_FileGeneration_snippet.R" to generate the file containning the focused PE-PE and PE-TSS cis-interactions, using this input files:
- peakMatrix_file (Peakmatrix_chicago_score_cutoff5_NPtransition_differentiation.txt), available on OSF.
- rmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.rmap), available on OSF.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt), available on OSF.
- TSS_file (2021.09_TSS_1toY_bedFormat.tsv), available in the "GeneOntology_OtherEnds_Figure2F" directory.
