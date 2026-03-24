This directory contains scripts to analyse CUT&Tag data following PRC2 perturbations (PROTAC degrader or UNC1999 catalytic inhibitor). See Methods in Nocente et al. for details.

For quantitative analysis, individual calibrated H3K27me3 bigWig files were converted into pileUp files using this script: "01-bwigToPileup_CUTnTag_PRC2perturbations.sh" and the baitmap file available on OSF.

These pileUp files (available on OSF) were tested for differential counts (each perturbation versus DMSO) using the R package DESeq2 with the script "DESeq2_CutnTag_PRC2_perturbations_final_RunSeq2_15-02-2026.Rmd".
The output files are:
- "CUTnTag_deseq2_res_PROTAC_DMSO_dmsoMean_RunSeq2.txt", available on OSF.
- "CUTnTag_deseq2_res_UNC1999_DMSO_dmsoMean_RunSeq2.txt", available on OSF.

To reproduce Figure 5B, run the script "02-Heatmaps_DeepTools_PRC2perturbationsv23032026.sh" with the input file:
- "regions_chr_PE_baits.bed", available in this directory.