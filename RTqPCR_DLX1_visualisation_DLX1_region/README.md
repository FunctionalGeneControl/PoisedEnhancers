This directory contains a script to reproduce Figures 6E–F and Figure S4E. These figures show normalised gene expression levels of DLX1, DLX2, and METAP1D following inducible CRISPR activation (Figure 6F), as well as visualisations of H3K27me3 and H3K4me1 CUT&Tag signals and the “gained-early” PECHi-C contact in the DLX1–DLX2 region across the naïve-to-primed transition (Figures 6E and S4E) (see Methods in Nocente et al. for details).

To reproduce Figures 6E-F and Figure S4E, run the following script "Figure6_SuppFigS4_Dlx1_CRISPRa_final_script_MN_17-02-2026.R" using these input files:
- Peakmatrix_file (Day5_PRC2totalPooled_NPtransition_NE_DE.txt), available on OSF.
- TssOE_file (OtherEnd_TSS_coordinates.txt), available in this directory and generated in the "GeneOntology_OtherEnds_Figure2F" directory.
- baitmap_file (Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap), available on OSF.