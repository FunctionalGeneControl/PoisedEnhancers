################################################################################
### Examples of loss of PECHi-C contacts following treatment with a PRC2 PROTAC or a PRC2 catalytic inhibitor, compared to DMSO - Figure 5D and Supplemental Figure S3F
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 5D and Supplemental Figure S3F:
###  Top panel: Gene locations, with arcs showing PE contacts in each condition. 
###  Middle and bottom panels: H3K27me3 (middle) and H2AK119ub (bottom) CUT&Tag signals in each condition. 
################################################################################

########################################
### Load required libraries
########################################
library(data.table)
library(GenomicRanges)
library(plotgardener)
library(grid)
library(rtracklayer)

########################################
### User configuration
########################################
base_dir <- "~/Documents/Bioinformatics"
working_dir <- file.path(base_dir, "Find_examples")

Peakmatrix_all_file <- file.path(base_dir, "PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
baitmap_file <- file.path(base_dir, "Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")

# H3K27me3 and H2AK119Ub CUT&Tag files
CUTnTag_dir <- file.path(base_dir,"CnT_PRC2_perturb_Day5/Run2_CUTnTagPrc2_June2025/")

DMSO_K27me3_CnT <- file.path(CUTnTag_dir, "Calibrated_Bigwig_from_Andrew/MergedBigwigs/DMSO_H3K27me3.merged.bw")
PROTAC_K27me3_CnT <- file.path(CUTnTag_dir, "Calibrated_Bigwig_from_Andrew/MergedBigwigs/PROTAC_H3K27me3.merged.bw")
UNC_K27me3_CnT <- file.path(CUTnTag_dir, "Calibrated_Bigwig_from_Andrew/MergedBigwigs/UNC_H3K27me3.merged.bw")

DMSO_K119Ub_CnT <- file.path(CUTnTag_dir, "Calibrated_Bigwig_from_Andrew/MergedBigwigs/DMSO_H2AK119Ub.merged.bw")
PROTAC_K119Ub_CnT <- file.path(CUTnTag_dir,"Calibrated_Bigwig_from_Andrew/MergedBigwigs/PROTAC_H2AK119Ub.merged.bw")
UNC_K119Ub_CnT <- file.path(CUTnTag_dir, "Calibrated_Bigwig_from_Andrew/MergedBigwigs/UNC_H2AK119Ub.merged.bw")


########################################
### Upload and prepare the PECHi-C data (after CHiCAGO = score of contact)
########################################
Peakmatrix_all <- fread(Peakmatrix_all_file)
dim(Peakmatrix_all) # 1322199      24
head(Peakmatrix_all)

# Keep only the PRC2 perturbation scores
Peakmatrix_score <- Peakmatrix_all[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", 
                                       "dist", "Day5_PRC2_DMSO_total_pooled", "Day5_PRC2_PROTAC_total_pooled", "Day5_PRC2_UNC1999_total_pooled", "day5")]
dim(Peakmatrix_score) # 1322199      15
head(Peakmatrix_score)

### Loop pairs
Peakmatrix_score_DMSO <- Peakmatrix_score[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd", "Day5_PRC2_DMSO_total_pooled")]
Peakmatrix_score_PROTAC <- Peakmatrix_score[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd", "Day5_PRC2_PROTAC_total_pooled")]
Peakmatrix_score_UNC <- Peakmatrix_score[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd", "Day5_PRC2_UNC1999_total_pooled")]

# Round the score column to the nearest integer (or specify digits if needed)
Peakmatrix_score_DMSO[, Day5_PRC2_DMSO_total_pooled := round(Day5_PRC2_DMSO_total_pooled)]
Peakmatrix_score_PROTAC[, Day5_PRC2_PROTAC_total_pooled := round(Day5_PRC2_PROTAC_total_pooled)]
Peakmatrix_score_UNC[, Day5_PRC2_UNC1999_total_pooled := round(Day5_PRC2_UNC1999_total_pooled)]

# Keep only rows where the score is not 0
Peakmatrix_score_DMSO <- Peakmatrix_score_DMSO[Day5_PRC2_DMSO_total_pooled != 0]
Peakmatrix_score_PROTAC <- Peakmatrix_score_PROTAC[Day5_PRC2_PROTAC_total_pooled != 0]
Peakmatrix_score_UNC <- Peakmatrix_score_UNC[Day5_PRC2_UNC1999_total_pooled != 0]

dim(Peakmatrix_score_DMSO) # 375829
dim(Peakmatrix_score_PROTAC) # 395261
dim(Peakmatrix_score_UNC) # 392871

PRC2perturb_loops_pairs_DMSO <- Peakmatrix_score_DMSO[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd")]
PRC2perturb_loops_pairs_DMSO[, baitChr := paste0("chr", baitChr)]
PRC2perturb_loops_pairs_DMSO[, oeChr := paste0("chr", oeChr)]

PRC2perturb_loops_pairs_PROTAC <- Peakmatrix_score_PROTAC[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd")]
PRC2perturb_loops_pairs_PROTAC[, baitChr := paste0("chr", baitChr)]
PRC2perturb_loops_pairs_PROTAC[, oeChr := paste0("chr", oeChr)]

PRC2perturb_loops_pairs_UNC <- Peakmatrix_score_UNC[,c("baitChr", "baitStart", "baitEnd", "oeChr", "oeStart", "oeEnd")]
PRC2perturb_loops_pairs_UNC[, baitChr := paste0("chr", baitChr)]
PRC2perturb_loops_pairs_UNC[, oeChr := paste0("chr", oeChr)]


############################
## Upload baitmap
############################
baitmap <- fread(baitmap_file)
setnames(baitmap, old = colnames(baitmap),
         new = c("chr", "start", "end", "baitID", "PEname"))


##################################################################
##################################################################
#### Plot for bait 128479 #### 
##################################################################
##################################################################
# 12 128267677 128269043 128479

PRC2perturb_loops_pairs_DMSO_filtered_128479 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr12" & baitStart == "128267677" & baitEnd == "128269043"
]

PRC2perturb_loops_pairs_UNC_filtered_128479 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr12" & baitStart == "128267677" & baitEnd == "128269043"
]

PRC2perturb_loops_pairs_PROTAC_filtered_128479 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr12" & baitStart == "128267677" & baitEnd == "128269043"
]

##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 127767677
chrend = 128769043

params <- pgParams(
  chrom = "chr12", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "12", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page
##################################################################
pdf(file.path(working_dir, "Figure5D_left.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_128479,
  params = params, linecolor = "#37a7db", fill="#37a7db",
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_128479,
  params = params, fill = "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_128479,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
bait_of_interest <- baitmap[
  chr == "12" & 
    start == 128267677 & 
    end == 128269043
  ]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 100)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)


dev.off()


##################################################################
##################################################################
#### Plot for bait 157748 #### 
##################################################################
##################################################################
# 14 61280073 61281698 baitID=157748

PRC2perturb_loops_pairs_DMSO_filtered_157748 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr14" & baitStart == "61280073" & baitEnd == "61281698"
]

PRC2perturb_loops_pairs_UNC_filtered_157748 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr14" & baitStart == "61280073" & baitEnd == "61281698"
]

PRC2perturb_loops_pairs_PROTAC_filtered_157748 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr14" & baitStart == "61280073" & baitEnd == "61281698"
]


##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 60780073
chrend = 61781698

params <- pgParams(
  chrom = "chr14", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "14", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page
##################################################################
pdf(file.path(working_dir, "Figure5D_right.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_157748,
  params = params, fill="#37a7db", linecolor = "#37a7db",
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_157748,
  params = params, fill = "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_157748,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
# 14 61280073 61281698 baitID=157748
bait_of_interest <- baitmap[
  chr == "14" & 
    start == 61280073 & 
    end == 61281698
]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 50)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()


##################################################################
##################################################################
#### Plot for bait 43615 #### 
##################################################################
##################################################################
# 1 228007191 228008105 baitID=43615

PRC2perturb_loops_pairs_DMSO_filtered_43615 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr1" & baitStart == "228007191" & baitEnd == "228008105"
]

PRC2perturb_loops_pairs_UNC_filtered_43615 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr1" & baitStart == "228007191" & baitEnd == "228008105"
]

PRC2perturb_loops_pairs_PROTAC_filtered_43615 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr1" & baitStart == "228007191" & baitEnd == "228008105"
]


##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 227507191
chrend = 228508105

params <- pgParams(
  chrom = "chr1", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "1", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page
##################################################################
pdf(file.path(working_dir, "Figure5D_middle.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_43615,
  params = params, fill="#37a7db",linecolor = "#37a7db",
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_43615,
  params = params, fill= "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_43615,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
# 1 228007191 228008105 baitID=43615
bait_of_interest <- baitmap[
  chr == "1" & 
    start == 228007191 & 
    end == 228008105
]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 100)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()


##################################################################
##################################################################
#### Plot for bait 208959 #### 
##################################################################
##################################################################
# 17 39217634 39219111 baitID=208959

PRC2perturb_loops_pairs_DMSO_filtered_208959 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr17" & baitStart == "39217634" & baitEnd == "39219111"
]

PRC2perturb_loops_pairs_UNC_filtered_208959 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr17" & baitStart == "39217634" & baitEnd == "39219111"
]

PRC2perturb_loops_pairs_PROTAC_filtered_208959 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr17" & baitStart == "39217634" & baitEnd == "39219111"
]


##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 38717634
chrend = 39719111

params <- pgParams(
  chrom = "chr17", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "17", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page 
##################################################################
pdf(file.path(working_dir, "SuppFigureS3F_middle.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_208959,
  params = params, fill="#37a7db", linecolor = "#37a7db", 
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_208959,
  params = params, fill = "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_208959,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
# 17 39217634 39219111 baitID=208959

bait_of_interest <- baitmap[
  chr == "17" & 
    start == 39217634 & 
    end == 39219111
]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 50)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()


##################################################################
##################################################################
#### Plot for bait 359979 #### 
##################################################################
##################################################################
# 3 184230037 184231477 baitID=359979

PRC2perturb_loops_pairs_DMSO_filtered_359979 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr3" & baitStart == "184230037" & baitEnd == "184231477"
]

PRC2perturb_loops_pairs_UNC_filtered_359979 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr3" & baitStart == "184230037" & baitEnd == "184231477"
]

PRC2perturb_loops_pairs_PROTAC_filtered_359979 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr3" & baitStart == "184230037" & baitEnd == "184231477"
]


##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 183730037
chrend = 184731477

params <- pgParams(
  chrom = "chr3", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "3", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page
##################################################################
pdf(file.path(working_dir, "SuppFigureS3F_left.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_359979,
  params = params, fill="#37a7db", linecolor = "#37a7db",
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_359979,
  params = params, fill = "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_359979,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
# 3 184230037 184231477 baitID=359979

bait_of_interest <- baitmap[
  chr == "3" & 
    start == 184230037 & 
    end == 184231477
]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 100)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()

##################################################################
##################################################################
#### Plot for bait 434991 #### 
##################################################################
##################################################################
# 5 176365129 176366766 baitID=434991

PRC2perturb_loops_pairs_DMSO_filtered_434991 <- PRC2perturb_loops_pairs_DMSO[
  baitChr == "chr5" & baitStart == "176365129" & baitEnd == "176366766"
]

PRC2perturb_loops_pairs_UNC_filtered_434991 <- PRC2perturb_loops_pairs_UNC[
  baitChr == "chr5" & baitStart == "176365129" & baitEnd == "176366766"
]

PRC2perturb_loops_pairs_PROTAC_filtered_434991 <- PRC2perturb_loops_pairs_PROTAC[
  baitChr == "chr5" & baitStart == "176365129" & baitEnd == "176366766"
]


##################################################################
# 1. Define genomic region and plotting parameters
##################################################################
chrstart = 175865129
chrend = 176866766

params <- pgParams(
  chrom = "chr5", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_Bigwig <- pgParams(
  chrom = "5", 
  chromstart = chrstart, chromend = chrend,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##################################################################
# 2. Create plotgardener page
##################################################################
pdf(file.path(working_dir, "SuppFigureS3F_right.pdf"), width = 6, height = 10)
pageCreate(width = 4, height = 9, default.units = "inches", showGuides = FALSE)

##################################################################
# 3. Plot gene models
##################################################################
genesPlot <- plotGenes(
  params = params, 
  bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 4. Plot interaction arches for each condition
##################################################################
plotPairsArches(
  data = PRC2perturb_loops_pairs_DMSO_filtered_434991,
  params = params, fill="#37a7db", linecolor = "#37a7db",
  y = 1.2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_PROTAC_filtered_434991,
  params = params, fill = "#ef8a62", linecolor = "#ef8a62",
  y = 2, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

plotPairsArches(
  data = PRC2perturb_loops_pairs_UNC_filtered_434991,
  params = params, fill = "#66c2a5", linecolor = "#66c2a5",
  y = 2.8, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches", flip = TRUE,
  alpha = 1
)

##################################################################
# 5. Highlight bait region of interest
##################################################################
# 5 176365129 176366766 baitID=434991

bait_of_interest <- baitmap[
  chr == "5" & 
    start == 176365129 & 
    end == 176366766
]

genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params_Bigwig,
  width = 3, height = 0.2,
  y = 0.7,
  linecolor = "red", fill = "red",
  collapse = TRUE,
  just = c("left", "top"), default.units = "inches"
)

##################################################################
# 6. Plot CUT&Tag
##################################################################
y_range_k27 <- c(0, 50)         # Adjust this to your data range
y_range_k119 <- c(0, 50)

plotSignal(
  data = DMSO_K27me3_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 3.5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 3.4,
  height = 0.5
)

plotSignal(
  data = PROTAC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 4.3,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 4.2,
  height = 0.5
)

plotSignal(
  data = UNC_K27me3_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 5,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k27,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC K27me3 CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 4.9,
  height = 0.5
)

plotSignal(
  data = DMSO_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#37a7db",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#37a7db",
  y = 6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "DMSO H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#37a7db",
  rot = 0, x = 0.6, y = 5.9,
  height = 0.5
)

plotSignal(
  data = PROTAC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#ef8a62",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#ef8a62",
  y = 6.8,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "PROTAC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#ef8a62",
  rot = 0, x = 0.6, y = 6.7,
  height = 0.5
)

plotSignal(
  data = UNC_K119Ub_CnT,
  params = params_Bigwig,
  fill = "#66c2a5",          # Or colors_k27me3[["DMSO"]]
  linecolor = "#66c2a5",
  y = 7.6,                   # Adjust based on your layout
  width = 3, height = 0.5,
  just = c("left", "top"),
  default.units = "inches",
  range = y_range_k119,       # Adjust based on your signal range
  scale = TRUE
)

plotText(
  label = "UNC H2AK119Ub CUT&Tag",
  fontsize = 9, fontcolor ="#66c2a5",
  rot = 0, x = 0.6, y = 7.5,
  height = 0.5
)

## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 8.4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()

