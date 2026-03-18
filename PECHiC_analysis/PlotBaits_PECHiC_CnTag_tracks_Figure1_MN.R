################################################################################
### 1- PlotBait from CHiCAGO and contact and CUT&Tag tracks for Figure 1
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 1 panels showing:
### - CHiCAGO PlotBaits
### - PE-CHiC significant contacts (score >= 5)
### - CUT&Tag tracks (H3K27me3, H3K4me1, H3K4me3)
################################################################################

########################################
### Load libraries
########################################
library(data.table)    # Fast data handling
library(plotgardener)  # Genome visualization
library(tidyverse)     # Data manipulation (select, rename, etc.)

library(Chicago)
library(argparser)
setDTthreads(4)

########################################
### User configuration
########################################
## Base directories (edit these to match your system)
## PECHiC refers to Poised Enhancer Capture Hi-C (PE-CHiC)
BASE_DIR <- "~/Documents/Bioinformatics"
PECHIC_DIR <- "/rds/general/user/mnocente/home/analysis/PECHiC"
CNT_DIR <- file.path(BASE_DIR, "CnT_transition")

## Output directory (edit these to match your system)
OUT_DIR <- file.path(BASE_DIR, "Find_examples")
setwd(OUT_DIR)

## CHiCAGO RDS files
# Files obtained after CHiCAGO analysis of the poised enhancer Capture Hi-C data
CHICAGO_FILES <- list(
  naive  = file.path(PECHIC_DIR, "Rds/naive_comb_noDown.Rds"),
  primed = file.path(PECHIC_DIR, "Rds/primed_comb_noDown.Rds")
)

## PeakMatrix file (PECHi-C scores)
PEAKMATRIX_FILE <- file.path(
  BASE_DIR,
  "PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt"
)

## Baitmap file
BAITMAP_FILE <- file.path(
  BASE_DIR,
  "CnT_PRC2_perturb_Day5/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt"
)

## CUT&Tag bigWig files
BIGWIGS <- list(
  naive = list(
    K27me3 = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Naive_273.bam.bw"),
    K4me1  = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Naive_41.bam.bw"),
    K4me3  = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Naive_43.bam.bw")
  ),
  primed = list(
    K27me3 = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Primed_273.bam.bw"),
    K4me1  = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Primed_41.bam.bw"),
    K4me3  = file.path(CNT_DIR, "bigwig_merged_just_sorted_indexedBamFiles/merged_Primed_43.bam.bw")
  )
)


########################################
### CHiCAGO PlotBaits
########################################
## Using the plotBaits() function from CHiCAGO tool
## Done on a cluster
## Load CHiCAGO interaction objects
interactions <- lapply(CHICAGO_FILES, readRDS)

## Bait IDs used for Figure 1
baits <- list(
  bait580854 = 580854,
  bait68054 = 68054,
  bait1340 = 1340
)

## Generate PlotBaits PDFs for each bait and condition
for (bait_name in names(baits)) {
  baitID <- baits[[bait_name]]
  
  for (condition in names(interactions)) {
    pdf_filename <- sprintf("PlotBait_%s_PECHiC_%s_500.pdf", bait_name, condition)
    pdf(pdf_filename)
    plotBaits(interactions[[condition]], baits = baitID, plotBprof= TRUE, maxD = 500000, ylim = c(0,30))
    dev.off()
  }
}

########################################
### Load and process PE-CHiC PeakMatrix
########################################
## Load peak matrix
Peakmatrix_all <- fread(PEAKMATRIX_FILE)
dim(Peakmatrix_all) # 1322199      24

## Extract naive and primed contact scores
Peakmatrix_naive <- Peakmatrix_all[, c(1:11,15)]
Peakmatrix_primed <- Peakmatrix_all[, c(1:11,22)]

## Keep only bait–OE coordinates and contact score
Peakmatrix_naive <- Peakmatrix_naive[,c(1:3,6:8,12)]
Peakmatrix_primed <- Peakmatrix_primed[,c(1:3,6:8,12)]

## Remove zero-score interactions
Peakmatrix_naive <- Peakmatrix_naive[naive != 0]
Peakmatrix_primed <- Peakmatrix_primed[primed != 0]

## Ensure chromosome names are prefixed with 'chr'
Peakmatrix_naive[,  baitChr := paste0("chr", baitChr)]
Peakmatrix_naive[,  oeChr   := paste0("chr", oeChr)]
Peakmatrix_primed[, baitChr := paste0("chr", baitChr)]
Peakmatrix_primed[, oeChr   := paste0("chr", oeChr)]

## Keep only significant CHiCAGO interactions (score >= 5)
Peakmatrix_naive_Chicago5 <- Peakmatrix_naive[naive >= 5]
Peakmatrix_primed_Chicago5 <- Peakmatrix_primed[primed >= 5]

dim(Peakmatrix_naive_Chicago5) # 48038     7
dim(Peakmatrix_primed_Chicago5) # 62210     7


########################################
### Load baitmap
########################################

baitmap <- fread(BAITMAP_FILE)
setnames(baitmap, c("V1", "V2", "V3", "V4"),
         c("chrom", "start", "end", "baitID"))

## Ensure consistent chromosome naming
baitmap[, chrom := ifelse(startsWith(chrom, "chr"),
                          chrom,
                          paste0("chr", chrom))]


################################################################################
### Representation of contacts and CUT&Tag tracks for example baits
################################################################################

###############################################################################
## IMPORTANT NOTE ON CHROMOSOME NAMING
##
## plotgardener::plotSignal() expects chromosome names WITHOUT the "chr" prefix
## for bigWig files. Therefore:
## - params        -> uses "chrN"
## - params_bigwig -> uses "N"
###############################################################################


########################################
### Bait 580854 (chr10:22,324,300–22,325,358)
########################################

bait_of_interest <- baitmap[
  chrom == "chr10" & 
    start == 22324300 & 
    end == 22325358
]

peak_naive <- Peakmatrix_naive_Chicago5[baitStart == 22324300 & baitEnd == 22325358]
peak_primed <- Peakmatrix_primed_Chicago5[baitStart == 22324300 & baitEnd == 22325358]

## Genomic window for plotting
params <- pgParams(
  chrom = "chr10", 
  chromstart = 21824300, chromend = 22825358,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "10", 
  chromstart = 21824300, chromend = 22825358,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Signal ranges (fixed to ensure comparability)
y_range <- list(
  K27me3 = c(0, 500),
  K4me1  = c(0, 90),
  K4me3  = c(0, 400)
)

########################################
### NAIVE
########################################

pdf("Figure1_plotBaits580854_NAIVE_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

genesPlot <- plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)


plotPairsArches(
  data = peak_naive,
  params = params, fill = "darkred", linecolor = "darkred",
  y = 1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)


signal_K27me3 <- plotSignal(
  data = BIGWIGS$naive$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$naive$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$naive$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Naive", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()


########################################
### PRIMED
########################################

pdf("Figure1_plotBaits580854_PRIMED_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

plotPairsArches(
  data = peak_primed,
  params = params, fill = "darkred", linecolor = "darkred",
  y =1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)

signal_K27me3 <- plotSignal(
  data = BIGWIGS$primed$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$primed$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$primed$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Primed", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()


########################################
### Bait 68054 (chr10:101,134,702–101,135,992)
########################################

bait_of_interest <- baitmap[
  chrom == "chr10" & 
    start == 101134702 & 
    end == 101135992
]

peak_naive <- Peakmatrix_naive_Chicago5[baitStart == 101134702 & baitEnd == 101135992]
peak_primed <- Peakmatrix_primed_Chicago5[baitStart == 101134702 & baitEnd == 101135992]


## Genomic window for plotting
params <- pgParams(
  chrom = "chr10", 
  chromstart = 100634702, chromend = 101635992,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Define genomic region and widths of all plots
params_bigwig <- pgParams(
  chrom = "10", 
  chromstart = 100634702, chromend = 101635992,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Signal ranges
y_range <- list(
  K27me3 = c(0, 1500),
  K4me1  = c(0, 150),
  K4me3  = c(0, 400)
)


########################################
### NAIVE
########################################

pdf("Figure1_plotBaits68054_NAIVE_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

plotPairsArches(
  data = peak_naive,
  params = params, fill = "darkred", linecolor = "darkred",
  y = 1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)

signal_K27me3 <- plotSignal(
  data = BIGWIGS$naive$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$naive$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$naive$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Naive", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()


########################################
### PRIMED
########################################

pdf("Figure1_plotBaits68054_PRIMED_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

plotPairsArches(
  data = peak_primed,
  params = params, fill = "darkred", linecolor = "darkred",
  y = 1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)

signal_K27me3 <- plotSignal(
  data = BIGWIGS$primed$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$primed$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$primed$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Primed", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()


########################################
### Bait 1340 (chr1:6,246,497–6,247,425)
########################################

bait_of_interest <- baitmap[
  chrom == "chr1" & 
    start == 6246497 & 
    end == 6247425
]

peak_naive <- Peakmatrix_naive_Chicago5[baitStart == 6246497 & baitEnd == 6247425]
peak_primed <- Peakmatrix_primed_Chicago5[baitStart == 6246497 & baitEnd == 6247425]

## Genomic window for plotting
params <- pgParams(
  chrom = "chr1", 
  chromstart = 5746497, chromend = 6747425,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Define genomic region and widths of all plots
params_bigwig <- pgParams(
  chrom = "1", 
  chromstart = 5746497, chromend = 6747425,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Signal ranges
y_range <- list(
  K27me3 = c(0, 400),
  K4me1  = c(0, 200),
  K4me3  = c(0, 500)
)


########################################
### NAIVE
########################################

pdf("Figure1_plotBaits1340_NAIVE_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

plotPairsArches(
  data = peak_naive,
  params = params, fill = "darkred", linecolor = "darkred",
  y = 1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)

signal_K27me3 <- plotSignal(
  data = BIGWIGS$naive$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$naive$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$naive$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Naive", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()


########################################
### PRIMED
########################################

pdf("Figure1_plotBaits1340_PRIMED_CUTTag_ARCS.pdf", width = 8, height = 7)

pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.8,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"

plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

plotPairsArches(
  data = peak_primed,
  params = params, fill = "darkred", linecolor = "darkred",
  y = 1.5, width = 3, height = 0.6, alpha = 1,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "black", rot = 0, x=0.5, y = 1.4, height = 0.5,
)

signal_K27me3 <- plotSignal(
  data = BIGWIGS$primed$K27me3,
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 2.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K27me3, scale = TRUE
)

plotText(
  label = "H3K27me3 signal", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 2.4, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = BIGWIGS$primed$K4me1,
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 3.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me1, scale = TRUE
)

plotText(
  label = "H3K4me1 signal", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 3.4, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = BIGWIGS$primed$K4me3,
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 4.5,
  width = 3,
  height = 0.8, 
  just = c("left", "top"), default.units = "inches",
  range = y_range$K4me3, scale = TRUE
)

plotText(
  label = "H3K4me3 signal", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 4.4, height = 0.5,
) 

annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Kb",
  y = 5.5,
  just = c("left", "top")
)

plotText(
  label = "Primed", 
  fontsize = 9, fontcolor = "Black", rot = 0, x=2, y = 0.3, height = 0.5,
) 

dev.off()
