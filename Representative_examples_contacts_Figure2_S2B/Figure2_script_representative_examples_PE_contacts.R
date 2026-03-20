################################################################################
### 2- Script for Figure 2 and supplemental Figure S2B - Representative examples of PE contacts
### Author: Dr. Marina Nocente
###
### This script reproduces Figures 2A panel and S2B showing:
### - Figure 2A: example of PE-CHiC significant contacts and CUT&Tag signal
### - Figures 2C and S2B: representative examples of PE contacts during the naive-to-primed transition
################################################################################

########################################
### Load libraries
########################################
library(data.table)
library(dplyr)
library(tidyr)

library(plotgardener)
library(grid)

########################################
### User configuration
########################################
BASE_DIR <- "~/Documents/Bioinformatics"
FINAL_DIR <- file.path(BASE_DIR, "Final_scripts_paper")
FIG_DIR <- file.path(BASE_DIR, "Find_examples")
CUTnTag_DIR <- file.path(BASE_DIR, "CnT_transition/bigwig_merged_just_sorted_indexedBamFiles")
setwd(FINAL_DIR)

dtw_clusters_file <- file.path(FINAL_DIR, "PoisedEnhancers/Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt")
Peakmatrix_raw_score_file <- file.path(BASE_DIR, "PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
baitmap_file <- file.path(BASE_DIR, "CnT_PRC2_perturb_Day5/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt")
tss_file <- file.path(BASE_DIR, "GO_OtherEnds/OtherEnd_TSS_coordinates.txt")

########################################
########################################
### Figure 2A: example of PE-CHiC significant contacts and CUT&Tag signal
########################################
########################################

########################################
### Load the full matrix of PE-CHiC contact scores (raw scores) and the temporal clusters of contact
### Merge the peak matrix with the temporal clusters
########################################
Peakmatrix_all <- fread(Peakmatrix_raw_score_file)[,c("baitID", "oeID", "naive", "day1", "day3", "day5", "day7", "day10", "day14", "baitChr", "baitStart", "baitEnd", "baitName", "oeChr", "oeStart", "oeEnd", "oeName")]
dim(Peakmatrix_all) # 1322199      17

dtw_clusters <- fread(dtw_clusters_file)[, .(baitID, oeID, merged_dtw_trans_cl6_2with6)]
dim(dtw_clusters) # 9998    3

peakMatrix <- merge(Peakmatrix_all, dtw_clusters, by = c("baitID", "oeID"))
dim(peakMatrix) # 9998   18

########################################
######  Plot using plotGardener
########################################
## Keep just this interaction: CXCL16 gene
# baitID = 201999 - 17 4787398 4787760
# oeID = 201982 - 17   4731783 4736881

### Set region parameters for the plot
region_chr <- "chr17"
region_start <- 4730000
region_end <- 4790000

########## Set parameters main and text
params_main <- pgParams(chrom = region_chr, chromstart = region_start, chromend = region_end, 
                        assembly = "hg38")

params_text <- pgParams(chrom = region_chr, chromstart = region_start, chromend = region_end, 
                        assembly = "hg38")

########## Selection of the contact of interest
peakMatrix_filtered <- peakMatrix[baitStart == 4787398 & oeStart == 4731783]
head(peakMatrix_filtered)

baitmap <- fread(baitmap_file)
setnames(baitmap, c("V1", "V2", "V3", "V4"), c("chrom", "start", "end", "baitID"))
baitmap[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]

bait_of_interest <- baitmap[
  chrom == "chr17" & 
    start == 4787398 & 
    end == 4787760
]

######## Interaction and score for day14
## Start from the original data.table
hic_df_day14 <- peakMatrix_filtered[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14
)]

## Define genomic region and widths of all plots
params <- pgParams(
  chrom = "chr17", 
  chromstart = 4700000, chromend = 4800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "17", 
  chromstart = 4700000, chromend = 4800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Create plotgardener page
pdf(file.path(FIG_DIR, "Example_Figure1_forSchematic.pdf"))
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Plot genes with background color highlighting size
genesPlot <- plotGenes(
  params = params, bg = "#f6f6f6",
  y = 0.5, width = 3, height = 0.5,
  just = c("left", "top"), default.units = "inches"
) #bg = "#f6f6f6"


# plot bait
genesBait <- plotRanges( 
  data = bait_of_interest,
  params = params,
  width = 3, 
  linecolor="red", 
  fill = "red", 
  collapse = T,
  y=0.7,
  height = 0.2
)

## Plot and align DNA loops in region
## day14
bedpeLoops_day14 <- plotPairsArches(
  data = hic_df_day14,
  params = params, fill = "green4", linecolor = "darkgreen",
  y = 1.2, width = 3, height = hic_df_day14$score/10,
  just = c("left", "top"), default.units = "inches", flip = TRUE
)

plotText(
  label = "PECHiC", 
  fontsize = 9, fontcolor = "darkblue", rot = 0, x=1.3, y = 1.1, height = 0.5,
)


## CUT&Tag signals bigwig
y_range_k27 <- c(0, 400)
y_range_k4 <- c(0, 100)
y_range_k4me3 <- c(0, 800)

signal_K27me3 <- plotSignal(
  data = file.path(CUTnTag_DIR, "merged_Day14_273.bam.bw"),
  params = params_bigwig, fill = "darkred", linecolor = "darkred",
  y = 1.9,
  width = 3,
  height = 0.5, 
  just = c("left", "top"), default.units = "inches",
  range = y_range_k27, scale = TRUE
)

plotText(
  label = "H3K27me3", 
  fontsize = 9, fontcolor = "darkred", rot = 0, x=0.8, y = 1.8, height = 0.5,
) 

signal_K4me1 <- plotSignal(
  data = file.path(CUTnTag_DIR, "merged_Day14_41.bam.bw"),
  params = params_bigwig, fill = "blue", linecolor = "blue",
  y = 2.6,
  width = 3,
  height = 0.5, 
  just = c("left", "top"), default.units = "inches",
  range = y_range_k4, scale = TRUE
)

plotText(
  label = "H3K4me1", 
  fontsize = 9, fontcolor = "blue", rot = 0, x=0.8, y = 2.5, height = 0.5,
) 

signal_K4me3 <- plotSignal(
  data = file.path(CUTnTag_DIR, "merged_Day14_43.bam.bw"),
  params = params_bigwig, fill = "orange", linecolor = "orange",
  y = 3.3,
  width = 3,
  height = 0.5, 
  just = c("left", "top"), default.units = "inches",
  range = y_range_k4me3, scale = TRUE
)

plotText(
  label = "H3K4me3", 
  fontsize = 9, fontcolor = "orange", rot = 0, x=0.8, y = 3.2, height = 0.5,
) 


## Annotate genome label
annoGenomeLabel(
  plot = genesPlot, 
  params = params, scale = "Mb",
  y = 4,
  just = c("left", "top")
)

## Labels
plotText(
  label = "Genes", 
  fontsize = 8, fontcolor = "black", rot = 0, x=0.2, y = 0.45, height = 0.5,
)

dev.off()

########################################
########################################
### Figures 2C and S2B : representative examples of PE contacts during the naive-to-primed transition
########################################
########################################

########################################
## Keep only TSS Other Ends
########################################
## How many other ends are poised enhancers?
sum(grepl("poised enhancer", peakMatrix$oeName, ignore.case = TRUE)) # 1626

## How many other ends are "other things" = TSS?
sum(!grepl("poised enhancer", peakMatrix$oeName, ignore.case = TRUE)) # 8372

peakMatrix_oeTSS <- peakMatrix[!grepl("poised enhancer", oeName, ignore.case = TRUE),]
dim(peakMatrix_oeTSS) # 8372   18

# Upload coordinates of Other Ends which are TSS (with gene_id and gene_names)
TssOE <- fread(tss_file)
TssOE <- TssOE[,c("oeID", "baitID", "TSSchr", "TSSstart", "TSSend", "TSSstrand", "gene_name", "gene_id", "chr", "start", "end")] 
dim(TssOE) # 7538   11

peakMatrix_oeTSS <- merge(peakMatrix_oeTSS, TssOE, by=c("baitID", "oeID"))
dim(peakMatrix_oeTSS) # 7486   27

## Keep only the row which have a temporal cluster:
peakMatrix_tempClusters_NtoP_oeTSS <- peakMatrix_oeTSS[merged_dtw_trans_cl6_2with6 %in% c("Constant", "Lost", "Transient", "Gained early", "Gained late")]
dim(peakMatrix_tempClusters_NtoP_oeTSS) # 7486   27


########################################
## Selection of baitID (poised enhancers) involved in different types of contact
########################################
## Example of Constant bait
bait_of_interest_constant <- baitmap[
  chrom == "chr17" &
    start == 4787398 &
    end == 4787760
]

## Example of Gained early bait
bait_of_interest_gainedEarly <- baitmap[
  chrom == "chr11" & 
    start == 124863652 & 
    end == 124864276
]

## Example of Gained late bait
bait_of_interest_gainedLate <- baitmap[
  chrom == "chr17" & 
    start == 48720806 & 
    end == 48720955
]

## Example of Lost bait
bait_of_interest_lost <- baitmap[
  chrom == "chr11" & 
    start == 67318611 & 
    end == 67318982
]

## Example of Transient bait
bait_of_interest_transient <- baitmap[
  chrom == "chr16" & 
    start == 22376217 & 
    end == 22377253
]

########################################
## Selection of the contacts of interest
########################################
## Constant contact with the TSS of CXCL16 gene
peakMatrix_clusters_NtoP_filtered_constant <- peakMatrix_tempClusters_NtoP_oeTSS[baitStart == 4787398 & oeStart == 4731783]
head(peakMatrix_clusters_NtoP_filtered_constant)

## Gained early contact with the TSS of MSANTD2 gene
peakMatrix_clusters_NtoP_filtered_gainedEarly <- peakMatrix_tempClusters_NtoP_oeTSS[baitStart == 124863652 & oeStart == 124763235]
head(peakMatrix_clusters_NtoP_filtered_gainedEarly)

## Gained late contact with the TSS of HOXB8 gene
peakMatrix_clusters_NtoP_filtered_gainedLate <- peakMatrix_tempClusters_NtoP_oeTSS[baitStart == 48720806 & oeStart == 48609861]
head(peakMatrix_clusters_NtoP_filtered_gainedLate)

## Transient contact with the TSS of EEF2K gene
peakMatrix_tempClusters_NtoP_filtered_transient <- peakMatrix_tempClusters_NtoP_oeTSS[baitStart == 22376217 & oeStart == 22202138]
head(peakMatrix_tempClusters_NtoP_filtered_transient)

## Lost contact with the TSS of GRK2 gene
peakMatrix_tempClusters_NtoP_filtered_lost <- peakMatrix_tempClusters_NtoP_oeTSS[baitStart == 67318611 & oeStart == 67263137]
head(peakMatrix_tempClusters_NtoP_filtered_lost)

########################################
## Interaction and score for each time point and contact type
########################################
## Constant contact
hic_df_naive_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  # Change to "naive" timepoint column
)]

hic_df_day1_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  # Change to "day1" timepoint column
)]

hic_df_day3_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  # Change to "day3" timepoint column
)]

hic_df_day5_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  # Change to "day5" timepoint column
)]

hic_df_day7_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7  # Change to "day7" timepoint column
)]

hic_df_day10_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10  # Change to "day10" timepoint column
)]

hic_df_day14_constant <- peakMatrix_clusters_NtoP_filtered_constant[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]


## Gained early contact
hic_df_naive_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  # Change to "naive" timepoint column
)]

hic_df_day1_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  # Change to "day1" timepoint column
)]

hic_df_day3_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  # Change to "day3" timepoint column
)]

hic_df_day5_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  # Change to "day5" timepoint column
)]

hic_df_day7_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7  # Change to "day7" timepoint column
)]

hic_df_day10_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10  # Change to "day10" timepoint column
)]

hic_df_day14_gainedEarly <- peakMatrix_clusters_NtoP_filtered_gainedEarly[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]


## Gained late contact
hic_df_naive_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  # Change to "naive" timepoint column
)]

hic_df_day1_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  # Change to "day1" timepoint column
)]

hic_df_day3_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  # Change to "day3" timepoint column
)]

hic_df_day5_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  # Change to "day5" timepoint column
)]

hic_df_day7_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7  # Change to "day7" timepoint column
)]

hic_df_day10_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10  # Change to "day10" timepoint column
)]

hic_df_day14_gainedLate <- peakMatrix_clusters_NtoP_filtered_gainedLate[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]


## Transient contact
hic_df_naive_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  # Change to "naive" timepoint column
)]

hic_df_day1_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  # Change to "day1" timepoint column
)]

hic_df_day3_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  # Change to "day3" timepoint column
)]

hic_df_day5_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  # Change to "day5" timepoint column
)]

hic_df_day7_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7  # Change to "day7" timepoint column
)]

hic_df_day10_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10  # Change to "day10" timepoint column
)]

hic_df_day14_Transient <- peakMatrix_tempClusters_NtoP_filtered_transient[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]


## Lost contact
hic_df_naive_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  # Change to "naive" timepoint column
)]

hic_df_day1_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  # Change to "day1" timepoint column
)]

hic_df_day3_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  # Change to "day3" timepoint column
)]

hic_df_day5_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  # Change to "day5" timepoint column
)]

hic_df_day7_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7  # Change to "day7" timepoint column
)]

hic_df_day10_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10  # Change to "day10" timepoint column
)]

hic_df_day14_lost <- peakMatrix_tempClusters_NtoP_filtered_lost[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]


######################################################################################################
### Plot of arcsinh-transformed Chicago scores with color scale (0 → max arcsinh score) for CONSTANT
######################################################################################################
timepoints_labels  <- c("Naive","Day1","Day3","Day5","Day7","Day10","Day14")
timepoints_columns <- c("naive","day1","day3","day5","day7","day10","day14")

#######################################
## General plot parameters
#######################################
params <- pgParams(
  chrom = "chr17",
  chromstart = 4700000, chromend = 4800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#######################################
## Build Hi-C list (Constant)
#######################################
datasets <- list(
  hic_df_naive_constant, hic_df_day1_constant, hic_df_day3_constant,
  hic_df_day5_constant, hic_df_day7_constant, hic_df_day10_constant,
  hic_df_day14_constant
)

names(datasets) <- timepoints_labels

#######################################
## Color: arcsinh + global scale
#######################################
score_palette <- c("white", "darkgreen")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$score_asinh <- asinh(df$score)
  datasets[[i]] <- df
}

# global max
all_vals <- unlist(lapply(datasets, function(x) x$score_asinh))
#  Naive     Day1     Day3     Day5     Day7    Day10    Day14 
# 1.537948 1.117696 1.680639 1.629600 1.729936 1.440054 1.835236 
global_max_asinh <- 2.86 # the five examples have the same maximum

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  datasets[[i]] <- df
}

#######################################
## Compute arcsinh scores and global maximum
#######################################
pdf(file.path(FIG_DIR, "Figure2C_Constant_example_PEcontact.pdf"), width = 8, height = 7)
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_constant, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
timepoint_spacing <- 0.5

#######################################
## Loop: all timepoints
#######################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 0.7 + (i-1)*timepoint_spacing
  
  ## ARCS
  df_arc <- datasets[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkgreen"))(256 )
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + 0.5
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.2,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=4.6, just = c("left","top")
)

## label
plotText(
  label = "constant (arcsinh Chicago score)",
  fontsize = 10,
  fontcolor = "darkgreen",
  rot = 0,
  x = 2.0, y = 0.2
)

plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkgreen"))(256) , alpha.f = 1) # currently no transparency as in the default,

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(5.5 * sf, "inches"),
  y = unit(3 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((3-5/2) * sf, (3+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((5.75 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score\n(0 – ", round(global_max_asinh,2), ")"),
  x = unit(6 * sf, "inches"),
  y = unit(3 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)
dev.off()

######################################################################################################
### Plot arcsinh-transformed Chicago scores with color scale (0 → max arcsinh score) for GAINED EARLY
######################################################################################################

#######################################
## General plot parameters
#######################################
params <- pgParams(
  chrom = "chr11", 
  chromstart = 124700000, chromend = 124900000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#######################################
## Dataset list (gained early)
#######################################
datasets <- list(
  hic_df_naive_gainedEarly, hic_df_day1_gainedEarly, hic_df_day3_gainedEarly,
  hic_df_day5_gainedEarly, hic_df_day7_gainedEarly, hic_df_day10_gainedEarly,
  hic_df_day14_gainedEarly
)

names(datasets) <- timepoints_labels

#######################################
## Color: arcsinh + global scale
#######################################
score_palette <- c("white", "darkred")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$score_asinh <- asinh(df$score)
  datasets[[i]] <- df
}

# global max
all_vals <- unlist(lapply(datasets, function(x) x$score_asinh))
#  Naive     Day1     Day3     Day5     Day7    Day10    Day14 
# 0.0000000 1.3719229 0.7960238 1.8718955 1.4698352 1.4918255 1.4468114
global_max_asinh <- 2.86 # the five examples have the same maximum

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  datasets[[i]] <- df
}

#######################################
## Compute arcsinh scores and global maximum
#######################################

pdf(file.path(FIG_DIR, "Figure2C_gainedEarly_example_PEcontact.pdf"), width = 8, height = 7)
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_gainedEarly, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
timepoint_spacing <- 0.5

#######################################
## Loop: all timepoints
#######################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 0.7 + (i-1)*timepoint_spacing
  
  ## ARCS
  df_arc <- datasets[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkred"))(256 )
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + 0.5
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.2,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=4.6, just = c("left","top")
)

## label
plotText(
  label = "gained early (arcsinh Chicago score)",
  fontsize = 10,
  fontcolor = "darkred",
  rot = 0,
  x = 2.0, y = 0.2
)

plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkred"))(256) , alpha.f = 1) # currently no transparency as in the default,

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(5.5 * sf, "inches"),
  y = unit(3 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((3-5/2) * sf, (3+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((5.75 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score\n(0 – ", round(global_max_asinh,2), ")"),
  x = unit(6 * sf, "inches"),
  y = unit(3 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)
dev.off()

######################################################################################################
### Plot arcsinh-transformed Chicago scores with color scale (0 → max arcsinh score) for GAINED LATE
######################################################################################################

#######################################
## General plot parameters
#######################################
params <- pgParams(
  chrom = "chr17",
  chromstart = 48580000, chromend = 48730000, 
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#######################################
## Color mapping function
#######################################
score_palette <- c("white", "darkorange")
color_fun <- colorRamp(score_palette, space = "Lab")

getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0          # NA becomes white
  scaled <- scores / max_val          # scale to 0–1 range
  scaled[scaled > 1] <- 1             # safety clamp
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255)
}

#######################################
## Dataset list (gained early)
#######################################
datasets <- list(
  hic_df_naive_gainedLate, hic_df_day1_gainedLate, hic_df_day3_gainedLate,
  hic_df_day5_gainedLate, hic_df_day7_gainedLate, hic_df_day10_gainedLate,
  hic_df_day14_gainedLate
)

names(datasets) <- timepoints_labels

#######################################
## Color: arcsinh + global scale
#######################################
score_palette <- c("white", "darkorange")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$score_asinh <- asinh(df$score)
  datasets[[i]] <- df
}

# global max
all_vals <- unlist(lapply(datasets, function(x) x$score_asinh))
#  Naive     Day1     Day3     Day5     Day7    Day10    Day14 
# 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.9610447 2.3816986 
global_max_asinh <- 2.86 # the five examples have the same maximum

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  datasets[[i]] <- df
}

#######################################
## Compute arcsinh scores and global maximum
#######################################

pdf(file.path(FIG_DIR, "Figure2C_gainedLate_example_PEcontact.pdf"), width = 8, height = 7)
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_gainedLate, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
timepoint_spacing <- 0.5

#######################################
## Loop: all timepoints
#######################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 0.7 + (i-1)*timepoint_spacing
  
  ## ARCS
  df_arc <- datasets[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkorange"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + 0.5
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.2,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=4.6, just = c("left","top")
)

## label
plotText(
  label = "gained late (arcsinh Chicago score)",
  fontsize = 10,
  fontcolor = "darkorange",
  rot = 0,
  x = 2.0, y = 0.2
)

plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkorange"))(256) , alpha.f = 1) # currently no transparency as in the default,

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(5.5 * sf, "inches"),
  y = unit(3 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((3-5/2) * sf, (3+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((5.75 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score\n(0 – ", round(global_max_asinh,2), ")"),
  x = unit(6 * sf, "inches"),
  y = unit(3 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)
dev.off()

######################################################################################################
### Plot arcsinh-transformed Chicago scores with color scale (0 → max arcsinh score) for TRANSIENT
######################################################################################################

#######################################
## General plot parameters
#######################################
params <- pgParams(
  chrom = "chr16", 
  chromstart = 22190000, chromend = 22400000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#######################################
## Color mapping function
#######################################
score_palette <- c("white", "purple")
color_fun <- colorRamp(score_palette, space = "Lab")

getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0          # NA becomes white
  scaled <- scores / max_val          # scale to 0–1 range
  scaled[scaled > 1] <- 1             # safety clamp
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255)
}

#######################################
## Dataset list (transient)
#######################################
datasets <- list(
  hic_df_naive_Transient, hic_df_day1_Transient, hic_df_day3_Transient,
  hic_df_day5_Transient, hic_df_day7_Transient, hic_df_day10_Transient,
  hic_df_day14_Transient
)

names(datasets) <- timepoints_labels

#######################################
## Color: arcsinh + global scale
#######################################
score_palette <- c("white", "purple")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$score_asinh <- asinh(df$score)
  datasets[[i]] <- df
}

# global max
all_vals <- unlist(lapply(datasets, function(x) x$score_asinh))
#  Naive     Day1     Day3     Day5     Day7    Day10    Day14 
# 0.0000000 0.8859033 2.3866830 1.8302047 0.8671292 0.8709328 0.0000000 
global_max_asinh <- 2.86 # the five examples have the same maximum

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  datasets[[i]] <- df
}

#######################################
## Compute arcsinh scores and global maximum
#######################################

pdf(file.path(FIG_DIR, "SuppFigureS2B_Transient_example_PEcontact.pdf"), width = 8, height = 7)
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_transient, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
timepoint_spacing <- 0.5

#######################################
## Loop: all timepoints
#######################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 0.7 + (i-1)*timepoint_spacing
  
  ## ARCS
  df_arc <- datasets[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","purple"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + 0.5
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.2,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=4.6, just = c("left","top")
)

## label
plotText(
  label = "Transient (arcsinh Chicago score)",
  fontsize = 10,
  fontcolor = "purple",
  rot = 0,
  x = 2.0, y = 0.2
)

plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "purple"))(256) , alpha.f = 1) # currently no transparency as in the default,

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(5.5 * sf, "inches"),
  y = unit(3 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((3-5/2) * sf, (3+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((5.75 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score\n(0 – ", round(global_max_asinh,2), ")"),
  x = unit(6 * sf, "inches"),
  y = unit(3 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)
dev.off()

######################################################################################################
### Plot arcsinh-transformed Chicago scores with color scale (0 → max arcsinh score) for LOST
######################################################################################################

#######################################
## General plot parameters
#######################################
params <- pgParams(
  chrom = "chr11", 
  chromstart = 67250000, chromend = 67350000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#######################################
## Color mapping function
#######################################
score_palette <- c("white", "#009ACD")
color_fun <- colorRamp(score_palette, space = "Lab")

getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0          # NA becomes white
  scaled <- scores / max_val          # scale to 0–1 range
  scaled[scaled > 1] <- 1             # safety clamp
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255)
}

#######################################
## Dataset list (gained early)
#######################################
datasets <- list(
  hic_df_naive_lost, hic_df_day1_lost, hic_df_day3_lost,
  hic_df_day5_lost, hic_df_day7_lost, hic_df_day10_lost,
  hic_df_day14_lost
)

names(datasets) <- timepoints_labels

#######################################
## Color: arcsinh + global scale
#######################################
score_palette <- c("white", "#009ACD")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$score_asinh <- asinh(df$score)
  datasets[[i]] <- df
}

# global max
all_vals <- unlist(lapply(datasets, function(x) x$score_asinh))
#  Naive     Day1     Day3     Day5     Day7    Day10    Day14 
# 2.8540514 2.0992570 1.9126000 2.1447315 1.7103808 0.6513029 0.0000000  
global_max_asinh <- 2.86 # the five examples have the same maximum

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(datasets)) {
  df <- datasets[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  datasets[[i]] <- df
}

#######################################
## Compute arcsinh scores and global maximum
#######################################

pdf(file.path(FIG_DIR, "SuppFigureS2B_Lost_example_PEcontact.pdf"), width = 8, height = 7)
pageCreate(width = 4, height = 5, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_lost, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
timepoint_spacing <- 0.5

#######################################
## Loop: all timepoints
#######################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 0.7 + (i-1)*timepoint_spacing
  
  ## ARCS
  df_arc <- datasets[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","#009ACD"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + 0.5
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.2,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=4.6, just = c("left","top")
)

## label
plotText(
  label = "Lost (arcsinh Chicago score)",
  fontsize = 10,
  fontcolor = "#009ACD",
  rot = 0,
  x = 2.0, y = 0.2
)

plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "#009ACD"))(256) , alpha.f = 1) # currently no transparency as in the default,

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(5.5 * sf, "inches"),
  y = unit(3 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((3-5/2) * sf, (3+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((5.75 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score\n(0 – ", round(global_max_asinh,2), ")"),
  x = unit(6 * sf, "inches"),
  y = unit(3 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)
dev.off()

