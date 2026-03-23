################################################################################
### Dynamics of enhancer poising during the naive-to-primed transition - Figures 3C and S2C
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 3C and Figure S2C showing:
###  H3K27me3 and H3K4me1 CUT&Tag and PE-Capture Hi-C contacts dynamics across the transition for different classes of contact
###  Example contacts: constant, gained early, gained late
################################################################################

########################################
### Load required libraries
########################################
library(data.table)
library(dplyr)
library(tidyr)

library(plotgardener)
library(grid)

########################################
### User configuration
########################################
root_dir <- "~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT"
CUTnTag_file_dir <- "~/Documents/Bioinformatics/CnT_transition/bigwig_merged_just_sorted_indexedBamFiles/"

## Custom functions
source(file.path(root_dir, "varistran.R"))

## Input files
Peakmatrix_raw_score_file <- file.path("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
dtw_clusters_file <- file.path("~/Documents/Bioinformatics/Final_scripts_paper/PoisedEnhancers/Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt")
baitmap_file <- file.path("~/Documents/Bioinformatics/CnT_PRC2_perturb_Day5/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt")

############################################################
############################################################
### Load the full matrix of PE-CHiC contact scores (raw scores) and the temporal clusters of contact
### Merge the peak matrix with the temporal clusters
############################################################
############################################################
#### Import contact data
Peakmatrix_all <- fread(Peakmatrix_raw_score_file)[,c("baitID", "oeID", "naive", "day1", "day3", "day5", "day7", "day10", "day14", "baitChr", "baitStart", "baitEnd", "baitName", "oeChr", "oeStart", "oeEnd", "oeName")]
dim(Peakmatrix_all) # 1322199      17

dtw_clusters <- fread(dtw_clusters_file)[, .(baitID, oeID, merged_dtw_trans_cl6_2with6)]
dim(dtw_clusters) # 9998    3

peakMatrix <- merge(Peakmatrix_all, dtw_clusters, by = c("baitID", "oeID"))
dim(peakMatrix) # 9998   18

baitmap <- fread(baitmap_file)
setnames(baitmap, c("V1", "V2", "V3", "V4"), c("chrom", "start", "end", "baitID"))
baitmap[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]

############################################################
############################################################
### General common parameters and files
############################################################
############################################################

########################################
## Timepoints
########################################
timepoints_labels  <- c("Naive","Day1","Day3","Day5","Day7","Day10","Day14")
timepoints_columns <- c("naive","day1","day3","day5","day7","day10","day14")

########################################
## bigWig files
########################################
k27me3_files <- c(
  naive="merged_Naive_273.bam.bw",
  day1 ="merged_Day1_273.bam.bw",
  day3 ="merged_Day3_273.bam.bw",
  day5 ="merged_Day5_273.bam.bw",
  day7 ="merged_Day7_273.bam.bw",
  day10="merged_Day10_273.bam.bw",
  day14="merged_Day14_273.bam.bw"
)

k4me1_files <- c(
  naive="merged_Naive_41.bam.bw",
  day1 ="merged_Day1_41.bam.bw",
  day3 ="merged_Day3_41.bam.bw",
  day5 ="merged_Day5_41.bam.bw",
  day7 ="merged_Day7_41.bam.bw",
  day10="merged_Day10_41.bam.bw",
  day14="merged_Day14_41.bam.bw"
)

############################################################
############################################################
### Plot an example for a CONSTANT contact
############################################################
############################################################

########################################
########## Selection of the contact of interest
########################################
peakMatrix_clusters_NtoP_filtered_Constant <- peakMatrix[
  baitStart == 46492257 & oeStart == 46545138
]

bait_of_interest_constant <- baitmap[
  chrom == "chr1" & start == 46492257 & end == 46492789
]

########################################
## Build Hi-C list
########################################
hic_list <- lapply(timepoints_columns, function(tp_col) {
  df <- peakMatrix_clusters_NtoP_filtered_Constant[, .(
    baitChr = paste0("chr", baitChr),
    baitStart, baitEnd,
    oeChr = paste0("chr", oeChr),
    oeStart, oeEnd,
    score = get(tp_col)
  )]
  return(df)
})
names(hic_list) <- timepoints_labels

########################################
## Parameters
########################################
params <- pgParams(
  chrom = "chr1", chromstart = 46400000, chromend = 46600000,
  assembly = "hg38", x = unit(0.5,"inches"), width = unit(2,"inches")
)

params_bigwig <- pgParams(
  chrom = "1", chromstart = 46400000, chromend = 46600000,
  assembly = "hg38", x = unit(0.5,"inches"), width = unit(2,"inches")
)

########################################
## Color: arcsinh + global scale
########################################
score_palette <- c("white", "darkgreen")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$score_asinh <- asinh(df$score)
  hic_list[[i]] <- df
}

# global max
all_vals <- unlist(lapply(hic_list, function(x) x$score_asinh))
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE)

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  hic_list[[i]] <- df
}

########################################
## PDF + PAGE
########################################
pdf(file.path(root_dir, "Figure3_CnT_PECHiC_NtoDay14_constant_arcsinhColor_bigger0.8_v3.pdf"), width = 7, height = 32)

pageCreate(width = 4, height = 15, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.8,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_constant, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

############################################################
## Layout constants
############################################################
y_range_k27me3 <- c(0,2000)
y_range_k4me1  <- c(0,300)

timepoint_spacing <- 3
k27me3_offset <-0.5
k4me1_offset  <- 1.5

############################################################
## LOOP: ALL TIMEPOINTS
############################################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.2 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_file_dir, k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k27me3, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_file_dir, k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k4me1, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkgreen"))(256 )
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k4me1_offset + 1
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.4,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    plotText(
      label=paste(tp_label,"PECHiC"),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=23
)

## label
plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkgreen"))(256) , alpha.f = 1) # currently no transparency as in the default,
# but leave it like this so we can vary alpha
# here and in plotPairsArcs()

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(6 * sf, "inches"),
  y = unit(6 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((6-5/2) * sf, (6+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((6.4 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score – ", round(global_max_asinh,2), ")"),
  x = unit(5.8 * sf, "inches"),
  y = unit(8 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()

############################################################
############################################################
### Plot an example for a  GAINED early contact
############################################################
############################################################

########################################
## Selection of the contact of interest
########################################
peakMatrix_clusters_NtoP_filtered_Gained_early <- peakMatrix[baitStart == 27109583 & oeStart == 27153301]  
head(peakMatrix_clusters_NtoP_filtered_Gained_early)  

bait_of_interest_gained_early <- baitmap[  
  chrom == "chr7" &  
    start == 27109583 &  
    end == 27110604  
]  

########################################
## Build Hi-C list
########################################
hic_list <- lapply(timepoints_columns, function(tp_col) {
  df <- peakMatrix_clusters_NtoP_filtered_Gained_early[, .(
    baitChr = paste0("chr", baitChr),
    baitStart,
    baitEnd,
    oeChr = paste0("chr", oeChr),
    oeStart,
    oeEnd,
    score = get(tp_col)
  )]
  return(df)
})
names(hic_list) <- timepoints_labels

########################################
## Parameters
########################################
params <- pgParams(
  chrom = "chr7",
  chromstart = 27050000, chromend = 27200000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)


params_bigwig <- pgParams(
  chrom = "7",
  chromstart = 27050000, chromend = 27200000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

########################################
## Color: arcsinh + global scale
########################################
score_palette <- c("white", "darkred")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$score_asinh <- asinh(df$score)
  hic_list[[i]] <- df
}

# global max
all_vals <- unlist(lapply(hic_list, function(x) x$score_asinh))
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE)

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  hic_list[[i]] <- df
}

########################################
## PDF + PAGE
########################################
pdf(file.path(root_dir, "Figure3_CnT_PECHiC_NtoDay14_gainedEarly_arcsinhColor_bigger0.8.pdf"), width = 7, height = 32)

pageCreate(width = 4, height = 15, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.8,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_gained_early, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

############################################################
## Layout constants
############################################################
y_range_k27me3 <- c(0,4000)
y_range_k4me1  <- c(0,300)

timepoint_spacing <- 3
k27me3_offset <- 0.5
k4me1_offset  <- 1.5

############################################################
## LOOP: ALL TIMEPOINTS
############################################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.2 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_file_dir, k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k27me3, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_file_dir, k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k4me1, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkred"))(256 )
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k4me1_offset + 1
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.4,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    
    plotText(
      label=paste(tp_label,"PECHiC"),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=23
)

## label
plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkred"))(256) , alpha.f = 1) # currently no transparency as in the default,
# but leave it like this so we can vary alpha
# here and in plotPairsArcs()

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(6 * sf, "inches"),
  y = unit(6 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((6-5/2) * sf, (6+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((6.4 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score – ", round(global_max_asinh,2), ")"),
  x = unit(5.8 * sf, "inches"),
  y = unit(8 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()


############################################################
############################################################
### Plot an example for a  GAINED late contact
############################################################
############################################################

########################################
## Selection of the contact of interest
########################################
peakMatrix_clusters_NtoP_filtered_Gained_late <- peakMatrix[baitStart == 37096272 & oeStart == 36878651]  
head(peakMatrix_clusters_NtoP_filtered_Gained_late)  

bait_of_interest_gained_late <- baitmap[  
  chrom == "chr15" &  
    start == 37096272 &  
    end == 37096646  
]  

########################################
## Build Hi-C list
########################################
# Prepare Hi-C data list with height for arcs
hic_list <- lapply(timepoints_columns, function(tp_col) {
  df <- peakMatrix_clusters_NtoP_filtered_Gained_late[, .(
    baitChr = paste0("chr", baitChr),
    baitStart,
    baitEnd,
    oeChr = paste0("chr", oeChr),
    oeStart,
    oeEnd,
    score = get(tp_col)
  )]
  return(df)
})
names(hic_list) <- timepoints_labels


########################################
## Parameters
########################################
params <- pgParams(
  chrom = "chr15",
  chromstart = 36800000, chromend = 37140000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)


params_bigwig <- pgParams(
  chrom = "15",
  chromstart = 36800000, chromend = 37140000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

########################################
## Color: arcsinh + global scale
########################################
score_palette <- c("white", "darkorange")
color_fun <- colorRamp(score_palette, space = "Lab")

# arcsinh transform
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$score_asinh <- asinh(df$score)
  hic_list[[i]] <- df
}

# global max
all_vals <- unlist(lapply(hic_list, function(x) x$score_asinh))
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE)

# mapping function
getScoreColors <- function(scores, max_val) {
  scores[is.na(scores)] <- 0
  scaled <- pmin(scores / max_val, 1)
  rgb_vals <- color_fun(scaled)
  rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue=255)
}

# apply colors
for(i in seq_along(hic_list)) {
  df <- hic_list[[i]]
  df$color_asinh <- getScoreColors(df$score_asinh, global_max_asinh)
  hic_list[[i]] <- df
}

########################################
## PDF + PAGE
########################################
pdf(file.path(root_dir, "Figure3_CnT_PECHiC_NtoDay14_gainedLate_arcsinhColor_bigger0.8.pdf"), width = 7, height = 32)

pageCreate(width = 4, height = 15, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.8,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_gained_late, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

############################################################
## Layout constants
############################################################
y_range_k27me3 <- c(0,1200)
y_range_k4me1  <- c(0,200)

timepoint_spacing <- 3
k27me3_offset <- 0.5
k4me1_offset  <- 1.5

############################################################
## LOOP: ALL TIMEPOINTS
############################################################
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.2 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_file_dir, k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k27me3, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_file_dir, k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.8,
    just=c("left","top"), range=y_range_k4me1, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkorange"))(256 )
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k4me1_offset + 1
    
    plotPairsArches(
      data=df_arc, params=params,
      fill=df_arc$color_asinh,
      linecolor=df_arc$color_asinh,
      y= y_arc, 
      width=3, height=0.4,
      just=c("left","top"), flip=TRUE, 
      alpha=1 # remove the default transparency mask - very important for colour-coded arcs!
    )
    
    
    plotText(
      label=paste(tp_label,"PECHiC"),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=23
)

## label
plotText(label="Genes", x=0.2, y=0.45, fontsize=8)

## color scale legend
## size reduction factor
sf <- 0.7

## color scale legend (smaller)
legend_colors <- adjustcolor(colorRampPalette(c("white", "darkorange"))(256) , alpha.f = 1) # currently no transparency as in the default,
# but leave it like this so we can vary alpha
# here and in plotPairsArcs()

grid.raster(
  matrix(rev(legend_colors), ncol = 1),
  x = unit(6 * sf, "inches"),
  y = unit(6 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((6-5/2) * sf, (6+5/2) * sf, length.out = length(breaks_labels))

for(i in seq_along(breaks_labels)){
  grid.text(
    as.character(breaks_labels[i]),
    x = unit((6.4 * sf), "inches"),
    y = unit(ylab_positions[i], "inches"),
    gp = gpar(cex = 0.6)   ## smaller text
  )
}

grid.text(
  paste0("arcsinh Chicago score – ", round(global_max_asinh,2), ")"),
  x = unit(5.8 * sf, "inches"),
  y = unit(8 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()
