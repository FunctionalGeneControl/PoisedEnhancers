######################################################################################################
### Figure 6 and Supplemental Figure S4 - PE contacts persist after ectopic CRISPR activation targeting and can mediate long-range gene induction
### Author: Dr. Marina Nocente
###
### This script reproduces some panels of Figure 6 and Supplemental Figure S4
###  Figure 6E and Supplemental Figure S4E - Visualisation of the H3K27me3 and H3K4me1 CUT&Tag signals and the ‘gained-early’ PECHi-C contact
###                                          between PE 280801 and the DLX1 promoter during the naive-to-primed transition
###
###  Figure 6F - Normalised gene expression levels of DLX1, DLX2 and METAP1D following inducible CRISPR activation (CRISPRa)
###              of PE 280801 versus controls using scrambled gRNAs
###
######################################################################################################

#######################################
## Libraries
#######################################
library(data.table)
library(dplyr)

library(ggplot2)
library(plotgardener)
library(grid)

#######################################
## File paths
#######################################
RTqPCR_dir <- file.path("~/Documents/Bioinformatics/Scripts_R_graphs")

working_dir <- file.path("~/Documents/Bioinformatics/Find_examples")

Peakmatrix_file <- file.path("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
TssOE_file <- file.path("~/Documents/Bioinformatics/GO_OtherEnds/OtherEnd_TSS_coordinates.txt")
baitmap_file <- file.path("~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")

CUTnTag_dir <- file.path("~/Documents/Bioinformatics/CnT_transition/bigwig_merged_just_sorted_indexedBamFiles/")


##############################################################################
##############################################################################
## Figure 6E and Supplemental Figure S4E
##############################################################################
##############################################################################

#######################################
## Upload data
#######################################
## Upload the peakmatrix with the Naive-to-primed Chicago score (PECHiC data)
Peakmatrix_all <- fread(Peakmatrix_file)
dim(Peakmatrix_all) # 1322199      24

## How many other ends are poised enhancers?
sum(grepl("poised enhancer", Peakmatrix_all$oeName, ignore.case = TRUE)) # 29385

Peakmatrix_all_oePE <- Peakmatrix_all[grepl("poised enhancer", oeName, ignore.case = TRUE),]
dim(Peakmatrix_all_oePE) # 29385    24

## How many other ends are "other things" = TSS?
sum(!grepl("poised enhancer", Peakmatrix_all$oeName, ignore.case = TRUE)) # 1292814

Peakmatrix_all_oeTSS <- Peakmatrix_all[!grepl("poised enhancer", oeName, ignore.case = TRUE),]
dim(Peakmatrix_all_oeTSS) # 1292814      24

## Upload coordinates of Other Ends which are TSS (with gene_id and gene_names)
TssOE <- fread(TssOE_file)
TssOE <- TssOE[,c("oeID", "baitID", "TSSchr", "TSSstart", "TSSend", "TSSstrand", "gene_name", "gene_id", "chr", "start", "end")] 
dim(TssOE) # 7538   11

Peakmatrix_all_oeTSS <- merge(Peakmatrix_all_oeTSS, TssOE, by=c("baitID", "oeID"))
dim(Peakmatrix_all_oeTSS) # 7486   35

## Upload the baitmap data
baitmap <- fread(baitmap_file)
setnames(baitmap, c("V1", "V2", "V3", "V4"), c("chrom", "start", "end", "baitID"))
baitmap[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]

#######################################
## Selection of the contact of interest:
## Contact between PE n°280801 and the DLX1 gene
#######################################
## Keep just this interaction: DLX1 gene
# baitID = 280801 - 2 172095943 172096312
# oeID = 280785 - 2   172082151 172084569

Peakmatrix_all_oeTSS_DLX1 <- Peakmatrix_all_oeTSS[baitStart == 172095943 & oeStart == 172082151]
head(Peakmatrix_all_oeTSS_DLX1)

bait_of_interest_upEarly_DLX1 <- baitmap[
  chrom == "chr2" & 
    start == 172095943 & 
    end == 172096312
]

## Interaction and score for each timepoint 
hic_df_naive_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = naive  
)]

hic_df_day1_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day1  
)]

hic_df_day3_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day3  
)]

hic_df_day5_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day5  
)]

hic_df_day7_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day7 
)]

hic_df_day10_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day10 
)]

hic_df_day14_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14 
)]

hic_df_primed_UpEarly_DLX1 <- Peakmatrix_all_oeTSS_DLX1[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = primed 
)]

#######################################
## Define variables
#######################################
## Timepoints
timepoints_labels  <- c("Naive","Day1","Day3","Day5","Day7","Day10","Day14", "Primed")
timepoints_columns <- c("naive","day1","day3","day5","day7","day10","day14", "primed")

## CUT&Tag bigWig files
k27me3_files <- c(
  naive="merged_Naive_273.bam.bw",
  day1 ="merged_Day1_273.bam.bw",
  day3 ="merged_Day3_273.bam.bw",
  day5 ="merged_Day5_273.bam.bw",
  day7 ="merged_Day7_273.bam.bw",
  day10="merged_Day10_273.bam.bw",
  day14="merged_Day14_273.bam.bw",
  primed="merged_Primed_273.bam.bw"
)

k4me1_files <- c(
  naive="merged_Naive_41.bam.bw",
  day1 ="merged_Day1_41.bam.bw",
  day3 ="merged_Day3_41.bam.bw",
  day5 ="merged_Day5_41.bam.bw",
  day7 ="merged_Day7_41.bam.bw",
  day10="merged_Day10_41.bam.bw",
  day14="merged_Day14_41.bam.bw",
  primed="merged_Day14_41.bam.bw"
)

## Build Hi-C list
hic_list <- list(
  hic_df_naive_UpEarly_DLX1,
  hic_df_day1_UpEarly_DLX1,
  hic_df_day3_UpEarly_DLX1,
  hic_df_day5_UpEarly_DLX1,
  hic_df_day7_UpEarly_DLX1,
  hic_df_day10_UpEarly_DLX1,
  hic_df_day14_UpEarly_DLX1,
  hic_df_primed_UpEarly_DLX1
)

names(hic_list) <- timepoints_labels

## Parameters for Plotgardener
# Genomic region
params <- pgParams(
  chrom = "chr2", 
  chromstart = 172070000, chromend = 172111000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "2", 
  chromstart = 172070000, chromend = 172111000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Color: arcsinh + global scale
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
#     Naive      Day1      Day3      Day5      Day7     Day10     Day14    Primed
# 0.2684735 1.5455407 0.5337224 0.7849851 1.4547205 1.9008881 1.6946927 2.4164895 
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 2.426489 when we include primed

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

#######################################
## Plot
#######################################
pdf(file.path(working_dir, "Figure6E_SuppFigS4E_CnT_PECHiC_NtoP_gainedEarly_DLX1_arcsinhColor.pdf"), width = 7, height = 20)

pageCreate(width = 7, height = 20, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_upEarly_DLX1, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
# Cut&Tag y-ranges
y_range_k27 <- c(0, 2000)
y_range_k4  <- c(0, 300)

timepoint_spacing <- 2
k27me3_offset <- 0.8
k4me1_offset  <- 1.6

## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.2 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.2, y=base_y+0.5)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkred"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k4me1_offset + 0.7
    
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
      label=paste(tp_label,"PECHiC"),
      fontsize=9, fontcolor="darkblue",
      x=0.8, y=y_arc - 0.1
    )
  }
}

## Genome axis
annoGenomeLabel(
  plot = genesPlot, params=params,
  scale="Mb", y=18
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


##############################################################################
##############################################################################
## Figure 6F
##############################################################################
##############################################################################

#######################################
## RT-qPCR normalised data
#######################################
# METAP1D_primer Scrmbl: 1, 1, 1, 1
# METAP1D_primer Dlx1: 1.39, 0.44, 1.12, 1.12
# Dlx1_primer Scrmbl: 1, 1, 1, 1
# Dlx1_primer Dlx1: 18.16, 8.90, 16.12, 5.41
# DLX2_primer Scrmbl: 1, 1, 1, 1
# DLX2_primer Dlx1: 0.39, 4.92, 1.06, 1.33


## Create a data frame containing the RT-qPCR normalised data
df_RTqPCR_Dlx1 <- data.frame(
  Primer = c("METAP1D Scramble", "METAP1D Scramble", "METAP1D Scramble", "METAP1D Scramble", "METAP1D Dlx1", "METAP1D Dlx1", 
             "METAP1D Dlx1", "METAP1D Dlx1", "DLX1 Scramble", "DLX1 Scramble", "DLX1 Scramble", "DLX1 Scramble", "DLX1 Dlx1", "DLX1 Dlx1", "DLX1 Dlx1", "DLX1 Dlx1", 
             "DLX2 Scramble", "DLX2 Scramble", "DLX2 Scramble", "DLX2 Scramble",
             "DLX2 Dlx1", "DLX2 Dlx1", "DLX2 Dlx1", "DLX2 Dlx1"),
  Value = c(1, 1, 1, 1, 1.39, 0.44, 1.12, 1.12, 1, 1, 1, 1, 18.16, 8.90, 16.12, 5.41, 1, 1, 1, 1, 0.39, 4.92, 1.06, 1.33)
)

df_RTqPCR_Dlx1$Primer <- factor(df_RTqPCR_Dlx1$Primer, levels = unique(df_RTqPCR_Dlx1$Primer))

## Means and sd calculation
means <- tapply(df_RTqPCR_Dlx1$Value, df_RTqPCR_Dlx1$Primer, mean, na.rm = TRUE)
sds <- tapply(df_RTqPCR_Dlx1$Value, df_RTqPCR_Dlx1$Primer, sd, na.rm = TRUE)

## Plot
pdf(file.path(RTqPCR_dir, "Figure6F_RTqPCR.pdf"), width = 4, height = 6)
ggplot(df_RTqPCR_Dlx1, aes(
  # Primer extraction (Dlx1, METAP1D, etc.)
  x = factor(gsub(" .*", "", Primer), levels = c("METAP1D", "DLX1", "DLX2")), 
  y = Value, 
  # Sample extraction (Scramble ou Dlx1)
  fill = factor(gsub(".* ", "", Primer), levels = c("Scramble", "Dlx1"))
)) +
  geom_bar(
    stat = "summary", 
    fun = "mean", 
    width = 0.7, 
    position = position_dodge(width = 0.8) 
  ) +
  geom_errorbar(
    stat = "summary", 
    fun.data = function(x) {c(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))}, 
    width = 0.2, 
    position = position_dodge(width = 0.8) 
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
    color = "black", 
    size = 3, 
    alpha = 0.7
  ) +
  labs(
    title = "",
    x = "", 
    y = "Normalised gene expression", 
    fill = "Samples"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.4, size = 14)
  ) +
  scale_fill_manual(values = c("Scramble" = "dodgerblue", "Dlx1" = "purple")) + 
  scale_y_continuous(breaks = seq(0, 20, by = 2), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0, 20))

dev.off()


#######################################
## Statistics
#######################################
# Separate the data in two groups : "DLX1 Scramble" and "DLX1 Dlx1"
dlx1_scramble <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "DLX1 Scramble"]
dlx1_dlx1 <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "DLX1 Dlx1"]

# Separate the data in two groups : "DLX2 Scramble" and "DLX2 Dlx1"
dlx2_scramble <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "DLX2 Scramble"]
dlx2_dlx1 <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "DLX2 Dlx1"]

# Separate the data in two groups : "METAP1D Scramble" and "METAP1D Dlx1"
metap1d_scramble <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "METAP1D Scramble"]
metap1d_dlx1 <- df_RTqPCR_Dlx1$Value[df_RTqPCR_Dlx1$Primer == "METAP1D Dlx1"]

### Wilcoxon-Mann-Whitney test
wilcox.test(dlx1_scramble, dlx1_dlx1, exact = FALSE) # p-value = 0.02107
wilcox.test(dlx2_scramble, dlx2_dlx1, exact = FALSE) # p-value = 0.2817
wilcox.test(metap1d_scramble, metap1d_dlx1, exact = FALSE) # p-value = 0.2784


