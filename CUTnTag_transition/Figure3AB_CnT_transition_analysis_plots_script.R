################################################################################
### Dynamics of the poising-associated histone modifications H3K27me3 (3A) and H3K4me1 (3B) during the naive-to-primed transition profiled by CUT&Tag - Figure 3
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 3 panels showing:
###  Mean RPKM-normalized signals across PEs per time point, normalised to the levels detected in the naive state. 
###  The same values averaged for PEs engaged in each temporal class of contacts. 
################################################################################

########################################
### Load required libraries
########################################
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(Hmisc)

########################################
### User configuration
########################################
## Root directory containing all input files and custom functions
root_dir <- "~/Documents/Bioinformatics"

## Input files
pileup_file  <- file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/MN_CnT_273_41_PileUp_NormComb_5kbRmap_v2.txt")
baitmap_file <- file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")
cluster_file <- file.path("~/Documents/Bioinformatics/Final_scripts_paper/pepm_transition_dtw_final_02-02-2026.txt")

########################################
### Load and prepare CUT&Tag data
########################################
## Upload normalised CUT&Tag pileup matrix (5 kb bins) for H3K27me3 and H3K4me1
PileUp_norm_41_273 <- fread(pileup_file)

## Upload PE-Capture Hi-C bait annotation
baitmap <- fread(baitmap_file)

## Merge CUT&Tag signal with bait information
## V4 corresponds to bait ID in PE-Capture Hi-C baitmap
baits.plp <- PileUp_norm_41_273 %>%
  inner_join(baitmap, by = c("ID" = "V4")) %>%
  mutate(isBait = TRUE)

## Load temporal (DTW) and chromatin bait clustering
dtwClustersBaitClust <- fread(cluster_file)
dtwClustersBaitClust <- dtwClustersBaitClust[,c("baitID", "oeID", "merged_dtw_trans_cl6_2with6")]

## Merge clustering information with CUT&Tag data
baitsplp_clust <- merge(
  baits.plp,
  dtwClustersBaitClust,
  by.x = "ID",
  by.y = "baitID"
)

## Keep one row per bait ID (when a bait is involved in several contacts just keep one)
baitsplp_clustUn <- unique(setDT(baitsplp_clust), by = "ID")


########################################
### RPKM calculation
########################################
div <- function(x, y = baitsplp_clustUn$fragLen) {x/y}

## Compute RPKM values for CUT&Tag signal columns
RPKM <- div(apply(baitsplp_clustUn[,6:25], 2, function(x)x))
RPKM <- setDT(as.data.frame(RPKM))

## Combine RPKM values with metadata
Ncol <- ncol(baitsplp_clustUn)
RPKM_cl <- cbind(
  baitsplp_clustUn[, c(1:5, 32), with = FALSE],
  RPKM
)


########################################
### Shared plotting parameters
########################################
custom_colors <- c(
  "Constant"  = "green",
  "Gained early"  = "red",
  "Lost"      = "blue",
  "Gained late"   = "orange",
  "Transient" = "purple"
)

################################################################################
### H3K27me3 CUT&Tag analysis
################################################################################
## Reshape RPKM values to long format (keep naive to day14)
RPKM.long.k27 <- as.data.table(
  RPKM_cl %>%
    pivot_longer(
      cols = 7:13,
      names_to = "timepoint",
      values_to = "RPKMvalue"
    ) %>%
    mutate(time = gsub("day(\\d+)_.+", "\\1", timepoint)) %>%
    dplyr::select(
      ID,
      merged_dtw_trans_cl6_2with6,
      timepoint,
      RPKMvalue,
      time
    )
)

## Define naive timepoint
RPKM.long.k27[time == "N_273", time := 0]
RPKM.long.k27[, time := as.numeric(time)]
RPKM.long.k27[, timeFac := factor(time, levels = c(0,1,3,5,7,10,14))]

## Remove baits without temporal/contact cluster assignment
RPKM.long.k27 <- RPKM.long.k27[merged_dtw_trans_cl6_2with6 != ""]


########################################
### H3K27me3 – mean relative to naive
########################################
naive_rpkm_k27 <- RPKM.long.k27[time == 0, mean(RPKMvalue)]

RPKM.long.k27[, RPKM_relative := RPKMvalue / naive_rpkm_k27]

pdf(file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/Figure3A_left_panel_K27me3.pdf"))
ggplot(RPKM.long.k27, aes(x = time, y = RPKM_relative)) +
  stat_summary(fun = mean, geom = "line", size = 1.5, color = "darkgrey") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  theme_minimal() +
  labs(
    y = "Mean RPKM relative to naive"
  ) +
  scale_x_continuous(breaks = c(0,1,3,5,7,10,14)) +
  scale_y_continuous(breaks = seq(0, 3.5, 0.5)) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme(legend.position = "none")
dev.off()


########################################
### H3K27me3 – per cluster relative to naive
########################################
naive_means_k27 <- RPKM.long.k27[
  time == 0,
  .(naive_mean = mean(RPKMvalue, na.rm = TRUE)),
  by = merged_dtw_trans_cl6_2with6
]

RPKM.long.k27 <- merge(
  RPKM.long.k27,
  naive_means_k27,
  by = "merged_dtw_trans_cl6_2with6"
)

RPKM.long.k27[, RPKM_relative := RPKMvalue / naive_mean]

pdf(file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/Figure3A_right_panel_K27me3.pdf"))
ggplot(
  RPKM.long.k27,
  aes(
    x = time,
    y = RPKM_relative,
    color = merged_dtw_trans_cl6_2with6,
    group = merged_dtw_trans_cl6_2with6
  )
) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal() +
  labs(
    y = "Mean RPKM H3K27me3 relative to naive",
    color = "Contact cluster"
  ) +
  scale_x_continuous(breaks = c(0,1,3,5,7,10,14)) +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  coord_cartesian(ylim = c(0, 4)) +
  scale_color_manual(values = custom_colors)
dev.off()


################################################################################
### H3K4me1 analysis
################################################################################
RPKM.long.k4 <- as.data.table(
  RPKM_cl %>%
    pivot_longer(
      cols = 17:23,
      names_to = "timepoint",
      values_to = "RPKMvalue"
    ) %>%
    mutate(time = gsub("day(\\d+)_.+", "\\1", timepoint)) %>%
    dplyr::select(
      ID,
      merged_dtw_trans_cl6_2with6,
      timepoint,
      RPKMvalue,
      time
    )
)

RPKM.long.k4[time == "N_41", time := 0]
RPKM.long.k4[, time := as.numeric(time)]
RPKM.long.k4[, timeFac := factor(time, levels = c(0,1,3,5,7,10,14))]

## Remove baits without temporal/contact cluster assignment
RPKM.long.k4 <- RPKM.long.k4[merged_dtw_trans_cl6_2with6 != ""]


########################################
### H3K4me1 – mean relative to naive
########################################
naive_rpkm_k4 <- RPKM.long.k4[time == 0, mean(RPKMvalue)]
RPKM.long.k4[, RPKM_relative := RPKMvalue / naive_rpkm_k4]

pdf(file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/Figure3B_left_panel_K4me1.pdf"))
ggplot(RPKM.long.k4, aes(x = time, y = RPKM_relative)) +
  stat_summary(fun = mean, geom = "line", size = 1.5, color = "darkgrey") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  theme_minimal() +
  labs(
    y = "Mean RPKM relative to naive"
  ) +
  scale_x_continuous(breaks = c(0,1,3,5,7,10,14)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(ylim = c(0, 2)) +
  theme(legend.position = "none")
dev.off()

########################################
### H3K4me1 – per cluster relative to naive
########################################
naive_means_k4 <- RPKM.long.k4[
  time == 0,
  .(naive_mean = mean(RPKMvalue, na.rm = TRUE)),
  by = merged_dtw_trans_cl6_2with6
]

RPKM.long.k4 <- merge(
  RPKM.long.k4,
  naive_means_k4,
  by = "merged_dtw_trans_cl6_2with6"
)

RPKM.long.k4[, RPKM_relative := RPKMvalue / naive_mean]

pdf(file.path(root_dir, "Bioinfo_REDO_figure_Monica_CnT/Figure3B_right_panel_K4me1.pdf"))
ggplot(
  RPKM.long.k4,
  aes(
    x = time,
    y = RPKM_relative,
    color = merged_dtw_trans_cl6_2with6,
    group = merged_dtw_trans_cl6_2with6
  )
) +
  stat_summary(fun = mean, geom = "line", size = 1.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal() +
  labs(
    y = "Mean RPKM H3K4me1 relative to naive",
    color = "Contact cluster"
  ) +
  scale_x_continuous(breaks = c(0,1,3,5,7,10,14)) +
  scale_y_continuous(breaks = seq(0, 2, 0.5)) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_color_manual(values = custom_colors)
dev.off()
