################################################################################
### Integration of PECHi-C and CUT&Tag data after PRC2 perturbation at Day5 of the transition
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 5 panels D and E, and Figure S3F showing:
###  Figure 5D: Changes in H3K27me3 signal at PEs after treatment with PRC2 PROTAC (red) or PRC2 catalytic inhibitor (blue) relative to DMSO, stratified by bins of log10-H3K27me3 signal in the DMSO condition. 
###  Figure 5E: Changes in connectivity at PEs after treatment with PRC2 PROTAC (red) or PRC2 catalytic inhibitor (blue) relative to DMSO, stratified by bins of log10-H3K27me3 signal in the DMSO condition. 
###  Figure S3F: Changes in the overall connectivity of PEs after treatment with PRC2 PROTAC or PRC2 catalytic inhibitor
################################################################################

###########################################
### Libraries
###########################################
library(data.table)
library(dplyr)
library(tidyr)

library(ggplot2)
library(rstatix)

###########################################
### File paths
###########################################
base_dir <- "~/Documents/Bioinformatics"
working_dir <- "~/Documents/Bioinformatics/CnT_PRC2_perturb_Day5"

CUTnTag_deseq2_results <- "~/Documents/Project_writing/Paper_PEs/Data_release/CUT&Tag_PRC2/final/"
PECHiC_deseq2_results <- "~/Documents/Project_writing/Paper_PEs/Data_release/PECHiC_PRC2/"

temporal_clusters_file <- file.path(base_dir, "Final_scripts_paper/pepm_transition_dtw_final_02-02-2026.txt")
baitmap_file <- file.path(base_dir, "CnT_PRC2_perturb_Day5/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt")

###########################################
### Upload and prepare the PECHi-C data after PRC2 perturbation and DESeq2 analysis
###########################################
## Upload the DESeq2 results on read counts per bait (PECHi-C results)
res_PROTAC_vs_DMSO_merged <- fread(file.path(PECHiC_deseq2_results, "res_deseq2PECHiC_PROTAC_DMSO_merged_all.txt"))
res_UNC1999_vs_DMSO_merged <- fread(file.path(PECHiC_deseq2_results, "res_deseq2PECHiC_UNC_DMSO_merged_all.txt"))

# Rename the columns
colnames(res_PROTAC_vs_DMSO_merged)[5] <- "stat_res_PROTAC_vs_DMSO_merged"
colnames(res_UNC1999_vs_DMSO_merged)[5] <- "stat_res_UNC1999_vs_DMSO_merged"

## Create a final object containing all what we need
PECHiC_deseq_PRCperturb_stats <- res_UNC1999_vs_DMSO_merged[,c("baitID", "stat_res_UNC1999_vs_DMSO_merged")]
PECHiC_deseq_PRCperturb_stats <- merge(PECHiC_deseq_PRCperturb_stats, 
                                       res_PROTAC_vs_DMSO_merged[, c("baitID", "stat_res_PROTAC_vs_DMSO_merged")],
                                       by = "baitID")

head(PECHiC_deseq_PRCperturb_stats)
dim(PECHiC_deseq_PRCperturb_stats) # 31754 baits with PROTAC and UNC stats


## Upload the temporal contact clustering data
TempClusters <- fread(temporal_clusters_file)
TempClusters <- TempClusters[,c("baitID", "oeID", "merged_dtw_trans_cl6_2with6")]


## Merge the PECHiC deseq2 results with the temporal clusters
PECHiC_deseq_PRCperturb_stats_TempClusters <- merge(PECHiC_deseq_PRCperturb_stats, 
                                                    TempClusters, 
                                                    by = "baitID")

dim(PECHiC_deseq_PRCperturb_stats_TempClusters) # 9998    5

length(PECHiC_deseq_PRCperturb_stats_TempClusters$baitID) # 9998 rows=baitID, several times the same baitID
length(unique(PECHiC_deseq_PRCperturb_stats_TempClusters$baitID)) # only 6445 unique baitID

fwrite(PECHiC_deseq_PRCperturb_stats_TempClusters,file.path("~/Documents/Project_writing/Paper_PEs/Data_release/PECHiC_deseq2_stats_PRC2-perturb_Tempclusters.txt"), sep = "\t", quote = FALSE)


## Add the coordinates of the bait
baitmap <- fread(baitmap_file)
colnames(baitmap) <- c("chr", "start", "end", "baitID")
dim(baitmap) # 32142     4

PECHiC_deseq_PRCperturb_stats_TempClusters <- merge(PECHiC_deseq_PRCperturb_stats_TempClusters, baitmap, by="baitID") 


###########################################
### Upload and prepare the CUT&Tag data after PRC2 perturbation and DESeq2 analysis
###########################################
# PROTAC vs DMSO (Sequencing Run2) normal DESeq2
res_PROTAC_DMSO_CnT_all_dmsoBin_normal = fread(file.path(CUTnTag_deseq2_results, "CUTnTag_deseq2_res_PROTAC_DMSO_dmsoMean_RunSeq2.txt"))
res_PROTAC_DMSO_CnT_all_dmsoBin_normal[, c("chr", "start", "end") := tstrsplit(region, "_")]
res_PROTAC_DMSO_CnT_all_dmsoBin_normal[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
res_PROTAC_DMSO_CnT_all_dmsoBin_normal <- res_PROTAC_DMSO_CnT_all_dmsoBin_normal[,-1] # 32142    11

res_UNC1999_DMSO_CnT_all_dmsoBin_normal = fread(file.path(CUTnTag_deseq2_results, "CUTnTag_deseq2_res_UNC1999_DMSO_dmsoMean_RunSeq2.txt"))
res_UNC1999_DMSO_CnT_all_dmsoBin_normal[, c("chr", "start", "end") := tstrsplit(region, "_")]
res_UNC1999_DMSO_CnT_all_dmsoBin_normal[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
res_UNC1999_DMSO_CnT_all_dmsoBin_normal <- res_UNC1999_DMSO_CnT_all_dmsoBin_normal[,-1] # 32142    11


###########################################
### Integrate PECHi-C and CUT&Tag data
###########################################
## For UNC and PROTAC all baits with a temporal class
PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal = merge(PECHiC_deseq_PRCperturb_stats_TempClusters, res_PROTAC_DMSO_CnT_all_dmsoBin_normal, by=c("chr", "start", "end")) 
# 9998

PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal = merge(PECHiC_deseq_PRCperturb_stats_TempClusters, res_UNC1999_DMSO_CnT_all_dmsoBin_normal, by=c("chr", "start", "end")) 
# 9998


###########################################
### Visualisation: changes in connectivity overall - 6,445 baits
###########################################
### Select the unique baitID (=PEs) to avoid multi-coumting
PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal_unique <- unique(PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal, by="baitID") 
dim(PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal_unique) # 6445   16

PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal_unique <- unique(PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal, by="baitID") 
dim(PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal_unique) # 6445   16


## stat PROTAC vs UNC1999 (no filter)
df_long <- 
  PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal_unique %>%
  select(chr, start, end,
         stat_res_PROTAC_vs_DMSO_merged,
         stat_res_UNC1999_vs_DMSO_merged) %>%
  pivot_longer(
    cols = c(stat_res_PROTAC_vs_DMSO_merged,
             stat_res_UNC1999_vs_DMSO_merged),
    names_to = "Treatment",
    values_to = "stat"
  )

df_long$Treatment <- recode(df_long$Treatment,
                            "stat_res_PROTAC_vs_DMSO_merged" = "PROTAC",
                            "stat_res_UNC1999_vs_DMSO_merged" = "UNC1999"
)

pdf(file.path(working_dir, "Total_connectivity_UNC1999_PROTAC_boxplot.pdf"), width = 6, height = 4)
ggplot(df_long,
       aes(x = Treatment, y = stat, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  #geom_jitter(width = 0.15, size = 0.2, alpha = 0.2, color = "gray") +  # Add all points
  coord_cartesian(ylim = c(-3, 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  labs(title = "Connectivity after PROTAC or catalytic inhibitor treatment",
       y = "Change in connectivity (t-value)") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(
    "UNC1999" = "skyblue",
    "PROTAC" = "lightcoral"
  ))
dev.off()

p_vals <- c(
  wilcox.test(df_long$stat[df_long$Treatment == "PROTAC"], mu = 0)$p.value, # p-value = 2.613694e-17
  wilcox.test(df_long$stat[df_long$Treatment == "UNC1999"], mu = 0)$p.value) # p-value = 0.3753731

median(df_long$stat[df_long$Treatment == "PROTAC"]) # -0.1332606
median(df_long$stat[df_long$Treatment == "UNC1999"]) # -0.001756495

# Wilcoxon rank sum test with continuity correction
wilcox.test(df_long$stat[df_long$Treatment == "PROTAC"], df_long$stat[df_long$Treatment == "UNC1999"])$p.value
# 1.013874e-16

######################################################################################
### As a function of meanDMSO
######################################################################################

###########################################
### Visualisation: changes in connectivity as a function of meanDMSO
### A plot to compare PROTAC and UNC1999 per DMOSmean_bin
###########################################
# Combine both datasets in long format
df_long <- bind_rows(
  PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal_unique %>%
    filter(DMSO_mean_bin != "") %>%
    select(DMSO_mean, DMSO_mean_bin, stat = stat_res_PROTAC_vs_DMSO_merged) %>%
    mutate(Treatment = "PROTAC"),
  PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal_unique %>%
    filter(DMSO_mean_bin != "") %>%
    select(DMSO_mean, DMSO_mean_bin, stat = stat_res_UNC1999_vs_DMSO_merged) %>%
    mutate(Treatment = "UNC1999")
)

# Compute correlations per treatment (optional, to annotate)
cor_values <- df_long %>%
  group_by(Treatment) %>%
  summarise(cor_value = cor(stat, as.numeric(as.character(DMSO_mean)), use = "complete.obs"))
# Treatment cor_value
# <chr>         <dbl>
# 1 PROTAC      -0.174 
# 2 UNC1999     -0.0725

# Plot
pdf(file.path(working_dir, "Connectivity_UNC1999_PROTAC_per_DMSO_mean_bin_boxplot.pdf"), width = 6, height = 4)
ggplot(df_long, aes(x = DMSO_mean_bin, y = stat, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(
    title = "Changes in connectivity by DMSO bins",
    x = "DMSO mean bin",
    y = "Change in connectivity (t-value)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("UNC1999" = "skyblue", "PROTAC" = "lightcoral")) +
  theme(legend.position = "top")
dev.off()

# Stat test: compare all to 0
wilcox_results <- df_long %>%
  group_by(DMSO_mean_bin, Treatment) %>%
  summarise(
    p_value = wilcox.test(stat, mu = 0)$p.value,
    median_stat = median(stat, na.rm = TRUE),
    n = sum(!is.na(stat))
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

wilcox_results

# Compare PROTAC vs UNC1999 in each DMSO_mean_bin. Test: Wilcoxon rank-sum test (non-parametric two-group comparison)
treatment_comparison <- df_long %>%
  group_by(DMSO_mean_bin) %>%
  summarise(
    p_value = wilcox.test(stat[Treatment == "PROTAC"],
                          stat[Treatment == "UNC1999"])$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

treatment_comparison


###########################################
### Visualisation: changes in CUT&Tag stat as a function of meanDMSO
###########################################

# Combine both datasets in long format
df_longCnT <- bind_rows(
  PECHiC_deseq_stats_clusters_res_PROTAC_DMSO_CnT_all_dmsoBin_normal_unique %>%
    filter(DMSO_mean_bin != "") %>%
    select(DMSO_mean, DMSO_mean_bin, stat, log2FoldChange) %>%
    mutate(Treatment = "PROTAC"),
  PECHiC_deseq_stats_clusters_res_UNC1999_DMSO_CnT_all_dmsoBin_normal_unique %>%
    filter(DMSO_mean_bin != "") %>%
    select(DMSO_mean, DMSO_mean_bin, stat, log2FoldChange) %>%
    mutate(Treatment = "UNC1999")
)

# Compute correlations per treatment (optional, to annotate)
cor_values <- df_longCnT %>%
  group_by(Treatment) %>%
  summarise(cor_value = cor(stat, as.numeric(as.character(DMSO_mean)), use = "complete.obs"))
# Treatment cor_value
# <chr>         <dbl>
#   1 PROTAC    -0.390
# 2 UNC1999     -0.432


# Plot stat per DMSO_mean bins
pdf(file.path(working_dir, "Stat_CUTnTag_UNC1999_PROTAC_per_DMSO_mean_bin_boxplot.pdf"), width = 6, height = 4)
ggplot(df_longCnT, aes(x = DMSO_mean_bin, y = stat, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, position = position_dodge(0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  coord_cartesian(ylim = c(-10, 3)) +
  labs(
    title = "Changes in stat H3K27me3 by DMSO bins",
    x = "DMSO mean bin",
    y = "stat CUT&Tag"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("UNC1999" = "skyblue", "PROTAC" = "lightcoral")) +
  theme(legend.position = "top")
dev.off()

# Stat test: compare all to 0
wilcox_results <- df_longCnT %>%
  group_by(DMSO_mean_bin, Treatment) %>%
  summarise(
    p_value = wilcox.test(stat, mu = 0)$p.value,
    median_stat = median(stat, na.rm = TRUE),
    n = sum(!is.na(stat))
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

wilcox_results

# Compare PROTAC vs UNC1999 in each DMSO_mean_bin. Test: Wilcoxon rank-sum test (non-parametric two-group comparison)
treatment_comparison <- df_longCnT %>%
  group_by(DMSO_mean_bin) %>%
  summarise(
    p_value = wilcox.test(stat[Treatment == "PROTAC"],
                          stat[Treatment == "UNC1999"])$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

treatment_comparison

