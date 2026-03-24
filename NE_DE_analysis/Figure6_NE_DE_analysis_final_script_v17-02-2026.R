######################################################################################################
### Figure 6 and Supplemental Figure S4 - PE contacts persist after developmental activation
### Author: Dr. Marina Nocente and Dr. Mikhail Spivakov
###
### This script reproduces some panels of Figure 6 and Supplemental Figure S4
###  Supplemental Figure S4A - Number of poised enhancers (PEs) and contacted other ends (promoters or other PEs) 
###                            for different chromatin state dynamics upon differentiation of hPSCs into definitive endoderm (DE)
###
###  Supplemental Figure S4B - Number of poised enhancers (PEs) and contacted other ends (promoters or other PEs) 
###                            for different chromatin state dynamics upon differentiation of hPSCs into neuroectoderm (NE)
###
###  Figure 6B - Log2 fold change (Log2FC) in contact strength (proxied by arcsinh-CHiCAGO scores) for PEs involved  
###              in different types of contact dynamics that became active during differentiation to DE
###
###  Figure 6C - Log2 fold change (Log2FC) in contact strength (proxied by arcsinh-CHiCAGO scores) for PEs involved
###              in different types of contact dynamics that became active during differentiation to NE
### 
###  Supplemental Figure S4C: Representative examples of PE contacts showing different temporal dynamics 
###                         upon the naive-to-primed transition on day 14 of the transition and after differentiation into DE
###
###  Supplemental Figure S4D: Representative examples of PE contacts showing different temporal dynamics 
###                         upon the naive-to-primed transition on day 14 of the transition and after differentiation into NE
######################################################################################################

#######################################
## Libraries
#######################################
library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(ppclust)

## Graphical representation
library(ggplot2)
library(scales)
library(plotgardener)
library(grid)


#######################################
## File paths
#######################################
base_dir <- file.path("~/Documents/Bioinformatics/NE_DE_analysis")

# CUT&Tag pileup
PileUp_273_file <- file.path("~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/PileUp_273_final.tsv")
PileUp_271_file <- file.path("~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/PileUp_271_final.tsv")
PileUp_41_file <- file.path("~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/PileUp_41_final.tsv")
PileUp_43_file <- file.path("~/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/PileUp_43_final.tsv")

# CUT&Tag bigwig directory
CUTnTag_files_dir <- file.path("~/Documents/Bioinformatics/CnT_transition/bigwig_merged_just_sorted_indexedBamFiles/")

# Focused peakmatrix CHiCAGO scores and temporal clusters
peakmatrixFull_file <- file.path("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
Focused_peakmatrix_file <- file.path("~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Focused_peakmatrix_PE_TSS_interaction.tsv")
cluster_file <- file.path("~/Documents/Bioinformatics/Final_scripts_paper/PoisedEnhancers/Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt")

# Coordinates of PEs
baitmap_file <- file.path("~/Documents/Bioinformatics/Network_plot/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt")


#######################################
## Prepare the CUT&Tag data for DESeq2 analysis between Day14 and DE or NE differentiation
#######################################
# Peakmatrix with contact score using arcsinh-CHiCAGO scores
pm <- fread(Focused_peakmatrix_file)
head(pm)
dim(pm) # 2712    12

# Contact classes
cluster <- fread(cluster_file)
head(cluster)

# Add the contact classes to the peakmatrix
pm <- merge(
  pm,
  cluster[, .(baitID, oeID, merged_dtw_trans_cl6_2with6)],
  by = c("baitID", "oeID"),
  all.x = TRUE
)
dim(pm) # 12712    13

# Add baits and otherEnds information
peakmatrixFull <- fread(peakmatrixFull_file)
pm <- pm %>%
  left_join(
    peakmatrixFull %>%
      select(baitID, oeID,
             baitChr, baitStart, baitEnd, baitName,
             oeChr, oeStart, oeEnd, oeName),
    by = c("baitID", "oeID")
  )
head(pm)
dim(pm)
fwrite(pm, file.path(base_dir, "Peakmatrix_focused_transition_NE_DE_TemporalClusters.tx"), sep = "\t", quote=F, col.names=T)


### DESeq2 analysis for CUT&Tag data with different histone marks across the differentiation
## H3K27me3
PileUp_273 <- fread(PileUp_273_file)
counts273 <- as.matrix(PileUp_273[,c("E14_273", "C14_273", "D14_273",
                                     "DE1_273", "DE2_273" ,"DE3_273",
                                     "NE1_273", "NE2_273", "NE3_273")])
rownames(counts273) <- PileUp_273$ID

coldata_273 <- data.frame(condition = rep("H3K27me3", ncol(counts273)) , # Add a column "condition": H3K27me3 for the 30 lignes
                      type = factor(rep(c("day14", "DE", "NE"), each = 3))) # Add a column "type": time course for the 30 lignes

row.names(coldata_273) <- colnames(counts273) # Name the rows of coldata as the columns of counts273 matrix
dds273 <- DESeqDataSetFromMatrix(countData = counts273,
                                 colData = coldata_273,
                                 design = ~ type)

dds273 <- DESeq(dds273)

# H3K27me3 DESeq2 results: DE versus day14
DEvsDay14_273 = results(dds273, contrast = c("type", "DE", "day14"))
nrow(DEvsDay14_273[which(DEvsDay14_273$padj<0.05),]) # 35467 / 577983

# H3K27me3 DESeq2 results: NE versus day14
NEvsDay14_273 = results(dds273, contrast = c("type", "NE", "day14"))
nrow(NEvsDay14_273[which(NEvsDay14_273$padj<0.05),]) # 62307 / 577983

# Convert to data frame
DEvsDay14_273_df <- as.data.frame(DEvsDay14_273)
DEvsDay14_273_df$fragID <- rownames(DEvsDay14_273_df)
head(DEvsDay14_273_df)

NEvsDay14_273_df <- as.data.frame(NEvsDay14_273)
NEvsDay14_273_df$fragID <- rownames(NEvsDay14_273_df)
head(NEvsDay14_273_df)

fwrite(DEvsDay14_273_df, file.path(base_dir, "deseq_DEvsDay14_273_vf.txt"), sep = "\t", quote=F, col.names=T)
fwrite(NEvsDay14_273_df, file.path(base_dir, "deseq_NEvsDay14_273_vf.txt"), sep = "\t", quote=F, col.names=T)

## H3K27ac
PileUp_271 <- fread(PileUp_271_file)
counts271 <- as.matrix(PileUp_271[,c("B14_271", "C14_271", "D14_271",
                                     "DE1_271", "DE2_271" ,"DE3_271",
                                     "NE1_271", "NE2_271", "NE3_271")])
rownames(counts273) <- PileUp_273$ID
coldata_271 <- data.frame(condition = rep("H3K27ac", ncol(counts271)) , 
                      type = factor(rep(c("day14", "DE", "NE"), each = 3))) 
row.names(coldata_271) <- colnames(counts271) 

dds271 <- DESeqDataSetFromMatrix(countData = counts271,
                                 colData = coldata_271,
                                 design = ~ type)
dds271 <- DESeq(dds271)

# H3K27ac DESeq2 results: DE versus day14
DEvsDay14_271 = results(dds271, contrast = c("type", "DE", "day14"))
nrow(DEvsDay14_271[which(DEvsDay14_271$padj<0.05),]) # 38004 / 577983

# H3K27ac DESeq2 results: NE versus day14
NEvsDay14_271 = results(dds271, contrast = c("type", "NE", "day14"))
nrow(NEvsDay14_271[which(NEvsDay14_271$padj<0.05),]) # 39922 / 577983

# Convert to data frame
DEvsDay14_271_df <- as.data.frame(DEvsDay14_271)
DEvsDay14_271_df$fragID <- rownames(DEvsDay14_271_df)
head(DEvsDay14_271_df)

NEvsDay14_271_df <- as.data.frame(NEvsDay14_271)
NEvsDay14_271_df$fragID <- rownames(NEvsDay14_271_df)
head(NEvsDay14_271_df)

fwrite(DEvsDay14_271_df, file.path(base_dir, "deseq_DEvsDay14_271_vf.txt"), sep = "\t", quote=F, col.names=T)
fwrite(NEvsDay14_271_df, file.path(base_dir, "deseq_NEvsDay14_271_vf.txt"), sep = "\t", quote=F, col.names=T)


## H3K4me1
PileUp_41 <- fread(PileUp_41_file)
counts41 <- as.matrix(PileUp_41[,c("B14_41", "C14_41", "D14_41",
                                     "DE1_41", "DE2_41" ,"DE3_41",
                                     "NE1_41", "NE2_41", "NE3_41")])
rownames(counts41) <- PileUp_41$ID
coldata_41 <- data.frame(condition = rep("H3K4me1", ncol(counts41)) , 
                          type = factor(rep(c("day14", "DE", "NE"), each = 3))) 
row.names(coldata_41) <- colnames(counts41) 

dds41 <- DESeqDataSetFromMatrix(countData = counts41,
                                 colData = coldata_41,
                                 design = ~ type)
dds41 <- DESeq(dds41)

# H3K4me1 DESeq2 results: DE versus day14
DEvsDay14_41 = results(dds41, contrast = c("type", "DE", "day14"))
nrow(DEvsDay14_41[which(DEvsDay14_41$padj<0.05),]) # 123233 / 577983

# H3K4me1 DESeq2 results: NE versus day14
NEvsDay14_41 = results(dds41, contrast = c("type", "NE", "day14"))
nrow(NEvsDay14_41[which(NEvsDay14_41$padj<0.05),]) # 106382 / 577983

# Convert to data frame
DEvsDay14_41_df <- as.data.frame(DEvsDay14_41)
DEvsDay14_41_df$fragID <- rownames(DEvsDay14_41_df)
head(DEvsDay14_41_df)

NEvsDay14_41_df <- as.data.frame(NEvsDay14_41)
NEvsDay14_41_df$fragID <- rownames(NEvsDay14_41_df)
head(NEvsDay14_41_df)

fwrite(DEvsDay14_41_df, file.path(base_dir, "deseq_DEvsDay14_41_vf.txt"), sep = "\t", quote=F, col.names=T)
fwrite(NEvsDay14_41_df, file.path(base_dir, "deseq_NEvsDay14_41_vf.txt"), sep = "\t", quote=F, col.names=T)


## H3K4me3
PileUp_43 <- fread(PileUp_43_file)
counts43 <- as.matrix(PileUp_43[,c("B14_43", "C14_43", "D14_43",
                                   "DE1_43", "DE2_43" ,"DE3_43",
                                   "NE1_43", "NE2_43", "NE3_43")])
rownames(counts43) <- PileUp_43$ID
coldata_43 <- data.frame(condition = rep("H3K4me3", ncol(counts43)) , # Add a column "condition": H3K4me3 for the 30 lignes
                         type = factor(rep(c("day14", "DE", "NE"), each = 3))) # Add a column "type": time course for the 30 lignes
row.names(coldata_43) <- colnames(counts43) # Name the rows of coldata as the columns of counts273 matrix

dds43 <- DESeqDataSetFromMatrix(countData = counts43,
                                colData = coldata_43,
                                design = ~ type)
dds43 <- DESeq(dds43)

# H3K4me3 DESeq2 results: DE versus day14
DEvsDay14_43 = results(dds43, contrast = c("type", "DE", "day14"))
nrow(DEvsDay14_43[which(DEvsDay14_43$padj<0.05),]) # 12804 / 577983

# H3K4me1 DESeq2 results: NE versus day14
NEvsDay14_43 = results(dds43, contrast = c("type", "NE", "day14"))
nrow(NEvsDay14_43[which(NEvsDay14_43$padj<0.05),]) # 15071 / 577983

# Convert to data frame
DEvsDay14_43_df <- as.data.frame(DEvsDay14_43)
DEvsDay14_43_df$fragID <- rownames(DEvsDay14_43_df)
head(DEvsDay14_43_df)

NEvsDay14_43_df <- as.data.frame(NEvsDay14_43)
NEvsDay14_43_df$fragID <- rownames(NEvsDay14_43_df)
head(NEvsDay14_43_df)

fwrite(DEvsDay14_43_df, file.path(base_dir, "deseq_DEvsDay14_43_vf.txt"), sep = "\t", quote=F, col.names=T)
fwrite(NEvsDay14_43_df, file.path(base_dir, "deseq_NEvsDay14_43_vf.txt"), sep = "\t", quote=F, col.names=T)


### Generate a Supplemental Table
## Define DESeq2 results
# H3K27ac
DE_H3K27ac <- DEvsDay14_271_df
NE_H3K27ac <- NEvsDay14_271_df

# H3K27me3
DE_H3K27me3 <- DEvsDay14_273_df
NE_H3K27me3 <- NEvsDay14_273_df

# H3K4me1
DE_H3K4me1 <- DEvsDay14_41_df
NE_H3K4me1 <- NEvsDay14_41_df

## Function to prepare supplementary table
prepare_supp_table <- function(DE_df, NE_df, mark_name) {
  
  # Make IDs character
  DE_df$fragID <- as.character(DE_df$fragID)
  NE_df$fragID <- as.character(NE_df$fragID)
  pm$baitID <- as.character(pm$baitID)
  
  # Merge DE results
  merged_df <- merge(pm,
                     DE_df[, c("fragID","log2FoldChange","pvalue","padj")],
                     by.x="baitID", by.y="fragID", all.x=TRUE)
  
  # Merge NE results
  merged_df <- merge(merged_df,
                     NE_df[, c("fragID","log2FoldChange","pvalue","padj")],
                     by.x="baitID", by.y="fragID", all.x=TRUE,
                     suffixes = c("_DE","_NE"))
  
  # Rename columns clearly
  colnames(merged_df)[which(colnames(merged_df) == "log2FoldChange_DE")] <- "log2FC_DE_day14"
  colnames(merged_df)[which(colnames(merged_df) == "pvalue_DE")] <- "pvalue_DE_day14"
  colnames(merged_df)[which(colnames(merged_df) == "padj_DE")] <- "padj_DE_day14"
  colnames(merged_df)[which(colnames(merged_df) == "log2FoldChange_NE")] <- "log2FC_NE_day14"
  colnames(merged_df)[which(colnames(merged_df) == "pvalue_NE")] <- "pvalue_NE_day14"
  colnames(merged_df)[which(colnames(merged_df) == "padj_NE")] <- "padj_NE_day14"
  
  # Reorder columns
  merged_df <- merged_df %>%
    select(baitID, oeID, day14, NE, DE, merged_dtw_trans_cl6_2with6,
           log2FC_DE_day14, pvalue_DE_day14, padj_DE_day14,
           log2FC_NE_day14, pvalue_NE_day14, padj_NE_day14)
  
  # Export
  out_file <- paste0("Supplementary_Table_", mark_name, ".tsv")
  fwrite(merged_df, file.path(base_dir, out_file), sep="\t", quote=F, col.names=T)
  message(paste("Table for", mark_name, "saved to", out_file))
}


## Generate tables for all marks
prepare_supp_table(DE_H3K27ac, NE_H3K27ac, "H3K27ac")
prepare_supp_table(DE_H3K27me3, NE_H3K27me3, "H3K27me3")
prepare_supp_table(DE_H3K4me1, NE_H3K4me1, "H3K4me1")


#######################################
# Generate all possible combinations of histone mark change states
# Each mark can be -1 (loss), 0 (not significant), or 1 (gain)
# For 3 marks → 3^3 = 27 combinations
#######################################
pepm <- pm

combs_sign <- as.matrix(expand.grid(rep(list(-1:1), 3)))
colnames(combs_sign)<- c("K27me3", "K4me1", "K27ac")

#######################################
# Initialize vectors to store counts for DE condition
# One entry per histone-mark combination (27 total)
#######################################

DEcounts_sign = vector("numeric", length=nrow(combs_sign))
DEcounts_baits_sign = vector("numeric", length=nrow(combs_sign))
DEcounts_oe_sign = vector("numeric", length=nrow(combs_sign))
DEcounts_oe_both = vector("numeric", length=nrow(combs_sign))

#######################################
# Loop over each histone-mark combination (DE case)
#######################################
# Build a logical expression dynamically depending on:
# - significance (pvalue < 0.05 if -1 or 1, >0.05 if 0)
# - direction of log2FoldChange (negative if -1, positive if 1)
# If 0 (not significant), fold change is unconstrained (<1e10 used as permissive condition)

for(i in 1:nrow(combs_sign)){
  # if not significant leaving the condition "XX$log2FoldChange<1e10" like this as a permissive plug as none of them are exactly 0
  expr_str <- paste0(
    "which(DEvsDay14_273$pvalue ", ifelse(combs_sign[i, 1], "<", ">"), " 0.05 & ",
    "DEvsDay14_41$pvalue ", ifelse(combs_sign[i, 2], "<", ">"), " 0.05 & ",
    "DEvsDay14_271$pvalue ", ifelse(combs_sign[i, 3], "<", ">"), " 0.05 & ",
    
    "DEvsDay14_273$log2FoldChange ",
    if (combs_sign[i, 1]) {
      if (combs_sign[i, 1] < 0) "<0" else ">0"
    } else {
      "<1e10"
    }, " & ",
    
    "DEvsDay14_41$log2FoldChange ",
    if (combs_sign[i, 2]) {
      if (combs_sign[i, 2] < 0) "<0" else ">0"
    } else {
      "<1e10"
    }, " & ",
    
    "DEvsDay14_271$log2FoldChange ",
    if (combs_sign[i, 3]) {
      if (combs_sign[i, 3] < 0) "<0" else ">0"
    } else {
      "<1e10"
    },
    ")"
  )
  
  res = eval(parse(text = expr_str))
  posIDs = PileUp_271$ID[res]
  
  DEcounts_sign[i] = length(posIDs)
  DEcounts_baits_sign[i] = length(posIDs[posIDs %in% pepm$baitID])
  DEcounts_oe_sign[i] = length(posIDs[posIDs %in% pepm$oeID])
  DEcounts_oe_both[i] = nrow(pepm[baitID%in%posIDs & oeID%in%posIDs])
  
}

# Combine combinations and DE counts into one matrix
DEcounts_sign_pepm = cbind(combs_sign, DEcounts_sign, DEcounts_baits_sign, DEcounts_oe_sign, DEcounts_oe_both)
DEcounts_sign_pepm

#######################################
# Define biologically meaningful subsets (DE case)
#######################################
unpoised = which(DEvsDay14_273$pvalue<0.05 & DEvsDay14_273$log2FoldChange<0 & DEvsDay14_271$pvalue<0.05 & DEvsDay14_271$log2FoldChange>0) 
gainK4me3 = which(DEvsDay14_43$pvalue<0.05 & DEvsDay14_43$log2FoldChange>0)
gainK4me1 = which(DEvsDay14_41$pvalue<0.05 & DEvsDay14_41$log2FoldChange>0)
lossK27me3 = which(DEvsDay14_273$pvalue<0.05 & DEvsDay14_273$log2FoldChange<0)

posIDs = PileUp_271$ID[unpoised]
length(posIDs[posIDs %in% pepm$baitID]) # 307
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[gainK4me3]]) # 33
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[c(gainK4me3, gainK4me1)]]) # 120
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[lossK27me3]]) # 173


#######################################
# Repeat the exact same combinatorial analysis for NE condition
#######################################

NEcounts_sign = vector("numeric", length=nrow(combs_sign))
NEcounts_baits_sign = vector("numeric", length=nrow(combs_sign))
NEcounts_oe_sign = vector("numeric", length=nrow(combs_sign))
NEcounts_oe_both = vector("numeric", length=nrow(combs_sign))

for(i in 1:nrow(combs_sign)){
  # if not significant leaving the condition "XX$log2FoldChange<1e10" like this as a permissive plug as none of them are exactly 0
  expr_str <- paste0(
    "which(NEvsDay14_273$pvalue ", ifelse(combs_sign[i, 1], "<", ">"), " 0.05 & ",
    "NEvsDay14_41$pvalue ", ifelse(combs_sign[i, 2], "<", ">"), " 0.05 & ",
    "NEvsDay14_271$pvalue ", ifelse(combs_sign[i, 3], "<", ">"), " 0.05 & ",
    
    "NEvsDay14_273$log2FoldChange ",
    if (combs_sign[i, 1]) {
      if (combs_sign[i, 1] < 0) "<0" else ">0"
    } else {
      "<1e10"
    }, " & ",
    
    "NEvsDay14_41$log2FoldChange ",
    if (combs_sign[i, 2]) {
      if (combs_sign[i, 2] < 0) "<0" else ">0"
    } else {
      "<1e10"
    }, " & ",
    
    "NEvsDay14_271$log2FoldChange ",
    if (combs_sign[i, 3]) {
      if (combs_sign[i, 3] < 0) "<0" else ">0"
    } else {
      "<1e10"
    },
    ")"
  )
  
  res = eval(parse(text = expr_str))
  posIDs = PileUp_271$ID[res]
  
  NEcounts_sign[i] = length(posIDs)
  NEcounts_baits_sign[i] = length(posIDs[posIDs %in% pepm$baitID])
  NEcounts_oe_sign[i] = length(posIDs[posIDs %in% pepm$oeID])
  NEcounts_oe_both[i] = nrow(pepm[baitID%in%posIDs & oeID%in%posIDs])
  
  
}

unpoised = which(NEvsDay14_273$pvalue<0.05 & NEvsDay14_273$log2FoldChange<0 & NEvsDay14_271$pvalue<0.05 &
                   NEvsDay14_271$log2FoldChange>0) 
gainK4me3 = which(NEvsDay14_43$pvalue<0.05 & NEvsDay14_43$log2FoldChange>0)
gainK4me1 = which(NEvsDay14_41$pvalue<0.05 & NEvsDay14_41$log2FoldChange>0)
lossK27me3 = which(NEvsDay14_273$pvalue<0.05 & NEvsDay14_273$log2FoldChange<0)
posIDs = PileUp_271$ID[unpoised]
length(posIDs[posIDs %in% pepm$baitID]) # 99
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[c(gainK4me3, gainK4me1)]]) # 54
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[lossK27me3]]) # 37
nrow(pepm[baitID%in%PileUp_271$ID[unpoised] & oeID%in%PileUp_271$ID[lossK27me3] &  oeID%in%PileUp_271$ID[gainK4me3]]) # 15


#######################################
# Combine DE and NE counts
#######################################
counts_sign_pepm = cbind(DEcounts_sign_pepm, NEcounts_sign, NEcounts_baits_sign, NEcounts_oe_sign, NEcounts_oe_both)
head(counts_sign_pepm)

counts_sign_pepm_annot <- as.data.frame(counts_sign_pepm)
fwrite(counts_sign_pepm_annot, "counts_sign_pepm_annot.csv")


#######################################
# Annotate chromatin state transitions into biological categories
#######################################
counts_sign_pepm_annot <- counts_sign_pepm_annot %>%
  mutate(State = case_when(
    K27me3 == 0 & K4me1 == 0 & K27ac == 0 ~ "remain poised",
    K27me3 == -1 & K4me1 == -1 & K27ac!=1 ~ "become neutral",
    K27me3 == -1 & K27ac == 1 ~ "become active",
    K4me1 == -1 & K27ac!=1 ~ "become inactive",
    TRUE ~ NA_character_
  ))

# Keep only biologically annotated states
counts_sign_pepm_annot_filtered <- counts_sign_pepm_annot %>% filter(!is.na(State))

# Aggregate counts per biological state
counts_sign_pepm_summary <- counts_sign_pepm_annot_filtered %>% group_by(State) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE))


#######################################
# Supplemental Figure S4A
#  Number of poised enhancers (PEs) (blue bars) and contacted other ends (promoters or other PEs; orange bars) 
#  for different chromatin state dynamics upon differentiation of hPSCs into definitive endoderm (DE)
#######################################
plot_data_DE <- counts_sign_pepm_summary %>%
  select(State, DEcounts_baits_sign, DEcounts_oe_sign) %>%
  pivot_longer(
    cols = c(DEcounts_baits_sign, DEcounts_oe_sign),
    names_to = "Type",
    values_to = "Count"
  )

my_colors <- c(
  "DEcounts_baits_sign" = "#0072B2",  # blue
  "DEcounts_oe_sign"    = "#D55E00"   # reddish
)

# Plot
pdf(file.path(base_dir, "SuppFigureS4A_Barplot_state_PE_after_DE_differentiation.pdf"))
ggplot(plot_data_DE, aes(x = State, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "DEcounts_baits_sign" = "#0072B2",
      "DEcounts_oe_sign"    = "#D55E00"
    ),
    labels = c("At the PE baits", "At the other ends")
  ) +
  labs(
    x = "PE state transition in DE",
    y = "Number of PEs per state after DE differentiation",
    fill = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
dev.off()

#######################################
# Supplemental Figure S4B
#  Number of poised enhancers (PEs) (blue bars) and contacted other ends (promoters or other PEs; orange bars) 
#  for different chromatin state dynamics upon differentiation of hPSCs into neuroectoderm (NE)
#######################################
plot_data_NE <- counts_sign_pepm_summary %>%
  select(State, NEcounts_baits_sign, NEcounts_oe_sign) %>%
  pivot_longer(
    cols = c(NEcounts_baits_sign, NEcounts_oe_sign),
    names_to = "Type",
    values_to = "Count"
  )

my_colors <- c(
  "NEcounts_baits_sign" = "#0072B2",  # blue
  "NEcounts_oe_sign"    = "#D55E00"   # reddish
)

# Plot
pdf(file.path(base_dir, "SuppFigureS4B_Barplot_state_PE_after_NE_differentiation.pdf"))
ggplot(plot_data_NE, aes(x = State, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "NEcounts_baits_sign" = "#0072B2",
      "NEcounts_oe_sign"    = "#D55E00"
    ),
    labels = c("At the PE baits", "At the other ends")
  ) +
  labs(
    x = "PE state transition in NE",
    y = "Number of PEs per state after NE differentiation",
    fill = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
dev.off()


#######################################
### Identifying poised chromatin regions that become active during differentiation (DE and NE)
#######################################
DEpcToActive = PileUp_271[which(DEvsDay14_273$pvalue<0.05 & DEvsDay14_273$log2FoldChange<0 &
                                  DEvsDay14_271$pvalue<0.05 &DEvsDay14_271$log2FoldChange>0)]$ID # 4285

NEpcToActive = PileUp_271[which(NEvsDay14_273$pvalue<0.05 & NEvsDay14_273$log2FoldChange<0 &
                                   NEvsDay14_271$pvalue<0.05 &NEvsDay14_271$log2FoldChange>0)]$ID # 9228
 

#######################################
### Figure 6B - Log2 fold change (Log2FC) in contact strength (proxied by arcsinh-CHiCAGO scores) 
###  for PEs involved in different types of contact dynamics that became active during differentiation to definitive endoderm (DE)
#######################################

### Temporal contact clusters
scores_TempClust <- pm[,c("baitID", "oeID", "day14","NE", "DE", "merged_dtw_trans_cl6_2with6")]
scores_TempClust[is.na(merged_dtw_trans_cl6_2with6), merged_dtw_trans_cl6_2with6 := "New"]

## For enhancers which become active after DE differentiation
DEpcToActive_new_class = list(Constant = log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Constant"]$DE + 0.1) - log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Constant"]$day14 + 0.1),
                              Gained_early = log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Gained early"]$DE + 0.1) - log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Gained early"]$day14 + 0.1),
                              Gained_late = log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Gained late"]$DE + 0.1) - log2(scores_TempClust[baitID %in% DEpcToActive & merged_dtw_trans_cl6_2with6=="Gained late"]$day14 + 0.1),
                              Other = log2(scores_TempClust[baitID %in% DEpcToActive & !merged_dtw_trans_cl6_2with6 %in% c("Gained early","Gained late","Constant")]$DE + 0.1) - log2(scores_TempClust[baitID %in% DEpcToActive & !merged_dtw_trans_cl6_2with6 %in% c("Gained early","Gained late","Constant")]$day14 + 0.1)
                              )
## Plot
pdf(file.path(base_dir, "Figure6B_DEpcToActive_boxplot.pdf"), width = 6, height = 5)
boxplot(DEpcToActive_new_class, 
        ylim = c(-5, 5), 
        main = "Poised enhancers activated during\n differentiation to definitive endoderm",
        ylab = "Log2FC in contact strength")
dev.off()

kruskal.test(DEpcToActive_new_class) # p-value < 2.2e-16

length(unique(scores_TempClust[baitID %in% DEpcToActive]$baitID)) # 307


#######################################
### Figure 6C - Log2 fold change (Log2FC) in contact strength (proxied by arcsinh-CHiCAGO scores) 
###  for PEs involved in different types of contact dynamics that became active during differentiation to neuroectoderm (NE)
#######################################

## For enhancers which become active after NE differentiation
NEpcToActive_new_class = list(Constant = log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Constant"]$NE + 0.1) - log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Constant"]$day14 + 0.1),
                              Gained_early = log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Gained early"]$NE + 0.1) - log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Gained early"]$day14 + 0.1),
                              Gained_late = log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Gained late"]$NE + 0.1) - log2(scores_TempClust[baitID %in% NEpcToActive & merged_dtw_trans_cl6_2with6=="Gained late"]$day14 + 0.1),
                              Other = log2(scores_TempClust[baitID %in% NEpcToActive & !merged_dtw_trans_cl6_2with6 %in% c("Gained early","Gained late","Constant")]$NE + 0.1) - log2(scores_TempClust[baitID %in% NEpcToActive & !merged_dtw_trans_cl6_2with6 %in% c("Gained early","Gained late","Constant")]$day14 + 0.1)
                              )

pdf(file.path(base_dir, "Figure6C_NEpcToActive_boxplot.pdf"), width = 6, height = 5)
boxplot(NEpcToActive_new_class, 
        ylim = c(-5, 5), 
        main = "Poised enhancers activated during\n differentiation to neuroectoderm",
        ylab = "Log2FC in contact strength")
dev.off()

kruskal.test(NEpcToActive_new_class) # p-value = 1.419e-07

length(unique(scores_TempClust[baitID %in% NEpcToActive]$baitID)) # 99


#######################################
# Baitmap file
#######################################
baitmap <- fread(baitmap_file)
setnames(baitmap, c("V1", "V2", "V3", "V4"), c("chrom", "start", "end", "baitID"))
baitmap[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]

#######################################
### Supplemental Figure 4C
###  Representative examples of PE contacts showing different temporal dynamics upon the naive-to-primed transition
###  on day 14 of the transition and after differentiation into DE
#######################################

#######################################
## Common parameters
#######################################
#### CUT&Tagg bigWig files
k27me3_files <- c(
  day14="merged_Day14_273.bam.bw",
  DE="merged_DE_273.bam.bw"
)

k4me1_files <- c(
  day14="merged_Day14_41.bam.bw",
  DE="merged_DE_41.bam.bw"
)

k27ac_files <- c(
  day14="merged_Day14_271.bam.bw",
  DE="merged_DE_271.bam.bw"
)

#### Timepoints
timepoints_labels  <- c("Day14", "DE")
timepoints_columns <- c("day14", "DE")

#######################################
#### DE becoming Active - CONSTANT contact of interest example
#######################################
#### Keep just this interaction: baitID=4034
# baitID = 4034 
# oeID = 4028 

scores_TempClust_bait4034 <- scores_TempClust[baitID == 4034 & oeID == 4028]
head(scores_TempClust_bait4034)

Peakmatrix <- pm
scores_TempClust_bait4034 <- merge(scores_TempClust_bait4034, 
                                   Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))


#### Interaction and score for each timepoint 
hic_df_day14_bait4034 <- scores_TempClust_bait4034[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]

hic_df_NE_bait4034 <- scores_TempClust_bait4034[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait4034 <- scores_TempClust_bait4034[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]

#### Selection of the contact of interest
bait_of_interest_4034 <- baitmap[
  chrom == "chr1" & 
    start == 18640186 & 
    end == 18640453
]

#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait4034,
  hic_df_DE_bait4034
)

names(hic_list) <- timepoints_labels

#### Parameters
params <- pgParams(
  chrom = "chr1", 
  chromstart = 18550000, chromend = 18800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "1", 
  chromstart = 18550000, chromend = 18800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#### Color: arcsinh + global scale
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
#    Day14       DE 
# 1.262336 1.261156 
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.272336

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


#### Plot PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4C_top_CnT_PECHiC_Day14toDE_Constant_arcsinhColor.pdf"), width = 6, height = 15)

pageCreate(width = 4, height = 15, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_4034, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)


## Layout constants
y_range_k27 <- c(0, 1000)
y_range_k4 <- c(0, 200)
y_range_k27ac <- c(0, 100)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9


## LOOP: ALL TIMEPOINTS

for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkgreen"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(14 * sf, "inches"), # note this is the *centre* position of the raster
  width = unit(0.2 * sf, "inches"),
  height = unit(5 * sf, "inches"),
  just = c("left", "center"),
  interpolate = TRUE
)

## label positions shrink too
breaks_labels <- round(seq(0, global_max_asinh, length.out = 7), 2)
ylab_positions <- seq((14-5/2) * sf, (14+5/2) * sf, length.out = length(breaks_labels))

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
  y = unit(14 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()


#######################################
#### DE becoming Active - GAINED LATE contact of interest example
#######################################
## Keep just this interaction: baitID=24129
# baitID = 24129 
# oeID = 24174

scores_TempClust_bait24129 <- scores_TempClust[baitID == 24129 & oeID == 24174]
head(scores_TempClust_bait24129)

scores_TempClust_bait24129 <- merge(scores_TempClust_bait24129, 
                                   Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))
head(scores_TempClust_bait24129)


#### Interaction and score for each timepoint 
hic_df_day14_bait24129 <- scores_TempClust_bait24129[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]

hic_df_NE_bait24129 <- scores_TempClust_bait24129[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait24129 <- scores_TempClust_bait24129[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]


#### Define genomic region and widths of all plots:
params <- pgParams(
  chrom = "chr1", 
  chromstart = 115700000, chromend = 116100000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "1", 
  chromstart = 115700000, chromend = 116100000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

##### Selection of the contact of interest
bait_of_interest_24129 <- baitmap[
  chrom == "chr1" & 
    start == 115826072 & 
    end == 115827569
]


#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait24129,
  hic_df_DE_bait24129
)

names(hic_list) <- timepoints_labels

#### Parameters
params <- pgParams(
  chrom = "chr1", 
  chromstart = 115700000, chromend = 116100000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "1", 
  chromstart = 115700000, chromend = 116100000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)


#### Color: arcsinh + global scale
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
#    Day14       DE 
# 1.5098929 0.5650333 
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.519893

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


#### Plot: PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4C_bottom_CnT_PECHiC_Day14toDE_GainedLate_arcsinhColor.pdf"), width = 6, height = 10)

pageCreate(width = 6, height = 10, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_24129, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
y_range_k27 <- c(0, 200)
y_range_k4 <- c(0, 200)
y_range_k27ac <- c(0, 200)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9

## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkorange"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(6 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)


dev.off()



#######################################
#### DE becoming Active - GAINED EARLY contact of interest example
#######################################
## Keep just this interaction: baitID=581761
# baitID = 581761 
# oeID = 581818 

scores_TempClust_bait581761 <- scores_TempClust[baitID == 581761 & oeID == 581818]
head(scores_TempClust_bait581761)

scores_TempClust_bait581761 <- merge(scores_TempClust_bait581761, 
                                    Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))
head(scores_TempClust_bait581761)


#### Interaction and score for each timepoint 
hic_df_day14_bait581761 <- scores_TempClust_bait581761[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]

hic_df_NE_bait581761 <- scores_TempClust_bait581761[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait581761 <- scores_TempClust_bait581761[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]


#### Selection of the contact of interest
bait_of_interest_581761 <- baitmap[
  chrom == "chrX" & 
    start == 137426685 & 
    end == 137427334
]

#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait581761,
  hic_df_DE_bait581761
)

names(hic_list) <- timepoints_labels


#### Parameters
params <- pgParams(
  chrom = "chrX", 
  chromstart = 137200000, chromend = 137700000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "X", 
  chromstart = 137200000, chromend = 137700000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)


#### Color: arcsinh + global scale
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
#    Day14       DE 
# 1.134010 1.264348  
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.274348

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


#### Plot: PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4C_middle_CnT_PECHiC_Day14toDE_GainedEarly_arcsinhColor.pdf"), width = 6, height = 10)

pageCreate(width = 6, height = 10, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_581761, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
y_range_k27 <- c(0, 500)
y_range_k4 <- c(0, 200)
y_range_k27ac <- c(0, 200)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9


## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkred"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(6 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()



#######################################
### Supplemental Figure 4D
###  Representative examples of PE contacts showing different temporal dynamics upon the naive-to-primed transition
###  on day 14 of the transition and after differentiation into NE
#######################################

#######################################
## Common parameters
#######################################
#### CUT&Tagg bigWig files
k27me3_files <- c(
  day14="merged_Day14_273.bam.bw",
  NE="merged_NE_273.bam.bw"
)

k4me1_files <- c(
  day14="merged_Day14_41.bam.bw",
  NE="merged_NE_41.bam.bw"
)

k27ac_files <- c(
  day14="merged_Day14_271.bam.bw",
  NE="merged_NE_271.bam.bw"
)

#### Timepoints
timepoints_labels  <- c("Day14", "NE")
timepoints_columns <- c("day14", "NE")


#######################################
#### NE becoming Active - GAINED EARLY contact of interest example
#######################################
## Keep just this interaction: baitID=103817
# baitID = 103817 
# oeID = 103761

scores_TempClust_bait103817 <- scores_TempClust[baitID == 103817 & oeID == 103761]
head(scores_TempClust_bait103817)

scores_TempClust_bait103817 <- merge(scores_TempClust_bait103817, 
                                     Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))
head(scores_TempClust_bait103817)


#### Interaction and score for each timepoint 
hic_df_day14_bait103817 <- scores_TempClust_bait103817[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]

hic_df_NE_bait103817 <- scores_TempClust_bait103817[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait103817 <- scores_TempClust_bait103817[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]


#### Selection of the contact of interest
bait_of_interest_103817 <- baitmap[
  chrom == "chr12" & 
    start == 5430195 & 
    end == 5430969
]


#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait103817,
  hic_df_NE_bait103817
)

names(hic_list) <- timepoints_labels


#### Parameters
params <- pgParams(
  chrom = "chr12", 
  chromstart = 5100000, chromend = 5500000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "12", 
  chromstart = 5100000, chromend = 5500000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)


#### Color: arcsinh + global scale
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
# Day14        NE 
# 0.9901653 1.7021311  
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.712131

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


#### Plot: PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4D_middle_CnT_PECHiC_Day14toNE_GainedEarly_arcsinhColor.pdf"), width = 6, height = 10)

pageCreate(width = 6, height = 10, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_103817, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
y_range_k27 <- c(0, 1000)
y_range_k4 <- c(0, 200)
y_range_k27ac <- c(0, 100)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9


## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkred"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(6 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()



#######################################
#### NE becoming Active - GAINED LATE contact of interest example
#######################################
## Keep just this interaction: baitID=174494
# baitID = 174494 
# oeID = 174392 

scores_TempClust_bait174494 <- scores_TempClust[baitID == 174494 & oeID == 174392]
head(scores_TempClust_bait174494)

scores_TempClust_bait174494 <- merge(scores_TempClust_bait174494, 
                                     Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))
head(scores_TempClust_bait174494)

#### Interaction and score for each timepoint 
hic_df_day14_bait174494 <- scores_TempClust_bait174494[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14 
)]

hic_df_NE_bait174494 <- scores_TempClust_bait174494[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait174494 <- scores_TempClust_bait174494[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]

#### Selection of the contact of interest
bait_of_interest_174494 <- baitmap[
  chrom == "chr15" & 
    start == 55742277 & 
    end == 55743902
]

#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait174494,
  hic_df_NE_bait174494
)

names(hic_list) <- timepoints_labels

#### Parameters
params <- pgParams(
  chrom = "chr15", 
  chromstart = 55200000, chromend = 55800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "15", 
  chromstart = 55200000, chromend = 55800000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#### Color: arcsinh + global scale
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
# Day14        NE 
# 1.584444 0.000000   
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.594444

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


#### Plot: PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4D_bottom_CnT_PECHiC_Day14toNE_GainedLate_arcsinhColor.pdf"), width = 6, height = 10)

pageCreate(width = 6, height = 10, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_174494, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
y_range_k27 <- c(0, 1200)
y_range_k4 <- c(0, 200)
y_range_k27ac <- c(0, 200)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9

## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkorange"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(6 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()



#######################################
#### NE becoming Active - CONSTANT contact of interest example
#######################################
## Keep just this interaction: baitID=159523
# baitID = 159523 
# oeID = 159562 

scores_TempClust_bait159523 <- scores_TempClust[baitID == 159523 & oeID == 159562]
head(scores_TempClust_bait159523)

scores_TempClust_bait159523 <- merge(scores_TempClust_bait159523, 
                                     Peakmatrix[,c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "oeID", "oeChr", "oeStart", "oeEnd", "oeName")], by=c("baitID", "oeID"))
head(scores_TempClust_bait159523)


#### Interaction and score for each timepoint 
hic_df_day14_bait159523 <- scores_TempClust_bait159523[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = day14  # Change to "day14" timepoint column
)]

hic_df_NE_bait159523 <- scores_TempClust_bait159523[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = NE
)]

hic_df_DE_bait159523 <- scores_TempClust_bait159523[, .(
  baitChr = paste0("chr", baitChr),  # Add 'chr' prefix if not already present
  baitStart,
  baitEnd,
  oeChr = paste0("chr", oeChr),
  oeStart,
  oeEnd,
  score = DE
)]


#### Selection of the contact of interest
bait_of_interest_159523 <- baitmap[
  chrom == "chr14" & 
    start == 69876518 & 
    end == 69877947
]


#### Build Hi-C list
hic_list <- list(
  hic_df_day14_bait159523,
  hic_df_NE_bait159523
)

names(hic_list) <- timepoints_labels

#### Parameters
params <- pgParams(
  chrom = "chr14", 
  chromstart = 69800000, chromend = 70150000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

params_bigwig <- pgParams(
  chrom = "14", 
  chromstart = 69800000, chromend = 70150000,
  assembly = "hg38",
  x = unit(0.5, "inches"), width = unit(2, "inches")
)

#### Color: arcsinh + global scale
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
# Day14        NE 
# 1.707487 1.665013   
global_max_asinh <- max(all_vals+0.01, na.rm = TRUE) # 1.717487

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


## Plot: PDF + PAGE
pdf(file.path(base_dir, "SuppFigure4D_top_CnT_PECHiC_Day14toNE_Constant_arcsinhColor.pdf"), width = 6, height = 10)

pageCreate(width = 6, height = 10, default.units = "inches", showGuides = FALSE)

## Gene track
genesPlot <- plotGenes(
  params=params, bg="#f6f6f6",
  y=0.5, width=3, height=0.5,
  just=c("left","top")
)

## Bait
plotRanges(
  data=bait_of_interest_159523, params=params,
  width=3, linecolor="red", fill="red",
  collapse=TRUE, y=0.7, height=0.2
)

## Layout constants
y_range_k27 <- c(0, 300)
y_range_k4 <- c(0, 300)
y_range_k27ac <- c(0, 150)

timepoint_spacing <- 2.8
k27me3_offset <- 0.3
k4me1_offset  <- 1.1
k27ac_offset  <- 1.9


## LOOP: ALL TIMEPOINTS
for(i in seq_along(timepoints_labels)) {
  
  tp_label <- timepoints_labels[i]
  tp_col   <- timepoints_columns[i]
  base_y <- 1.5 + (i-1)*timepoint_spacing
  
  ## label
  plotText(label=tp_label, fontsize=9, x=0.5, y=base_y)
  
  ## K27me3
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27me3_files[tp_col]),
    params=params_bigwig, fill="darkred", linecolor="darkred",
    y=base_y + k27me3_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27, scale=TRUE
  )
  plotText(label="H3K27me3 signal", fontsize=9,
           x=0.8, y=base_y+k27me3_offset-0.1, fontcolor="darkred")
  
  ## K4me1
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k4me1_files[tp_col]),
    params=params_bigwig, fill="blue", linecolor="blue",
    y=base_y + k4me1_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k4, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k4me1_offset-0.1, fontcolor="blue")
  
  ## K27ac
  plotSignal(
    data=paste0(CUTnTag_files_dir,
                k27ac_files[tp_col]),
    params=params_bigwig, fill="orange", linecolor="orange",
    y=base_y + k27ac_offset, width=3, height=0.5,
    just=c("left","top"), range=y_range_k27ac, scale=TRUE
  )
  plotText(label="H3K4me1 signal", fontsize=9,
           x=0.8, y=base_y+k27ac_offset-0.1, fontcolor="orange")
  
  ## ARCS
  df_arc <- hic_list[[tp_label]]
  df_arc <- df_arc[!is.na(score) & score > 0]
  
  colors = colorRampPalette(c("white","darkgreen"))(256)
  scorebins = seq(from=0, to=global_max_asinh, length.out=257)
  whichBin = findInterval(df_arc$score_asinh, scorebins)
  df_arc$color_asinh = colors[whichBin]
  
  if(nrow(df_arc) > 0) {
    
    y_arc <- base_y + k27ac_offset + 0.7
    
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
  scale="Mb", y=7.5
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
  y = unit(6 * sf, "inches"),
  rot = 90,
  gp = gpar(cex = 0.65)  ## smaller title
)

dev.off()

