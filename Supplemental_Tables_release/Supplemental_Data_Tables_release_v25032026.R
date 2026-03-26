################################################################################
## Supplemental Data and Tables release
## Dr Marina Nocente, Dr Mikhial Spivakov
################################################################################

library(data.table)
library(openxlsx)
library(dplyr)

################################################################################
### Peak Matrix with CHiCAGO scores for data releasing on OSF 
################################################################################

#######################
### Peakmatrix with CHiCAGO scores with a threshold of 5 for contact score in at least one time point of the transition + primded + NE and DE = DATA S1
#######################
fullMatrix <- fread("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
transitionMatrix <- fullMatrix[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist",
                                   "naive", "day1", "day3", "day5", "day7", "day10", "day14", "primed", "NE", "DE")] # 1322199      21

score_cols <- c("naive","day1","day3","day5","day7","day10","day14","primed","NE","DE")

# Keep only the rows of transitionMatrix where at least one column in score_cols has a value greater than 5.
transitionMatrix_filt <- transitionMatrix[
  transitionMatrix[, rowSums(.SD > 5) > 0, .SDcols = score_cols]
]

head(transitionMatrix_filt)
dim(transitionMatrix_filt) # 766,138     21


TSS = fread("~/Documents/Bioinformatics/GO_OtherEnds/2021.09_TSS_1toY_bedFormat.tsv")
setkey(TSS, V1, V2, V3)
transitionMatrix_filt_tss = foverlaps(transitionMatrix_filt, TSS, by.x=c("oeChr", "oeStart", "oeEnd"),
                    by.y=c("V1", "V2", "V3"), 
                    nomatch = NA)
# V4 is TSS name
transitionMatrix_filt_tss_unique = transitionMatrix_filt_tss[, .(baitChr=baitChr[1], baitStart=baitStart[1], 
                                                                 baitEnd=baitEnd[1], 
                                                                 baitName=baitName[1], oeChr=oeChr[1], 
                                                                 oeStart=oeStart[1], oeEnd=oeEnd[1], 
                                                                 oeName=oeName[1], dist=dist[1],
                                                                 naive=naive[1], day1=day1[1], day3=day3[1], 
                                                                 day5=day5[1], day7=day7[1], day10=day10[1], day14=day14[1], 
                                                                 primed=primed[1], NE=NE[1], DE=DE[1], 
                                                                 V4=paste(V4[!is.na(V4)], collapse=";")), 
                                                                 by=c("baitID", "oeID")]
head(transitionMatrix_filt_tss)
dim(transitionMatrix_filt_tss) # 768341     24

head(transitionMatrix_filt_tss_unique)
dim(transitionMatrix_filt_tss_unique) # 766138     22

setcolorder(transitionMatrix_filt_tss_unique, c(3:5,1,6:9,2,10:ncol(transitionMatrix_filt_tss_unique)))

# Adding PE IDs
transitionMatrix_filt_tss_unique[, oeNamePE:=oeName]
transitionMatrix_filt_tss_unique[, oeName:=apply(.SD, 1, function(x){
  if(x["V4"]!="" & x["oeNamePE"]!=".") {return(paste(x["V4"], x["oeNamePE"], sep=";"))}
  else{
    if(x["oeNamePE"]!="."){ return(x["oeNamePE"]) }
    if(x["V4"]!=""){ return(x["V4"]) }
  }
}), .SDcols=c("V4","oeNamePE")] # V4 is TSS name

transitionMatrix_filt_tss_unique[, oeNamePE:=NULL]
transitionMatrix_filt_tss_unique[, oeName1:= sapply(oeName, function(x)paste(unlist(x)[unlist(x)!=""], collapse=";"))]
transitionMatrix_filt_tss_unique[, oeName:=oeName1]
transitionMatrix_filt_tss_unique[, oeName1:=NULL]
transitionMatrix_filt_tss_unique[, V4:=NULL]
dim(transitionMatrix_filt_tss_unique) # 766138     21

fwrite(transitionMatrix_filt_tss_unique, "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Peakmatrix_chicago_score_cutoff5_NPtransition_differentiation.txt", sep = "\t", quote=F, col.names=T)


#######################
### Peakmatrix with CHiCAGO scores with a threshold of 5 for contact score in at least one PRC2 condition = DATA S2
#######################
fullMatrix <- fread("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/peakMatrix/PerChicagoScores/Day5_PRC2totalPooled_NPtransition_NE_DE.txt")
MatrixPRC2 <- fullMatrix[, c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist",
                             "Day5_PRC2_DMSO_total_pooled", "Day5_PRC2_PROTAC_total_pooled", "Day5_PRC2_UNC1999_total_pooled")] # 1322199      14

score_cols <- c("Day5_PRC2_DMSO_total_pooled","Day5_PRC2_PROTAC_total_pooled","Day5_PRC2_UNC1999_total_pooled")

# Keep only the rows of MatrixPRC2 where at least one column in score_cols has a value greater than 5.
MatrixPRC2_filt <- MatrixPRC2[
  MatrixPRC2[, rowSums(.SD > 5) > 0, .SDcols = score_cols]
]

head(MatrixPRC2_filt)
dim(MatrixPRC2_filt) # 434,524     14

TSS = fread("~/Documents/Bioinformatics/GO_OtherEnds/2021.09_TSS_1toY_bedFormat.tsv")
setkey(TSS, V1, V2, V3)
MatrixPRC2_filt_tss = foverlaps(MatrixPRC2_filt, TSS, by.x=c("oeChr", "oeStart", "oeEnd"),
                                      by.y=c("V1", "V2", "V3"), 
                                      nomatch = NA)
head(MatrixPRC2_filt_tss)
dim(MatrixPRC2_filt_tss) # 435,558     17

MatrixPRC2_filt_tss_unique = MatrixPRC2_filt_tss[, .(baitChr=baitChr[1], baitStart=baitStart[1], 
                                                                 baitEnd=baitEnd[1], 
                                                                 baitName=baitName[1], oeChr=oeChr[1], 
                                                                 oeStart=oeStart[1], oeEnd=oeEnd[1], 
                                                                 oeName=oeName[1], dist=dist[1],
                                                     Day5_PRC2_DMSO_total_pooled=Day5_PRC2_DMSO_total_pooled[1], 
                                                     Day5_PRC2_PROTAC_total_pooled=Day5_PRC2_PROTAC_total_pooled[1],
                                                     Day5_PRC2_UNC1999_total_pooled=Day5_PRC2_UNC1999_total_pooled[1],
                                                                 V4=paste(V4[!is.na(V4)], collapse=";")), 
                                                             by=c("baitID", "oeID")]


dim(MatrixPRC2_filt_tss_unique) # 434524     15

MatrixPRC2_filt_tss_unique[, oeNamePE:=oeName]
MatrixPRC2_filt_tss_unique[, oeName:=apply(.SD, 1, function(x){
  if(x["V4"]!="" & x["oeNamePE"]!=".") {return(paste(x["V4"], x["oeNamePE"], sep=";"))}
  else{
    if(x["oeNamePE"]!="."){ return(x["oeNamePE"]) }
    if(x["V4"]!=""){ return(x["V4"]) }
  }
}), .SDcols=c("V4","oeNamePE")] # V4 is TSS name

dim(MatrixPRC2_filt_tss_unique) # 434524     16

MatrixPRC2_filt_tss_unique[, oeNamePE:=NULL]
MatrixPRC2_filt_tss_unique[, oeName1:= sapply(oeName, function(x)paste(unlist(x)[unlist(x)!=""], collapse=";"))]
MatrixPRC2_filt_tss_unique[, oeName:=oeName1]
MatrixPRC2_filt_tss_unique[, oeName1:=NULL]
MatrixPRC2_filt_tss_unique[, V4:=NULL]

dim(MatrixPRC2_filt_tss_unique) # 434524     14

fwrite(MatrixPRC2_filt_tss_unique, "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Peakmatrix_chicago_score_cutoff5_PRC2perturbation.txt", sep = "\t", quote=F, col.names=T)


################################################################################
### Peak Matrix with raw CHiCAGO scores for PE-promoter and PE-PE contacts transition + temporal clusters (dynamic of contacts) - Table S1
################################################################################
Focused_pm <- fread("~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Focused_peakmatrix_PE_TSS_interaction.tsv")
head(Focused_pm)
dim(Focused_pm) #  12712    12


Focused_pm_coor <- merge(Focused_pm, transitionMatrix_filt_tss_unique, by=c("baitID", "oeID"))
Focused_pm_coor <- Focused_pm_coor[,c("baitChr", "baitStart", "baitEnd", "baitID", "baitName", "oeChr", "oeStart", "oeEnd", "oeID", "oeName", "dist", "naive.x", "day1.x", "day3.x", "day5.x", "day7.x", "day10.x", "day14.x", "primed.x", "NE.x", "DE.x")]

colnames(Focused_pm_coor) <- c(
  "baitChr", "baitStart", "baitEnd", "baitID", "baitName", 
  "oeChr", "oeStart", "oeEnd", "oeID", "oeName", 
  "dist", "naive", "day1", "day3", "day5", "day7", "day10", "day14", 
  "primed", "NE", "DE"
)

dim(Focused_pm_coor) # 12712    21
head(Focused_pm_coor) 

cluster = fread("~/Documents/Bioinformatics/Final_scripts_paper/PoisedEnhancers/Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt")

Focused_pm_coor_clusters = merge(Focused_pm_coor, cluster[, c("baitID", "oeID", "merged_dtw_trans_cl6_2with6")], by=c("baitID", "oeID"), all.x=TRUE)

dim(Focused_pm_coor_clusters) # 12712    22
head(Focused_pm_coor_clusters) 

write.xlsx(Focused_pm_coor_clusters, file = "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Table_S1.xlsx", sheetName = "Supplemental_Table_S1", overwrite = TRUE)



################################################################################
### The PE clusters - Table S2
################################################################################

PEclusters <- fread("~/Documents/Bioinformatics/PECHiC_PRC2_DPPA/Chromatin_baits_clustering_temporal_clustering_test2_26-02-2025/clust_5_full_H3K27me3_spreading.txt")
head(PEclusters)
setnames(PEclusters, "V1", "baitID")
setnames(
  PEclusters,
  old = c("Dppa2_naive", "Dppa2_primed", "Dppa4_naive", "Dppa4_primed"),
  new = c("DPPA2_naive", "DPPA2_primed", "DPPA4_naive", "DPPA4_primed")
)

dim(PEclusters) # 32063    24
head(PEclusters)


write.xlsx(PEclusters, file = "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Table_S2.xlsx", overwrite = TRUE)



################################################################################
### PRC2 analysis - Table S3
################################################################################

pechic_dese2_stats <- fread("~/Documents/Project_writing/Paper_PEs/Data_release/PECHiC_deseq2_stats_PRC2-perturb_Tempclusters.txt")
head(pechic_dese2_stats)
dim(pechic_dese2_stats) # 9998    5
pechic_dese2_stats[,  merged_dtw_trans_cl6_2with6:=NULL]
pechic_dese2_stats[,  oeID:=NULL]
pechic_dese2_stats = unique(pechic_dese2_stats, by="baitID")
dim(pechic_dese2_stats) # 6445    3

# Rename columns
colnames(pechic_dese2_stats)[colnames(pechic_dese2_stats) == "stat_res_UNC1999_vs_DMSO_merged"] <- "tval_PECHiC_UNC1999_vs_DMSO"
colnames(pechic_dese2_stats)[colnames(pechic_dese2_stats) == "stat_res_PROTAC_vs_DMSO_merged"] <- "tval_PECHiC_PROTAC_vs_DMSO"

baitmap <- fread("/Users/mnocente/Documents/Project_writing/Paper_PEs/Data_release/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")
colnames(baitmap) <- c("chr", "start", "end", "baitID", "baitName")
dim(baitmap) # 32142     5

pechic_dese2_stats <- merge(pechic_dese2_stats, baitmap, by="baitID") 
dim(pechic_dese2_stats) # 6445    7

# Reorder columns
pechic_dese2_stats <- pechic_dese2_stats[, c("chr", "start", "end", "baitID", "baitName", "tval_PECHiC_PROTAC_vs_DMSO", "tval_PECHiC_UNC1999_vs_DMSO")]
head(pechic_dese2_stats)

# CUT&Tag data
res_PROTAC_DMSO_CnT_all_dmsoBin_normal = fread("~/Documents/Project_writing/Paper_PEs/Data_release/CUT&Tag_PRC2/final/CUTnTag_deseq2_res_PROTAC_DMSO_dmsoMean_RunSeq2.txt")
res_PROTAC_DMSO_CnT_all_dmsoBin_normal[, c("chr", "start", "end") := tstrsplit(region, "_")]
res_PROTAC_DMSO_CnT_all_dmsoBin_normal[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
res_PROTAC_DMSO_CnT_all_dmsoBin_normal <- res_PROTAC_DMSO_CnT_all_dmsoBin_normal[,-1] # 32142    11

res_UNC1999_DMSO_CnT_all_dmsoBin_normal = fread("~/Documents/Project_writing/Paper_PEs/Data_release/CUT&Tag_PRC2/final/CUTnTag_deseq2_res_UNC1999_DMSO_dmsoMean_RunSeq2.txt")
res_UNC1999_DMSO_CnT_all_dmsoBin_normal[, c("chr", "start", "end") := tstrsplit(region, "_")]
res_UNC1999_DMSO_CnT_all_dmsoBin_normal[, c("start", "end") := lapply(.SD, as.integer), .SDcols = c("start", "end")]
res_UNC1999_DMSO_CnT_all_dmsoBin_normal <- res_UNC1999_DMSO_CnT_all_dmsoBin_normal[,-1] # 32142    11

## For UNC and PROTAC all baits with a temporal class
pechic_dese2_stats_H3K27me3_CmTstats = merge(pechic_dese2_stats, res_PROTAC_DMSO_CnT_all_dmsoBin_normal[,c("chr", "start", "end", "stat", "DMSO_mean", "DMSO_mean_bin")], by=c("chr", "start", "end")) 
dim(pechic_dese2_stats_H3K27me3_CmTstats) # 6445   10
head(pechic_dese2_stats_H3K27me3_CmTstats)

# Remove a small number of PEs with no detectable H3K27me3 at day 5
pechic_dese2_stats_H3K27me3_CmTstats = pechic_dese2_stats_H3K27me3_CmTstats[DMSO_mean>0]
dim(pechic_dese2_stats_H3K27me3_CmTstats) # 6204   10

colnames(pechic_dese2_stats_H3K27me3_CmTstats)[colnames(pechic_dese2_stats_H3K27me3_CmTstats) == "stat"] <- "tval_H3K27me3_PROTAC_vs_DMSO"
colnames(pechic_dese2_stats_H3K27me3_CmTstats)[colnames(pechic_dese2_stats_H3K27me3_CmTstats) == "DMSO_mean"] <- "mean_DMSO_H3K27me3"
colnames(pechic_dese2_stats_H3K27me3_CmTstats)[colnames(pechic_dese2_stats_H3K27me3_CmTstats) == "DMSO_mean_bin"] <- "DMSO_H3K27me3_bin"

pechic_dese2_stats_H3K27me3_CmTstats = merge(pechic_dese2_stats_H3K27me3_CmTstats, res_UNC1999_DMSO_CnT_all_dmsoBin_normal[,c("chr", "start", "end", "stat")], by=c("chr", "start", "end")) 
dim(pechic_dese2_stats_H3K27me3_CmTstats) # 6204 11
head(pechic_dese2_stats_H3K27me3_CmTstats)

colnames(pechic_dese2_stats_H3K27me3_CmTstats)[colnames(pechic_dese2_stats_H3K27me3_CmTstats) == "stat"] <- "tval_H3K27me3_UNC1999_vs_DMSO"

# Reorder
pechic_dese2_stats_H3K27me3_CmTstats <- pechic_dese2_stats_H3K27me3_CmTstats[, c("chr", "start", "end", "baitID", "baitName", "tval_PECHiC_PROTAC_vs_DMSO", "tval_H3K27me3_PROTAC_vs_DMSO", "tval_PECHiC_UNC1999_vs_DMSO", "tval_H3K27me3_UNC1999_vs_DMSO", "mean_DMSO_H3K27me3", "DMSO_H3K27me3_bin")]
head(pechic_dese2_stats_H3K27me3_CmTstats)
dim(pechic_dese2_stats_H3K27me3_CmTstats) # 6204   11


write.xlsx(pechic_dese2_stats_H3K27me3_CmTstats, file = "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Table_S3.xlsx", overwrite = TRUE)



################################################################################
### NE_DE_analysis DESeq analysis results - Table S4
################################################################################

H3K27me3_NE_DE <- fread("~/Documents/Bioinformatics/NE_DE_analysis/Supplementary_Table_H3K27me3.tsv")
H3K27ac_NE_DE <- fread("~/Documents/Bioinformatics/NE_DE_analysis/Supplementary_Table_H3K27ac.tsv")
H3K4me1_NE_DE <- fread("~/Documents/Bioinformatics/NE_DE_analysis/Supplementary_Table_H3K4me1.tsv")

setnames(H3K27me3_NE_DE, "day14", "day14_arcsinh_score")
setnames(H3K27me3_NE_DE, "NE" , "NE_arcsinh_score")
setnames(H3K27me3_NE_DE, "DE" , "DE_arcsinh_score")
setnames(H3K27me3_NE_DE, "merged_dtw_trans_cl6_2with6" , "PE_contact_clusters")
setnames(H3K27me3_NE_DE, "log2FC_DE_day14" , "H3K27me3_log2FC_DE_day14")
setnames(H3K27me3_NE_DE, "pvalue_DE_day14" , "H3K27me3_pvalue_DE_day14")
setnames(H3K27me3_NE_DE, "padj_DE_day14" , "H3K27me3_padj_DE_day14")
setnames(H3K27me3_NE_DE, "log2FC_NE_day14" , "H3K27me3_log2FC_NE_day14")
setnames(H3K27me3_NE_DE, "pvalue_NE_day14" , "H3K27me3_pvalue_NE_day14")
setnames(H3K27me3_NE_DE, "padj_NE_day14" , "H3K27me3_padj_NE_day14")

H3K27me3_NE_DE <- H3K27me3_NE_DE[,c("baitID", "oeID", "day14_arcsinh_score", "DE_arcsinh_score", "NE_arcsinh_score", "PE_contact_clusters",
                                    "H3K27me3_log2FC_DE_day14", "H3K27me3_pvalue_DE_day14", "H3K27me3_padj_DE_day14", "H3K27me3_log2FC_NE_day14", "H3K27me3_pvalue_NE_day14", "H3K27me3_padj_NE_day14")]

dim(H3K27me3_NE_DE) # 12712    12

setnames(H3K27ac_NE_DE, "day14", "day14_arcsinh_score")
setnames(H3K27ac_NE_DE, "NE" , "NE_arcsinh_score")
setnames(H3K27ac_NE_DE, "DE" , "DE_arcsinh_score")
setnames(H3K27ac_NE_DE, "merged_dtw_trans_cl6_2with6" , "PE_contact_clusters")
setnames(H3K27ac_NE_DE, "log2FC_DE_day14" , "H3K27ac_log2FC_DE_day14")
setnames(H3K27ac_NE_DE, "pvalue_DE_day14" , "H3K27ac_pvalue_DE_day14")
setnames(H3K27ac_NE_DE, "padj_DE_day14" , "H3K27ac_padj_DE_day14")
setnames(H3K27ac_NE_DE, "log2FC_NE_day14" , "H3K27ac_log2FC_NE_day14")
setnames(H3K27ac_NE_DE, "pvalue_NE_day14" , "H3K27ac_pvalue_NE_day14")
setnames(H3K27ac_NE_DE, "padj_NE_day14" , "H3K27ac_padj_NE_day14")


H3K27ac_NE_DE <- H3K27ac_NE_DE[,c("baitID", "oeID", "day14_arcsinh_score", "DE_arcsinh_score", "NE_arcsinh_score", "PE_contact_clusters",
                                  "H3K27ac_log2FC_DE_day14", "H3K27ac_pvalue_DE_day14", "H3K27ac_padj_DE_day14", "H3K27ac_log2FC_NE_day14", "H3K27ac_pvalue_NE_day14", "H3K27ac_padj_NE_day14")]

dim(H3K27ac_NE_DE) # 12712    12

setnames(H3K4me1_NE_DE, "day14", "day14_arcsinh_score")
setnames(H3K4me1_NE_DE, "NE" , "NE_arcsinh_score")
setnames(H3K4me1_NE_DE, "DE" , "DE_arcsinh_score")
setnames(H3K4me1_NE_DE, "merged_dtw_trans_cl6_2with6" , "PE_contact_clusters")
setnames(H3K4me1_NE_DE, "log2FC_DE_day14" , "H3K4me1_log2FC_DE_day14")
setnames(H3K4me1_NE_DE, "pvalue_DE_day14" , "H3K4me1_pvalue_DE_day14")
setnames(H3K4me1_NE_DE, "padj_DE_day14" , "H3K4me1_padj_DE_day14")
setnames(H3K4me1_NE_DE, "log2FC_NE_day14" , "H3K4me1_log2FC_NE_day14")
setnames(H3K4me1_NE_DE, "pvalue_NE_day14" , "H3K4me1_pvalue_NE_day14")
setnames(H3K4me1_NE_DE, "padj_NE_day14" , "H3K4me1_padj_NE_day14")

H3K4me1_NE_DE <- H3K4me1_NE_DE[,c("baitID", "oeID", "day14_arcsinh_score", "DE_arcsinh_score", "NE_arcsinh_score", "PE_contact_clusters",
                                  "H3K4me1_log2FC_DE_day14", "H3K4me1_pvalue_DE_day14", "H3K4me1_padj_DE_day14", "H3K4me1_log2FC_NE_day14", "H3K4me1_pvalue_NE_day14", "H3K4me1_padj_NE_day14")]

dim(H3K4me1_NE_DE) # 12712    12


dyn_histModif <- merge(H3K27me3_NE_DE, H3K4me1_NE_DE, by = c("baitID", "oeID", "day14_arcsinh_score", "DE_arcsinh_score", "NE_arcsinh_score", "PE_contact_clusters"))
dyn_histModif <- merge(dyn_histModif, H3K27ac_NE_DE, by = c("baitID", "oeID", "day14_arcsinh_score", "DE_arcsinh_score", "NE_arcsinh_score", "PE_contact_clusters"))
dim(dyn_histModif) # 12712    24
head(dyn_histModif)

write.xlsx(dyn_histModif[,-grep("pvalue", colnames(dyn_histModif)), with=F], file = "~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Table_S4.xlsx", overwrite = TRUE)


