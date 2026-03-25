######################################################################################################
### Focus on PE-TSS and PE-PE contacts
### Author: Dr. Monica Della Rosa
###
### This script generates the file "Focused_peakmatrix_PE_TSS_interaction.tsv"
###
######################################################################################################

#######################################
## Libraries
#######################################
library(data.table)
library(tidyverse)
library(GenomicRanges)

#######################################
## File paths
#######################################
peakMatrix_file <- file.path("~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/Peakmatrix_chicago_score_cutoff5_NPtransition_differentiation.txt")
rmap_file <- file.path("~/Documents/Project_writing/Paper_PEs/Data_release/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.rmap")
baitmap_file <- file.path("~/Documents/Bioinformatics/GO_OtherEnds/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt")
TSS_file <- file.path("~/Documents/Bioinformatics/GO_OtherEnds/2021.09_TSS_1toY_bedFormat.tsv")

#######################################
## Prepare the file
#######################################
#repeat PCA on peakMatrix of combined, non-downsampled data
peakMatrix <- fread(peakMatrix_file)
setcolorder(peakMatrix, c("baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist", "naive", "day1", "day3", "day5", "day7", "day10", "day14", "primed", "DE", "NE"))
pm.cis <- peakMatrix %>% filter(is.na(dist) == FALSE)
score <- pm.cis[,c(4,9,12:21)]
asinhScore <- apply(score[,3:12], 2, asinh)
asinhScore <- cbind(score[,1:2], asinhScore)
setDT(asinhScore)

rmap <- fread(rmap_file)
rmap.gr <- makeGRangesFromDataFrame(
  rmap[ ,list(seqname=paste(V1), start=paste(V2), end=paste(V3))]
)
baitmap <- fread(baitmap_file)
TSS <- fread(TSS_file)
(TSS.gr <- makeGRangesFromDataFrame(
  TSS[ ,list(seqname=paste(V1), start=paste(V2), end=paste(V3), geneName=paste(V4))], keep.extra.columns = T
)
)
b2b <- pm.cis[pm.cis$oeName!=".",]
#retrieve fragment ID for TSS
OETss <- findOverlaps(TSS.gr, rmap.gr, minoverlap = 10, select = "all")
rfIDtss <- bind_cols(TSS[OETss@from,], rmap[OETss@to, ])

setkey(asinhScore, oeID)
asinh.b2bTss <- asinhScore[oeID %in% c(b2b$oeID, rfIDtss$V4...8)]
dim(asinh.b2bTss)
#This is were only b2b interactions and interactions with annotated TSSs are included. 
#Note that when computing the overlap between restriction fragments and TSSs, in order to retrieve the fragment IDs, 
#I've requested a minimum overlap of 10bp. The nrow(asinh.b2bTss) == 12,712. 

fwrite(asinh.b2bTss, "~/Documents/Project_writing/Paper_PEs/Data_release/Focused_peakmatrix_PE_TSS_interaction.tsv", sep = "\t")
