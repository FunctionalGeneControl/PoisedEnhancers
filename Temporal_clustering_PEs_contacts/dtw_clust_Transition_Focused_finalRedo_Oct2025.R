library(data.table)
library(dtwclust)
setDTthreads(4)

## Load the PECHi-C contact matrix from the Monica'selection
# Each row corresponds to a genomic contact with columns representing signal strength (asinh-transformed) across 
# the naive-to-primed transition time points, H9 primed cells and NE and DE cells.
setwd("~/Documents/Bioinformatics/Final_scripts_paper/OSF_data/")
pepm = fread("Focused_peakmatrix_PE_TSS_interaction.tsv")
dim(pepm) # 12712    12

## Filter contacts based on contact score thresholds and time point
# trying dtwclust without Primed and those contacts that only show up in Primed, NE or DE...
# using the chicago cutoff of ~3 (asinh cutoff of 1.8) for the other time points
pepm_transition = pepm[naive>1.8|day1>1.8|day3>1.8|day5>1.8|day7>1.8|day10>1.8|day14>1.8] 
dim(pepm_transition) # retained 9998 rows, 78.7% of total

## Run Dynamic Time Warping (DTW) clustering
# Perform partitional clustering on temporal contact patterns using DTW as the distance metric. 
# Test k = 3 to 8 clusters.
hc_dtw_transition_MonicaSelect <- tsclust(pepm_transition[, 3:9], # select columns corresponding to time series
                        type = "partitional", 
                        k = c(3,4,5,6,7,8), distance = "dtw_basic",
                        trace = TRUE, seed=123)


## Add cluster assignments to the data table
pepm_transition[, dtw_trans_cl3:=hc_dtw_transition_MonicaSelect[[1]]@cluster]
pepm_transition[, dtw_trans_cl4:=hc_dtw_transition_MonicaSelect[[2]]@cluster]
pepm_transition[, dtw_trans_cl5:=hc_dtw_transition_MonicaSelect[[3]]@cluster]
pepm_transition[, dtw_trans_cl6:=hc_dtw_transition_MonicaSelect[[4]]@cluster]
pepm_transition[, dtw_trans_cl7:=hc_dtw_transition_MonicaSelect[[5]]@cluster]
pepm_transition[, dtw_trans_cl8:=hc_dtw_transition_MonicaSelect[[6]]@cluster]


## Save output
fwrite(pepm_transition, "~/Documents/Bioinformatics/Final_scripts_paper/Temporal_clustering_PEs_contacts/pepm_transition_dtwCluster_Focused_Oct2025.txt", sep="\t")
