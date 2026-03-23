################################################################################
### Association between PE chromatin features and their contact dynamics during the naive-to-primed transition - Figure 4 and Supplemental Figure S2 D-F
### Author: Dr. Marina Nocente
###
### This script reproduces Figure 4 panels:
###  Heatmap showing the extent of H3K27me3 spreading around PEs in primed cells, as well as the enrichment signals of transcription factors and cofactors at all analysed PEs in naive and primed hPSCs. 
###  The extent of H3K27me3 spreading around the PEs assigned to different chromatin clusters in primed cells
###  Normalised signal intensities of DPPA4, CTCF, and OCT4 across PE chromatin clusters in naive and primed cells
###  Heatmap showing the enrichment (or depletion) of PE chromatin clusters for each class of temporal contact dynamics
###
### This script reproduces some Supplemental Figure S2 panels:
###  Normalised signal intensities of SOX2, NANOG, DPPA2, and Mediator across PE chromatin clusters in naive and primed cells
###  Proportion of PEs containing CpG islands in the different classes of poised enhancers
###  Spreading of H3K27me3 around PEs involved in different temporal contact dynamics
################################################################################

########################################
### Load required libraries
########################################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ppclust)
library(pheatmap)

set.seed(123)

########################################
### User configuration
########################################
base_dir <- "~/Documents/Bioinformatics"
working_dir <- file.path(base_dir,"PECHiC_PRC2_DPPA/Chromatin_baits_clustering_temporal_clustering_test2_26-02-2025")

## Input files
H3K27me3_matrix_primed  <- file.path(base_dir, "CnT_transition/K27me3_heatmap_matrix.txt")
sorted_regions_primed <- file.path(base_dir, "CnT_transition/K27me3_heatmap_matrix_sorted_regions.txt")

baitmap_file <- file.path(base_dir, "Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")

peakMatrix <- file.path(base_dir, "Final_scripts_paper/OSF_data/Focused_peakmatrix_PE_TSS_interaction.tsv")
countsNorm <- file.path(working_dir, "Counts_DESeq2norm_fragLenNorm_ALL.txt")
lassoPred <- file.path(working_dir, "2209_LassoPredictorsMatrix_FragLenNorm_Scaled.tsv")

Temporalcluster_file <- file.path(base_dir, "Final_scripts_paper/PoisedEnhancers/Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt")


################################################################################
### Figure 4A
################################################################################

########################################
### Load H3K27me3 primed CUT&Tag matrix
########################################
# Columns 1–6: region metadata
# Columns 7+: per-bin signal
# Each row of mat_primed is a region (PE).
mat_primed <- read.table(H3K27me3_matrix_primed, header = FALSE, comment.char = "@", sep = "\t")
dim(mat_primed) # 32075  2006

all_vals <- unlist(mat_primed[, 7:ncol(mat_primed)]) 

############################
## Compute H3K27me3 spreading width (excluding PE center)
############################
threshold <- 25   # adjust based on histogram hist(all_vals)
bin_size  <- 10  # because I used --binSize 10 in computeMatrix

# Identify the central PE bin
num_bins <- ncol(mat_primed) - 6   # number of signal columns
center_bin <- 6 + ceiling(num_bins / 2)  # column index of PE center

#Compute width_kb without the center bin
# Count number of bins above threshold, excluding PE bin
signal_cols_noPE <- c(7:(center_bin-1), (center_bin+1):ncol(mat_primed))
mat_primed$width_kb_noPE <- apply(mat_primed[, signal_cols_noPE], 1, function(x) sum(x > threshold) * 10 / 1000) 

############################
## Load H3K27me3 sorted regions (primed)
############################
### Read sorted regions
sorted_regions_primed <- read.table(sorted_regions_primed, header = FALSE, skip = 1,   # skip metadata/header line
                             sep = "\t", stringsAsFactors = FALSE)

# Assign proper column names
colnames(sorted_regions_primed) <- c("chrom","start","end","name","score","strand",
                              "thickStart","thickEnd","itemRGB","blockCount",
                              "blockSizes","blockStart","deepTools_group")
head(sorted_regions_primed) # 32075    13

### Match the width to the sorted order
sorted_regions_primed$width_kb_noPE <- mat_primed$width_kb_noPE
head(sorted_regions_primed)

# Keep only the columns of interest and remove the Mitochondria (which are NA)
sorted_regions_primed_short <- sorted_regions_primed[, c("chrom", "start", "end", "width_kb_noPE")]
sorted_regions_primed_short <- sorted_regions_primed_short %>%
  filter(chrom != "MT")
head(sorted_regions_primed_short)
dim(sorted_regions_primed_short) # 32063     4


############################
## Add baitID
############################
baitmap <- fread(baitmap_file)
setnames(baitmap, old = colnames(baitmap),
         new = c("chrom", "start", "end", "baitID", "PEname"))

sorted_regions_short_baitID_primed <- sorted_regions_primed_short %>%
  left_join(baitmap,
            by = c("chrom" = "chrom",
                   "start" = "start",
                   "end"   = "end"))
head(sorted_regions_short_baitID_primed) # 32063     6

fwrite(sorted_regions_short_baitID_primed, file.path(working_dir, "sorted_regions_short_baitID_primed.txt"))

############################
## Load clustering input matrices
############################
pepm <- fread(peakMatrix)
countsNorm <- fread(countsNorm) # contains enrichment signals of transcription factors (DPPA2, DPPA4, CTCF, OCT4, NANOG, SOX2, TFAP2C) and cofactors (mediator, BRM, BRG1, BAF155) 
lassoPred <- fread(lassoPred)

############################
## Merge spreading width with normalised counts and prepare the clustering matrix
############################
countsNorm <- merge(countsNorm,
                    sorted_regions_short_baitID_primed[, c("baitID", "width_kb_noPE")],
                    by.x = "5KbID",
                    by.y = "baitID")

counts_norm_mat = as.matrix(countsNorm[, c(6:ncol(countsNorm)), with=F])
rownames(counts_norm_mat) = countsNorm$`5KbID`


## Select the useful columns, order and rename
counts_norm_mat_Wkb_noPE = counts_norm_mat[,c("width_kb_noPE", "Naive_Dppa2_sum" , 
                                              "Primed_Dppa2_sum","Naive_Dppa4_sum" , "Primed_Dppa4_sum",
                                              "CTCF_naive_sum", "CTCF_primed_sum","Mediator_naive_sum", "Mediator_primed_sum", "POUF51_naive_sum","POU5F1_primed_sum",
                                              "NANOG_naive_sum","NANOG_primed_sum",  "BRM_naive_sum", "BRG1_naive_sum", "BRG1_primed_sum",
                                              "BAF155_naive_sum", "BAF155_primed_sum", "TFAP2C_naive_sum","TFAP2C_primed_sum", "SOX2_naive_sum","Sox2_primed_sum")]

colnames(counts_norm_mat_Wkb_noPE) <- c("H3K27me3_spreading", "Dppa2_naive" , 
                                        "Dppa2_primed","Dppa4_naive" , "Dppa4_primed", "CTCF_naive", "CTCF_primed","Mediator_naive", "Mediator_primed", "OCT4_naive","OCT4_primed",
                                        "NANOG_naive","NANOG_primed",  "BRM_naive", "BRG1_naive", "BRG1_primed",
                                        "BAF155_naive", "BAF155_primed", "TFAP2C_naive","TFAP2C_primed", "SOX2_naive","Sox2_primed")


############################
## The clustering function
############################
#### This function performs k-means clustering on a numeric matrix (asinh-transformed normalised counts) and produces:
#### a row-level heatmap (all elements grouped by cluster)
#### a centroid heatmap (mean profile per cluster)
#### a returned matrix containing cluster assignments with descriptive labels

clusterCounts <- function(mat,
                          nclust,
                          excludeCols = NULL,
                          useAsinh = TRUE,
                          plotPerRow = TRUE,
                          plotCentroid = TRUE,
                          normToMax = TRUE,
                          maxQuantile = 0.999,
                          perRowCols = colorRampPalette(
                            c("white","lightblue","skyblue","steelblue","darkblue","darkblue")
                          )(256),
                          perRowBreaks = seq(0, 1, length.out = 257),
                          centroidCols = colorRampPalette(
                            c("white","lightblue","skyblue","steelblue","darkblue","darkblue")
                          )(256),
                          centroidBreaks = seq(-1, 1, length.out = 257),
                          # Descriptive cluster labels
                          cluster_labels = c(
                            "1" = "CTCF + TFs",
                            "2" = "Broad H3K27me3",
                            "3" = "High Occup",
                            "4" = "Intermediate",
                            "5" = "Low Occup"
                          ),
                          # Y-axis / legend order
                          cluster_order = c("2","3","1","4","5")
) {
  set.seed(123)
  
  # Exclude columns if needed
  if (length(excludeCols)) {
    mat <- mat[, -which(colnames(mat) %in% excludeCols)]
  }
  
  # Optional asinh transform
  if (useAsinh) {
    mat <- asinh(mat)
  }
  
  # K-means clustering
  ekm_counts <- ekm(scale(mat), centers = nclust)
  
  # Compute and print cluster sizes with descriptive names
  cluster_sizes <- table(cluster_labels[as.character(ekm_counts$cluster)])
  print(cluster_sizes)
  
  # Combine matrix with cluster assignments
  mat_clust <- cbind(mat, ekm_counts$cluster)
  colnames(mat_clust)[ncol(mat_clust)] <- "cluster"
  
  # Reorder rows according to cluster_order
  mat_clust <- mat_clust[order(factor(mat_clust[, "cluster"], levels = cluster_order)), ]
  
  # Optional normalization
  if (normToMax) {
    for (i in 1:(ncol(mat_clust) - 1)) {
      mat_clust[, i] <- mat_clust[, i] / quantile(mat_clust[, i], maxQuantile)
    }
  }
  
  # Per-row heatmap
  if (plotPerRow) {
    annotation_row <- data.frame(
      Cluster = factor(
        cluster_labels[as.character(mat_clust[, "cluster"])],
        levels = cluster_labels[cluster_order]
      )
    )
    rownames(annotation_row) <- rownames(mat_clust)
    
    annotation_colors <- list(
      Cluster = setNames(
        rainbow(length(cluster_order)),
        cluster_labels[cluster_order]
      )
    )
    
    pheatmap::pheatmap(
      mat_clust[, -ncol(mat_clust)],
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      show_rownames = FALSE,
      show_colnames = TRUE,
      color = perRowCols,
      breaks = perRowBreaks,
      raster = TRUE
    )
  }
  
  # Centroid heatmap
  if (plotCentroid) {
    clustMeans <- as.data.table(mat_clust)[,
                                           lapply(.SD, mean),
                                           by = "cluster",
                                           .SDcols = 1:(ncol(mat_clust) - 1)
    ]
    print(names(clustMeans))
    pheatmap::pheatmap(
      clustMeans[, 2:(ncol(clustMeans)), with = FALSE],
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = centroidCols,
      breaks = centroidBreaks
    )
  }
  
  # Replace numeric cluster IDs with descriptive names
  mat_clust[, "cluster"] <- cluster_labels[as.character(mat_clust[, "cluster"])]
  colnames(mat_clust)[ncol(mat_clust)] <- "Cluster_Name"
  
  # Attach cluster sizes as an attribute for backward compatibility
  attr(mat_clust, "cluster_sizes") <- cluster_sizes
  
  # Return matrix/data.frame directly
  invisible(mat_clust)
}


############################
## Run clustering (5 clusters) and plot
############################
pdf(file.path(working_dir, "clust_5_full_H3K27me3_spreading_Fig4A.pdf"))
clust_5_full_Wkb_noPE = clusterCounts(counts_norm_mat_Wkb_noPE, 5, plotCentroid = FALSE)
dev.off()

clust_5_full_Wkb_noPE <- as.data.frame(clust_5_full_Wkb_noPE)
head(clust_5_full_Wkb_noPE)
write.table(clust_5_full_Wkb_noPE, file.path(working_dir, "clust_5_full_H3K27me3_spreading.txt"), sep = "\t", quote=F, col.names=T, row.names=T)


################################################################################
### Figure 4B (top part) - The extent of H3K27me3 spreading around the PEs assigned to different chromatin clusters in primed cells
################################################################################

# Matrix of PEs cluster classes (normalised signals and PE cluster for each baitID)
clust_5_full_Wkb_noPE
head(clust_5_full_Wkb_noPE)
clust_5_full_Wkb_noPE$baitID <- rownames(clust_5_full_Wkb_noPE)

# In previous steps, baitID was converted to character.
# Since sorted_regions_short_baitID_primed$baitID is integer, we convert this column back to integer to avoid join errors.
clust_5_full_Wkb_noPE <- clust_5_full_Wkb_noPE %>%mutate(baitID = as.integer(baitID))
fwrite(clust_5_full_Wkb_noPE, file.path(working_dir, "clust_5_full_H3K27me3_spreading_baitID.txt"), sep = "\t", quote=F, col.names=T, row.names=T)


# Merge by baitID to add cluster_name to sorted_regions_short_baitID_primed
sorted_regions_with_PEcluster_primed <- sorted_regions_short_baitID_primed %>%
  left_join(clust_5_full_Wkb_noPE %>% select(baitID, Cluster_Name),
            by = "baitID")

head(sorted_regions_with_PEcluster_primed) # 32063     7

# Define biologically meaningful cluster order
sorted_regions_with_PEcluster_primed$Cluster_Name <- factor(
  sorted_regions_with_PEcluster_primed$Cluster_Name,
  levels = c("Broad H3K27me3", "High Occup", "CTCF + TFs", "Intermediate", "Low Occup")  
)

# Generate boxplot of H3K27me3 spreading per cluster
pdf(file.path(working_dir, "Figure_4Btop_Boxplots_H3K27me3_spreading_primed_VS_PEs_classes.pdf"))
ggplot(sorted_regions_with_PEcluster_primed, 
       aes(x = Cluster_Name, y = width_kb_noPE)) +
  geom_boxplot(outlier.shape = NA) +   # hide outliers
  theme_bw() +
  labs(
    title = "H3K27me3 spreading",
    x = "Cluster",
    y = "Spreading around PEs (kb)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



################################################################################
### Figure 4B (bottom part) and Supplemental Figure S2D - Boxplots of TFs and cofactors normalised signal intensity in naive and primed cells
################################################################################
# Ensure Cluster_Name is a factor with desired order
clust_5_full_Wkb_noPE$Cluster_Name <- factor(
  clust_5_full_Wkb_noPE$Cluster_Name,
  levels = c("Broad H3K27me3", "High Occup", "CTCF + TFs", 
             "Intermediate", "Low Occup")
)

# TF markers to plot
markers <- c("CTCF", "Dppa4", "Dppa2", "OCT4", "Mediator", "NANOG", "SOX2")

# Special column names for markers that differ
col_fixes <- list(
  "SOX2" = c("SOX2_naive", "Sox2_primed")
)

# Loop over markers
for (marker in markers) {
  
  # Determine Naive/Primed column names
  if (marker %in% names(col_fixes)) {
    naive_col <- col_fixes[[marker]][1]
    primed_col <- col_fixes[[marker]][2]
  } else {
    naive_col <- paste0(marker, "_naive")
    primed_col <- paste0(marker, "_primed")
  }
  
  # Reshape to long format
  df_long <- clust_5_full_Wkb_noPE %>%
    select(Cluster_Name, all_of(c(naive_col, primed_col))) %>%
    pivot_longer(
      cols = all_of(c(naive_col, primed_col)),
      names_to = "Condition",
      values_to = "Signal"
    ) %>%
    mutate(
      Signal = as.numeric(Signal)    # convert to numeric
    ) %>%
    filter(!is.na(Signal))           # remove NA values
  
  # Factor for Condition to control order
  df_long$Condition <- factor(
    df_long$Condition,
    levels = c(naive_col, primed_col),
    labels = c("Naive", "Primed")
  )
  
  # Plot boxplot
  pdf(file.path(working_dir, paste0(marker, "_cluster_boxplot2_Fig4B_SuppFigS2D.pdf")), width = 8, height = 6)
  print(
    ggplot(df_long, aes(x = Cluster_Name, y = Signal, fill = Condition)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = c("Naive" = "forestgreen", "Primed" = "darkgoldenrod1")) +
      labs(x = "Cluster", y = paste("Normalized", marker, "signal"), fill = "Condition") +
      ylim(0,1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
}



################################################################################
### Supplemental Figure S2E - Proportion of PEs containing CpG islands in the different classes of poised enhancers.
################################################################################

############################
## Prepare CpG annotation
############################
# In lassoPred, the baitID are the baitID from pepm in the same order
# Add the baitID from pepm in a new column in lassoPred
# lassoPred contains multiple rows per bait (contacts).
lassoPred$baitID <- pepm$baitID
dim(lassoPred) # 12712    

# Keep one row per baitID for CGI annotation.
# Subset clustered PEs to those present in lassoPred
clust_5_CGsubsetted_noPE = clust_5_full_Wkb_noPE %>% filter(baitID %in% lassoPred$baitID)

## Merge the CGI column of lassoPred with clust_5_CGsubsetted
clust_5_CGsubsetted_merged_noPE <- merge(clust_5_CGsubsetted_noPE, lassoPred[, c("baitID", "CGIass")], by = "baitID", all.x = TRUE)

############################
## Compute proportion per cluster
############################
summary_df_noPE <- clust_5_CGsubsetted_merged_noPE %>%
  group_by(Cluster_Name) %>%
  summarise(
    Proportion_with_CGI = mean(CGIass, na.rm = TRUE)
  )

############################
## Plot
############################
pdf(file.path(working_dir, "SuppFigureS2E.pdf"))
ggplot(summary_df_noPE, aes(x = Cluster_Name, y = Proportion_with_CGI)) +
  geom_bar(stat = "identity", fill = "gray", color = "black", width = 0.7) +
  labs(
    x = "Cluster",
    y = "Proportion of PEs with CpG islands",
    title = "Proportion of Baits with CGI per Cluster"
  ) + 
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

############################
## Statistics
############################
pos = tapply(clust_5_CGsubsetted_merged_noPE$CGIass, clust_5_CGsubsetted_merged_noPE$Cluster_Name, sum)
neg = tapply(clust_5_CGsubsetted_merged_noPE$CGIass, clust_5_CGsubsetted_merged_noPE$Cluster_Name, function(x)length(x)-sum(x))
cont_table = matrix(c(pos, neg), nrow=2, byrow=T)
chisq.test(cont_table)$p.value # p-value  = 2.273276e-27


################################################################################
### Supplemental Figure S2F - Spreading of H3K27me3 around PEs involved in different temporal contact dynamics.
################################################################################
## Load the temporal contact clusters
temporalContactClusters <- fread(Temporalcluster_file)
head(temporalContactClusters)
dim(temporalContactClusters) # 9998   13

#### Merge the H3K27me3 spreading at sorted PEs with the contact information (containing the cluster dynamics)
contacts_merged <- temporalContactClusters %>%
  left_join(sorted_regions_short_baitID_primed,
            by ="baitID")

head(contacts_merged) # 9998   18


# Boxplot
pdf(file.path(working_dir, "SuppFigureS2F_Boxplots_H3K27me3_spreading_VS_contact_classes.pdf"))
ggplot(contacts_merged[contacts_merged$merged_dtw_trans_cl6_2with6!=''], aes(x = merged_dtw_trans_cl6_2with6, y = width_kb_noPE)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  ylab("H3K27me3 spreading around PEs") +
  xlab("PE temporal class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

kruskal.test(width_kb_noPE ~ merged_dtw_trans_cl6_2with6, data = contacts_merged) # p-value = 1.593e-12



################################################################################
### Figure 4C - Heatmap showing the enrichment (or depletion) of PE chromatin clusters for each class of temporal contact dynamics 
################################################################################

############################
## The function
############################
# This function, clustEnrichment, is used to assess and visualize the enrichment or association between two clustering assignments 
# using a contingency table and Fisher's exact test
# The function takes a dataset x and two clustering column names (clust1col and clust2col). It performs the following analyses:
# 1. Construct a Contingency Table: tbl that counts occurrences of pairs of cluster labels.
# 2. Perform Chi-squared test on the Whole Table: tests the overall association between the two clusterings using chisq.test on the entire table (ct), storing the p-value in chisq_p and the standardized residuals in stdres
# 3. Calculate Log-Odds Ratios (LOR)
# 4. Plot a Heatmap (if plotHeatmap is TRUE): Colors represent the log-odds ratio (LOR) and asterisks (*) mark significant associations (stdres > 1.96)
# 5. Return Results: the function returns a list with:
# ct: The contingency table
# stdres: standardized residuals
# lor: The log-odds ratio matrix

clustEnrichment = function(x, clust1col="Cluster_Name", clust2col = "merged_dtw_trans_cl6_2with6", plotHeatmap=T, plotResiduals=F){
  # Create contingency table
  tbl = table(x[, get(clust1col)], x[, get(clust2col)])
  ct = as.matrix(tbl)
  
  # Pearson's Chi-squared test
  chisq <- chisq.test(ct)
  stdres <- chisq$stdres  # standardized residuals
  chisq_p <- chisq$p.value
  
  if(plotResiduals==TRUE){
    lor_df <- as.data.frame(as.table(stdres))
  }else{
    
    # Compute LOR
    lor = ct
    for(i in 1:nrow(ct)){
      for(j in 1:ncol(ct)){
        num1 = ct[i,j]
        denom1 = sum(ct[i, -j])
        num2 = sum(ct[-i,j])
        denom2 = sum(ct[-i,-j])
        lor[i,j] = log(num1/denom1) - log(num2/denom2)
      }
    }
    lor_df <- as.data.frame(as.table(lor))  # Works for matrices
    
  }
  
  
  if(plotHeatmap){
    # Convert matrices to data frame
    colnames(lor_df) <- c(clust1col, clust2col, "LOR")
    stdres_df <- as.data.frame(as.table(stdres))
    colnames(stdres_df) <- c(clust1col, clust2col, "stdres")
    lor_df$stdres <- stdres_df$stdres
    lor_df$signif <- ifelse(abs(lor_df$stdres) > 1.96, "*", "") 
    
    head(lor_df)
    # note that the colour key legend will say "LOR" irrespectively of whether LOR or stdres are plotted on the heatmap
    g=ggplot(lor_df, aes_string(x = clust1col, y = clust2col, fill = "LOR")) +
      geom_tile() +
      geom_text(aes(label = signif), color = "black", size = 8) + 
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(title = paste("Chi-squared test p-value =", sprintf("%.2e", chisq_p),")"),
           x = "Chromatin clusters of PEs",
           y = "Contact classes during the transition",
           fill = "LOR") +
      theme_minimal()
    print(g)
  }
  
  if (plotResiduals){
    res = list(ct=ct, stdres = stdres, chisq_p=chisq_p)
  }else{
    res = list(ct=ct, lor=lor, stdres = stdres, chisq_p=chisq_p)
  }
  
  res
  
}

############################
## Data preparation
############################
# The temporal contact clusters 
temporalContactClusters <- fread(Temporalcluster_file)
head(temporalContactClusters)

# The chromatin PEs clusters
baitclust <- clust_5_full_Wkb_noPE
head(baitclust)

# Merge the two clusters
pepm_total <- merge(temporalContactClusters, baitclust, by="baitID")
pepm_total_filtered <- pepm_total[,c("baitID", "oeID", "merged_dtw_trans_cl6_2with6", "Cluster_Name")]

# Set desired Y-axis order
pepm_total_filtered[, merged_dtw_trans_cl6_2with6 := factor(
  merged_dtw_trans_cl6_2with6,
  levels = rev(c("Constant", "Gained early", "Gained late", "Transient", "Lost"))
)]

# Set desired X-axis order
pepm_total_filtered[, Cluster_Name := factor(
  Cluster_Name,
  levels = c("Broad H3K27me3", "High Occup", "CTCF + TFs", "Intermediate", "Low Occup")
)]

############################
## Plot heatmap
############################
pdf(file.path(working_dir, "Figure4C_ClustEnrichment.pdf"), width = 7, height = 5)
res <- clustEnrichment(pepm_total_filtered, "Cluster_Name", "merged_dtw_trans_cl6_2with6") # temporal clustering without P, NE, DE and threshold, new clustering
dev.off()

############################
## Inspect statistical output
############################
## Chi test
chisq <- chisq.test(res$ct) # p-value = 1.972e-14

# Pearson residuals
chisq$residuals

# Standardized residuals
chisq$stdres

