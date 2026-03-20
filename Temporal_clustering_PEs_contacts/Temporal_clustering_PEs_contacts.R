################################################################################
### Temporal classes from contacts
### by Dr. Marina Nocente
###
### Investigate dynamic changes in chromatin contacts during the HNES1 naive-to-primed transition experiment
### Figure 2B: Dynamics of PE chromosomal contacts with gene promoters and other PEs during the naive-to-primed transition.
### Figure 2D: UpSet plot showing the numbers of PEs involved in different classes of contact dynamics
################################################################################

#################################
## Load the useful libraries
#################################
library(data.table)
library(dtwclust)
library(ggplot2)
library(UpSetR)

#################################
## Paths
#################################
working_dir <- file.path("~/Documents/Bioinformatics/Final_scripts_paper")

pepm_transition_clusters_files <- file.path(working_dir, "Temporal_clustering_PEs_contacts/pepm_transition_dtwCluster_Focused_Oct2025.txt")

#################################
## Run the following R script for the Dynamic Time Warping (DTW) clustering based on PECHi-C contact scores
####################################
## dtw_clust_Transition_MonicaSelection_finalRedo_Oct2025.R

#################################
## Upload the clustering results and create the final clustering of PECHi-C contacts
####################################
pepm_transition <- fread(pepm_transition_clusters_files)

# Merge clusters 2 and 6 (they have a similar behavior) in the clustering "dtw_trans_cl6" (containing initially 6 clusters) and create a new clustering "merged_dtw_trans_cl6_2with6"
pepm_transition[, merged_dtw_trans_cl6_2with6 := ifelse(dtw_trans_cl6 %in% c(2, 6), "2_6", as.character(dtw_trans_cl6))]

## Order and select columns
pepm_transition_final <- pepm_transition[, .(baitID, oeID, naive, day1, day3, day5, day7, day10, day14, primed, DE, NE, merged_dtw_trans_cl6_2with6)]
head(pepm_transition_final)
dim(pepm_transition_final) # 9998   13

pepm_transition_final[, merged_dtw_trans_cl6_2with6 := fcase(
  merged_dtw_trans_cl6_2with6 == 1, "Gained early",
  merged_dtw_trans_cl6_2with6 == "2_6", "Gained late",
  merged_dtw_trans_cl6_2with6 == 3, "Constant",
  merged_dtw_trans_cl6_2with6 == 4, "Lost",
  merged_dtw_trans_cl6_2with6 == 5, "Transient"
)]

fwrite(pepm_transition_final, file.path(working_dir, "Temporal_clustering_PEs_contacts/pepm_transition_dtw_final_02-02-2026.txt"), sep="\t")


#################################
## Function to plot the different temporal clusters
#################################
## Define cluster color:
custom_colors <- c("Constant" = "green3", "Gained early" = "red", "Gained late" = "orange", "Transient" = "purple", "Lost" = "blue")

## Define Time labels:
custom_time_labels <- c("1" = "naive", "2" = "day1", "3" = "day3", "4" = "day5", "5" = "day7", "6" = "day10", "7" = "day14", "8" = "primed", "9" = "DE", "10" = "NE")

## Run the function
plotTempClust = function(pepm_tot, clustcol, scorecols=3:10, plotBoxplots=T, plotMeans=T, cluster_colors=NULL, time_labels=NULL, y_limits=NULL){
  if(!clustcol %in% names(pepm_tot)){ 
    message("clustcol=", clustcol, " doesn't exist in the peak matrix")
    stop() 
  }
  
  if(plotBoxplots){
    pepm_long <- tidyr::pivot_longer(pepm_tot,   
                                     cols = scorecols, 
                                     names_to = "Time_Point", 
                                     values_to = "Value")
    pepm_long$Time_Point = factor(pepm_long$Time_Point, levels=names(pepm_tot)[scorecols])
    cl = vector("list")
    
    for(clust in sort(unique(pepm_tot[[clustcol]]))){
      cl[[clust]] = ggplot(pepm_long[pepm_long[[clustcol]]==clust,], 
                           aes(x = Time_Point, y = Value, fill = factor(!!sym(clustcol)))) +
        geom_boxplot(color = "black") + # Keep black borders for visibility
        stat_summary(fun = "median", geom = "point", color = "black", size = 2) +
        labs(title = paste("Cluster", clust), x = "Time Points", y = "Value") +
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis text
      
      # Apply custom colors if provided
      if(!is.null(cluster_colors)){
        cl[[clust]] = cl[[clust]] + scale_fill_manual(name = "Cluster", values = cluster_colors)
      }
      
      # Apply custom x-axis labels if provided
      if(!is.null(time_labels)){
        cl[[clust]] = cl[[clust]] + scale_x_discrete(labels = time_labels)
      }
      
      # Apply fixed Y-axis limits if provided
      if (!is.null(y_limits)) {
        cl[[clust]] = cl[[clust]] + coord_cartesian(ylim = y_limits)
      }
      
    }
    
    message("Plotting the boxplot")
    print(cowplot::plot_grid(plotlist=cl, align="v", ncol=1))
  }
  
  if(plotMeans){
    nclust = length(unique(pepm_tot[[clustcol]]))
    means = pepm_tot[, apply(.SD,2,mean), by=clustcol, .SDcols = scorecols]
    means[, timepoint:=rep(1:length(scorecols), nclust)]
    means[, timepoint:=as.factor(timepoint)]
    means[, c(clustcol):=as.factor(get(clustcol))]
    
    g = ggplot(means, aes_string(x="timepoint", group=clustcol, colour=clustcol, y="V1")) + 
      geom_line() +  
      geom_point(size = 2) +
      labs(y = "Mean of the contact score per category of contacts engaging PEs") +
      theme_minimal(base_size = 14)  
    
    # Apply the same custom colors
    if(!is.null(cluster_colors)){
      g = g + scale_color_manual(name = "Cluster", values = cluster_colors)
    }
    
    # Apply custom x-axis labels if provided
    if(!is.null(time_labels)){
      g = g + scale_x_discrete(labels = time_labels)
    }
    
    message("Plotting the means")
    print(g)
  }
}


#################################
## Plot the different temporal clusters of PE contacts between naive and day14 - Figure 2B
#################################
plotTempClust(pepm_transition_final, "merged_dtw_trans_cl6_2with6", plotBoxplots = F, plotMeans=T, 
              scorecols = 3:9, cluster_colors=custom_colors, time_labels=custom_time_labels)

# Number of contact per class
pepm_transition_final[ , .N, by = merged_dtw_trans_cl6_2with6]

#merged_dtw_trans_cl6_2with6     N
#                    Constant  1999
#                   Transient  2847
#                        Lost  2326
#                 Gained late  2045
#                Gained early   781


## Plot with gaps between the time points
long_filtered_pepm <- melt(
  pepm_transition_final,
  measure.vars = 3:9,
  variable.name = "timepoint",
  value.name = "value"
)

# Map your timepoint names to numeric time values
time_map <- data.table(
  timepoint = c("naive", "day1", "day3", "day5", "day7", "day10", "day14"),
  time = c(0, 1, 3, 5, 7, 10, 14)
)

long_filtered_pepm <- merge(long_filtered_pepm, time_map, by = "timepoint")


#pdf(file.path(working_dir, "Figure2B_Contact_temporal_clusters_NtoP_gaps_between_timepoints.pdf"))
ggplot(
  long_filtered_pepm, 
  aes(x = time, y = value, color = merged_dtw_trans_cl6_2with6, group = merged_dtw_trans_cl6_2with6)
) +
  # Mean line with thinner weight
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  
  # Mean ± SE error bars
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, size = 0.8) +
  
  # Color/fill scale
  scale_color_manual(values = custom_colors, name = "Clusters") +
  scale_fill_manual(values = custom_colors, guide = "none") +
  
  # Proper time axis
  scale_x_continuous(
    breaks = c(0, 1, 3, 5, 7, 10, 14), 
    labels = c("naive", "day 1", "day 3", "day 5", "day 7", "day 10", "day 14")
  ) +
  
  # Axis labels
  labs(
    y = expression("Mean of arcsinh-transformed contact scores per category of contacts engaging PEs"), 
    color = "Clusters"
  ) +
  
  # Clean theme with big text
  theme_classic(base_size = 14) +
  
  # Finer theme tweaks
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90")
  )

#dev.off()


########################################
### UpsetPlot - Figure 2D
########################################

pepm_transition_final_filtered <- pepm_transition_final[,c("baitID", "oeID", "merged_dtw_trans_cl6_2with6")]

membership_matrix <- dcast(
  pepm_transition_final_filtered[, c("baitID", "merged_dtw_trans_cl6_2with6")],
  baitID ~ merged_dtw_trans_cl6_2with6, 
  # the final table should have baitID and each value of merged_dtw_trans_cl6_2with6 as separate columns
  value.var = "merged_dtw_trans_cl6_2with6", # look at the value of merged_dtw_trans_cl6_2with6 for each baitID
  fun.aggregate = function(x) 1L, 
  # aggregate over all merged_dtw_trans_cl6_2with6 for each baitID
  # and return the presence of a cluster for a given baitID as 1
  fill = 0L # absence of a cluster for a given baitID = 0
) 

rownames(membership_matrix) <- membership_matrix$baitID
membership_matrix_bin = membership_matrix[, -1]
membership_matrix_bin = membership_matrix_bin[, c("Lost", "Transient", "Gained late", "Gained early", "Constant")]

upset(
  membership_matrix_bin,
  nsets = ncol(membership_matrix_bin), # number of classes to show
  empty.intersections = "on", 
  sets = names(membership_matrix_bin),
  keep.order = T, order.by="degree", decreasing=F,
  sets.bar.color = "#56B4E9",
  main.bar.color = "#D55E00",
  matrix.color = "#009E73"
)

pdf("Upsetplot_contacts_rewiring.pdf")
upset(
  membership_matrix_bin,
  nsets = ncol(membership_matrix_bin), # number of classes to show
  empty.intersections = "on", 
  sets.bar.color = "#56B4E9",
  main.bar.color = "#D55E00",
  matrix.color = "#009E73"
)
dev.off()
