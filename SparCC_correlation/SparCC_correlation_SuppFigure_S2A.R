################################################################################
### Sparse correlation analysis of PECHi-C signals - Supplemental Figure S2A
### Author: Dr. Marina Nocente
###
### This script reproduces Supplemental Figure S2A:
###  Heatmap showing SparCC values (sparse correlations for composition data) 
###  for pairwise relationships between PECHi-C contact scores in transitioning and primed cells
###  The grey colour indicates a correlation of 1
################################################################################

########################################
### Load required libraries
########################################
library(data.table)
library(reshape2)
library(ggplot2)

########################################
### User configuration
########################################
base_dir <- "~/Documents/Bioinformatics"
working_dir <- file.path(base_dir, "/Final_scripts_paper/SparCC_correlation")

peakMatrix <- file.path(base_dir, "/Final_scripts_paper/OSF_data/Focused_peakmatrix_PE_TSS_interaction.tsv")
spoutT_file <- file.path(working_dir, "sparcc_n_p_transition_corr.txt")

########################################
### Load PECHi-C contact scores data
########################################
pepm <- fread(peakMatrix)

# Keep the naive to day14 transition and primed contact scores
scores <- pepm[,3:10]
head(scores)

########################################
### Calculation of the SparCC values
### Done on a cluster because SpiecEasi didn't install on Mac
########################################
# library(SpiecEasi)
# spoutT = sparcc(scores)
# write.table(spoutT$Cor, "sparcc_n_p_transition_corr.txt", quote=F, sep="\t")

########################################
### Load SparCC values for naive-to-day 14 and primed time points
########################################
spoutT <- fread(spoutT_file, drop = 1)
head(spoutT)

spoutT = as.matrix(spoutT)
rownames(spoutT) = colnames(spoutT) = colnames(scores)

spoutT_long <- melt(spoutT)

########################################
### Plot SparCC values
########################################
pdf(file.path(working_dir, "SuppFigure_S2A.pdf"))
ggplot(spoutT_long, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +   # colored squares with white borders
  scale_fill_gradient2(low = "white", high = "red", 
                       limit = c(-0.01, max(spoutT[spoutT!=1])), space = "Lab",
                       name="SparCC\ncorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1)) +
  coord_fixed() +
  labs(x = "", y = "")
dev.off()
