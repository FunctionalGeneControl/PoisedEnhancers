################################################################################
### Script for Supplementary Figure S1
### Author: Dr. Marina Nocente
###
### This script reproduces Supplemental Figure S1 panels showing:
### - Number of hPSC PEs found in the “active” state in human adult tissues according to the ENCODE project. 
### - Number of PEs with in vivo tissue-specific activity overlapping poised enhancers according to VISTA Enhancer Browser, grouped by annotated tissue activity.
### - Enrichment of PE-contacted regions for histone marks versus distance-matched random regions in naive (yellow bars) and primed cells (blue bars). 
################################################################################

#######################################
## Libraries
#######################################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

#######################################
## File paths
#######################################
base_dir <- "~/Documents/Bioinformatics"
annotation_dir <- file.path(base_dir, "Annotation_enhancers/")
chicagoEnrich_dir <- file.path(base_dir, "Enrichment_CHiCAGO/")

baitmap_file <- file.path(base_dir, "Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap")
decoration_matrix_ENCODE_file <- file.path(annotation_dir, "cCRE_decoration.matrix.1")
encode_UCSC_file <- file.path(annotation_dir, "encodeCcreCombined.bed")
vista_enhancers_file <- file.path(annotation_dir, "vista_locus.tsv")

#######################################
## Load baitmap (Poised Enhancers (PE) coordinates)
#######################################
baitmap <- fread(baitmap_file)

setnames(baitmap, c("V1", "V2", "V3", "V4"), c("chr", "bait_start", "bait_end", "baitID"))
baitmap[, chr := gsub("chr", "", chr)]
setkey(baitmap, chr, bait_start, bait_end)

#######################################
## Supplemental Figure S1B
#######################################

### En-tex 'enhancer decoration matrix' from an ENCODE subproject
# https://downloads.wenglab.org/cCRE_decoration.matrix.1.gz
# Rows are cCREs (candidate cis-Regulatory Elements).
# Columns: One for the ID and hundreds of tissue-specific annotations (like whether the enhancer is active in a given tissue).

### Load enhancer decoration matrix
edec <- fread(decoration_matrix_ENCODE_file)

# filter the fields for "active" and the rows for at least one "1" in such a field
aedec = edec[, c(1, grep("active", names(edec))), with=F]
ncol(edec) # 723
ncol(aedec) # 337 # From 723 total columns, only 337 are related to active enhancer marks

# Filter rows with at least one active signal
aedecpos = apply(aedec[, -1], 1, sum)
aedecf = aedec[which(aedecpos>0),] 
nrow(aedec) # 890906
nrow(aedecf) # 635906 # ~635k cCREs have at least one "active" tissue call


### Merge with enhancer coordinate annotations
# Note: the matrix edec only has IDs, not genomic positions
# Use the encode v3 track at UCSC (downloaded as bigBed, converted using bigBedToBed)
enc_ucsc = fread(encode_UCSC_file)
nrow(enc_ucsc) # 926535

nrow(merge(aedecf, enc_ucsc, by.x="cCRE_id", by.y="V4"))  # 635906
# All 635k active cCREs were matched to their genomic coordinates.

# Now, enc_enh contains enhancer IDs and their genomic coordinates (chr, start, end).
enc_enh = merge(aedecf, enc_ucsc[, 1:4], by.x="cCRE_id", by.y="V4")
setnames(enc_enh, c("V1", "V2", "V3"), c("chr", "enh_start", "enh_end"))
enc_enh[, chr:=gsub("chr", "", chr)]

### Find overlaps between enhancers and baits
encanno = foverlaps(enc_enh, baitmap, by.x=c("chr","enh_start", "enh_end"), by.y=c("chr", "bait_start", "bait_end"), nomatch=NULL)

# Check how many baits and overlaps were found
nrow(baitmap) # 32142 baits
nrow(encanno) # 46819 overlaps

length(unique(encanno$baitID)) # 21950 unique baits with at least one enhancer overlap
21950/32142 #  0.6829071 (~68% of baits are associated with at least one enhancer)- could be better but good enough

# Identify tissue names from column names
tissues = unique(gsub("\\S+\\-(\\S+)","\\1", names(encanno)))[-c(1:7)]
tissues = tissues[!tissues%in%c("enh_start", "enh_end")] # 25

# Create binary indicators for enhancer activity in each tissue
for(tt in tissues){
  encanno[, c(tt):=apply(.SD,1,function(x)sum(x)>0), .SDcols = grep(tt, names(encanno))]
}

# Count in how many tissues each enhancer is active
enhm = encanno[, c(1:4,6, 343:ncol(encanno)), with=F]
enhm[, no_of_tissues:=apply(.SD,1,sum), .SDcols=8:ncol(enhm)]

### Distribution of enhancer activity across tissues
pdf(file.path(annotation_dir, "SuppFigure_S1B.pdf"), width = 8, height = 7)
hist(enhm$no_of_tissues, 
     xlim = c(0, 28), 
     ylim = c(0, 14000), 
     xaxt = "n",             # suppress default X-axis
     xlab = "Number of Tissues", 
     ylab = "Frequency")

# Add custom X-axis ticks, e.g., every 1 tissue
axis(side = 1, at = seq(0, 28, by = 1))
dev.off()


### Calculations
# Number of enhancers active in one human adult organ:
sum(enhm$no_of_tissues == 1) # 4455

# Number of enhancers active in 27 human adult organs:
sum(enhm$no_of_tissues == 27) # 13487

# How many enhancers are active in at least one tissue?
sum(enhm$no_of_tissues > 0) # 46260



#######################################
## Supplemental Figure S1C
#######################################

### Load and filter VISTA enhancers
vista <- fread(vista_enhancers_file)

# Filter to only positive enhancers (curation_status=="positive") with human IDs (^hs)
vistapos <- vista[curation_status=="positive"] # 2117
vistaposh <- vistapos[grep("^hs", vista_id)] # 1221

# Parse coordinate_hg38 into columns
vistaposh[, c("chr", "start", "end"):=tstrsplit(coordinate_hg38, "\\:|\\-")]

# Strip "chr" prefix and convert coordinates to numeric for overlap
vistaposh[, chrnochr:=gsub("chr", "", chr)]
vistaposh[, start:=as.numeric(start)]
vistaposh[, end:=as.numeric(end)]

### Overlap VISTA enhancers with baitmap
baitmap_vista = foverlaps(vistaposh, baitmap, nomatch = NULL, by.x=c("chrnochr", "start", "end"), by.y=c("chr", "bait_start", "bait_end")) # 261
sort(table(baitmap_vista$tissue)) # Summarize number of overlaps by tissue.

# total number of distinct VISTA enhancers overlapping (regardless of tissue):
length(unique(baitmap_vista$vista_id)) # 133

### VISTA tissue statistics
# Get a list of unique tissue names
vista_tissues = paste(unique(vista$tissue), collapse=";")
vista_tissues1 = unlist(strsplit(vista_tissues, ";"))
vista_tissues_u = unique(vista_tissues1)
vista_tissues_u = vista_tissues_u[vista_tissues_u!=""]


pe_vtm = vector("numeric")
for(tissue in vista_tissues_u){
  pe_vtm[tissue] = length(grep(tissue, baitmap_vista$tissue))
}
sort(pe_vtm, decreasing = T)

# total number across tissues (one enhancer can be active in multiple tissues):
sum(pe_vtm) # 647
sum(unique(pe_vtm)) # 638

## Named vector: abbreviation -> full name
tissue_names <- c(
  ba = "Branchial arch",
  hb = "Hindbrain",
  eye = "Eye",
  fb = "Forebrain",
  nt = "Neural tube",
  mb = "Midbrain",
  nose = "Nose",
  lb = "Limb",
  cn = "Cranial nerve",
  drg = "Dorsal root ganglion",
  ht = "Heart",
  tri = "Trigeminal ganglion",
  som = "Somite",
  gen = "Genital tubercle",
  mel = "Melanocytes",
  tail = "Tail",
  fm = "Facial mesenchyme",
  other = "Other",
  ear = "Ear",
  lv = "Liver",
  bv = "Blood vessels",
  pan = "Pancreas"
)


### plot the count of PECHiC baits overlapping with VISTA enhancers, per tissue

# counts per tissue:
counts_df <- data.frame(
  tissue = names(pe_vtm),
  count = pe_vtm
)

# Add full tissue names
counts_df$tissue_full <- tissue_names[counts_df$tissue]

### final figure
count_color <- "#4E79A7"      # steel blue

theme_genesdev <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

# Create the Counts plot
p_counts <- ggplot(counts_df, aes(x = reorder(tissue_full, count), y = count)) +
  geom_bar(stat = "identity", fill = count_color, width = 0.7) +
  coord_flip() +
  labs(
    title = "Counts of PECHi-C baits overlapping with VISTA enhancers",
    y = "Count"
  ) +
  theme_genesdev

ggsave(
  file.path(annotation_dir, "SuppFigure_S1C.pdf"),
  plot = p_counts,
  width = 12, height = 7
)



####################################
###Supplemental Figure S1D
####################################

### Original datasets (obtained with peakEnrichment4Features() function in CHiCAGO)
df_naive <- data.frame(
  Feature = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K27ac"),
  OLwithSI = c(5226, 1810, 3580, 3703),
  Mean = c(4505.95, 1641.33, 3463.78, 3205.13),
  LowerCI = c(4401.50, 1567.03, 3370.62, 3114.00),
  UpperCI = c(4610.40, 1715.63, 3556.94, 3296.26)
)

df_primed <- data.frame(
  Feature = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K27ac"),
  OLwithSI = c(7392, 1669, 2935, 3490),
  Mean = c(5956, 1293.94, 1926.23, 2845.4),
  LowerCI = c(5847, 1233.68, 1850.46, 2757.07),
  UpperCI = c(6065, 1354.20, 2002.00, 2933.73)
)

# Filter features of interest
features <- c("H3K4me1", "H3K27me3", "H3K4me3")
df_naive <- df_naive %>% filter(Feature %in% features)
df_primed <- df_primed %>% filter(Feature %in% features)

# Add cell type column
df_naive$Cell <- "Naive"
df_primed$Cell <- "Primed"

# Combine datasets
df_combined <- bind_rows(df_naive, df_primed)

# Calculate fold enrichment and proper error bars
df_combined <- df_combined %>%
  mutate(FoldEnrichment = OLwithSI / Mean,
         # Approximate SD of mean from CI
         SD_Mean = (UpperCI - LowerCI) / (2 * 1.96),
         # SE for fold enrichment
         SE_FE = FoldEnrichment * (SD_Mean / Mean),
         # 95% CI for fold enrichment
         FE_Lower = FoldEnrichment - 1.96 * SE_FE,
         FE_Upper = FoldEnrichment + 1.96 * SE_FE)

# Plot fold enrichment with error bars
pdf(file.path(chicagoEnrich_dir, "SuppFigure_S1D.pdf"), width = 7, height = 8)

ggplot(df_combined, aes(x = Feature, y = FoldEnrichment, fill = Cell)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = FE_Lower, ymax = FE_Upper),
                width = 0.2,
                position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("Naive" = "gold", "Primed" = "steelblue")) +
  coord_cartesian(ylim = c(0, 1.7)) + 
  theme_minimal(base_size = 14) +
  labs(y = "Fold enrichment of overlaps with feature (Observed / Random)",
       fill = "Cell type") +
  theme(legend.position = "right")

dev.off()

