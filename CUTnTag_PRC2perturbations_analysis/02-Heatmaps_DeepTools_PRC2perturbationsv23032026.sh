
#PBS -S /bin/bash
#PBS -N deepTools
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=80:mem=80gb
#PBS -o /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/Run2_June2025/computeMatrix.out
#PBS -e /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/Run2_June2025/computeMatrix.err


source /rds/general/user/mnocente/home/anaconda3/etc/profile.d/conda.sh
conda activate deeptools_env

cd /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/Run2_June2025/Calibrated_Bigwig_from_Andrew


###### Step 1: Prepare the regions of interest

## I will use my poised enhancer coordinates (baimtap file) (32142 PE regions)
# This file is available on GitHub "Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap"

## Check how are "called" the chr in the bigwig file (ex chr1 or 1)
#multiBigwigSummary bins -b ./Calibrated_Bigwig_from_Andrew/MergedBigwigs/DMSO_H2AK119Ub.merged.bw -o test_results.npz
# result: Number of bins found: 308839
# so it's "chr"

## Convert the chromosome name in the baitmap file
#awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' /Users/mnocente/Documents/Bioinformatics/Bioinfo_REDO_figure_Monica_CnT/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap > /Users/mnocente/Documents/Bioinformatics/CnT_PRC2_perturb_Day5/Run2_CUTnTagPrc2_June2025/regions_chr_PE_baits.bed



###### Step 2: Define variables

# Set variables
BED_FILE="regions_chr_PE_baits.bed"
BIGWIG_DIR="./MergedBigwigs"
OUTPUT_MATRIX_K27me3="matrix_PE_K27me3.gz"
OUTPUT_MATRIX_K119ub="matrix_PE_K119ub.gz"
THREADS=8

# List of bigWigs
BIGWIGS_K27me3=(
    "$BIGWIG_DIR/DMSO_H3K27me3.merged.bw"
    "$BIGWIG_DIR/PROTAC_H3K27me3.merged.bw"
    "$BIGWIG_DIR/UNC_H3K27me3.merged.bw"
)

BIGWIGS_K119ub=(
    "$BIGWIG_DIR/DMSO_H2AK119Ub.merged.bw"
    "$BIGWIG_DIR/PROTAC_H2AK119Ub.merged.bw"
    "$BIGWIG_DIR/UNC_H2AK119Ub.merged.bw"
)

###### Step 3: Compute the signal matrix

# Run computeMatrix for H3K27me3
computeMatrix reference-point \
    -S "${BIGWIGS_K27me3[@]}" \
    -R $BED_FILE \
    --referencePoint center \
    -a 10000 -b 10000 \
    -out $OUTPUT_MATRIX_K27me3 \
    --skipZeros \
    -p $THREADS

# Run computeMatrix for H2AK119ub
computeMatrix reference-point \
    -S "${BIGWIGS_K119ub[@]}" \
    -R $BED_FILE \
    --referencePoint center \
    -a 10000 -b 10000 \
    -out $OUTPUT_MATRIX_K119ub \
    --skipZeros \
    -p $THREADS


# Notes:
# --referencePoint center: aligns each region at its midpoint
# -a 10000 -b 10000 : 10 kb upstream and downstream of the center.
# -R: my PEs
# -S: my calibrated bigWigs

# computeMatrixOperations info -m matrix_PE.gz
# zcat matrix_PE.gz | grep -v "^@" | wc -l # to find the number of PEs represented on the heatmaps


###### Step 4: Plot heatmaps and profiles

# Heatmap for H3K27me3
plotHeatmap \
    -m matrix_PE_K27me3.gz \
    -out heatmap_PE_K27me3.pdf \
    --colorMap Reds \
    --zMin 0 \
    --zMax 5 \
    --heatmapHeight 10 \
    --heatmapWidth 5 \
    --samplesLabel "DMSO_H3K27me3" \
                   "PROTAC_H3K27me3" \
                   "UNC_H3K27me3" \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions descend

# Heatmap for H2AK119ub
plotHeatmap \
    -m matrix_PE_K119ub.gz \
    -out heatmap_PE_K119ub.pdf \
    --colorMap Reds \
    --zMin 0 \
    --zMax 10 \
    --heatmapHeight 10 \
    --heatmapWidth 5 \
    --samplesLabel "DMSO_H2AK119Ub" \
                   "PROTAC_H2AK119Ub" \
                   "UNC_H2AK119Ub" \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions descend


# Profile for H3K27me3
plotProfile \
    -m matrix_PE_K27me3.gz \
    -out profile_PE_K27me3.pdf \
    --perGroup \
    --colors blue purple black

# Profile for H2AK119ub
plotProfile \
    -m matrix_PE_K119ub.gz \
    -out profile_PE_K119ub.pdf \
    --perGroup \
    --colors red orange gray

