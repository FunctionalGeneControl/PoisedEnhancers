#PBS -S /bin/bash
#PBS -N deepTools
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -o /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_transition/computeMatrix.out
#PBS -e /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_transition/computeMatrix.err


#### Run computeMatrix and plotHeatmap for CUT&Tag poised enhancers

source /rds/general/user/mnocente/home/anaconda3/etc/profile.d/conda.sh
conda activate deeptools_env

cd /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_transition

# Set variables
BED_FILE="/rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/Run2_June2025/Calibrated_Bigwig_from_Andrew/regions_chr_PE_baits.bed"
BIGWIG_DIR="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bigwig/bigwig_merged_just_sorted_indexedBamFiles"
OUTPUT_MATRIX="matrix_K27me3_PE.gz"
OUTPUT_HEATMAP="heatmap_K27me3_PE.pdf"
THREADS=8

# List of bigWigs
BIGWIGS=("$BIGWIG_DIR/merged_Primed_273.bam.bw")

# Run computeMatrix
computeMatrix reference-point \
     -S "${BIGWIGS[@]}" \
     -R $BED_FILE \
     --referencePoint center \
     -a 10000 -b 10000 \
     -out $OUTPUT_MATRIX \
     --skipZeros \
     -p $THREADS

plotHeatmap \
    -m ${OUTPUT_MATRIX} \
    -out ${OUTPUT_HEATMAP} \
    --colorMap Blues \
    --zMin 0 \
    --zMax 20 \
    --outFileNameMatrix K27me3_heatmap_matrix.txt \
    --outFileSortedRegions K27me3_heatmap_matrix_sorted_regions.txt \
    --heatmapHeight 10 \
    --heatmapWidth 5 \
    --samplesLabel "merged_Primed_273" \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions descend

