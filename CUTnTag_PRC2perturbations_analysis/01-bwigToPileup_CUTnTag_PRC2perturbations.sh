#PBS -S /bin/bash
#PBS -N bigwig2PileUp
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -o /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations
#PBS -e /rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations

module load anaconda3/personal

source activate bwtool


#### Define paths ####
BigWiginputPATH="/rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/bigwiv"
PileUpoutputPATH="/rds/general/user/mnocente/home/analysis/Marina/CnTag/CnT_PRC2_perturbations/PileUp_from_bigwig"

# Note: The bigwig files are already calibrated (see the Method section "PRC2 perturbation experiments - CUT&Tag data") 

# Create output directories if they don’t exist
mkdir -p ${PileUpoutputPATH}


#### List of samples ####
samples=("DMSO_Rep1_H3K27me3" "PROTAC_Rep1_H3K27me3" "PROTAC_Rep2_H3K27me3" "UNC_Rep1_H3K27me3" "UNC_Rep2_H3K27me3"
    "lane8849_TAAGGCGA_TAGATCGC_DMSO_Rep1_H3K27me3_L001" "lane8849_GCTACGCT_AGAGTAGA_DMSO_Rep2_H3K27me3_L001"
    "lane8849_GCTACGCT_TATCCTCT_UNC_Rep1_H3K27me3_L001" "lane8849_TGCTGGGT_ACTGCATA_UNC_Rep2_H3K27me3_L001"
    "lane8849_TGCTGGGT_GTAAGGAG_PROTAC_Rep2_H3K27me3_L001")


#### Convert Bigwig to PileUp ####
echo "Converting Bigwig to PileUp"
for sample in "${samples[@]}"; do
    bigwig_file="${BigWiginputPATH}/${sample}.normalized.bw"
    PileUp_file="${PileUpoutputPATH}/${sample}.pileup"

    echo "Processing $sample..."

    bwtool summary /rds/general/user/mnocente/home/analysis/MDR/DesignDir/Human_PEcHiC_DpnII_bin5k_sol_baits/Human_PEcHiC_hg38_DpnII_bin5K_sol_baits.baitmap_4col.txt "$bigwig_file" "$PileUp_file" -header -with-sum
    
done

echo "Pipeline completed successfully!"
conda deactivate
