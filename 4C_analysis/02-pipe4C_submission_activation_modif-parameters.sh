#PBS -S /bin/bash
#PBS -N pipeC
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=80:mem=80gb
#PBS -o /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/pipe4C/modified_parameters/
#PBS -e /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/pipe4C/modified_parameters/


module load anaconda3/personal
source activate 4C_env

cd /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/pipe4C/modified_parameters

echo "start pipeC R script"

Rscript /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Aug2024/pipeC/pipe4C.R --vpFile /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/pipe4C/vpFile.tsv --fqFolder /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/fastp --outFolder /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/pipe4C/modified_parameters --cores 8 --wig --plot --genomePlot --mapUnique TRUE --qualityCutoff 1 --tsv


echo "done"

conda deactivate


### Notes:
# The pipe4C.R script (Krijger et al. 2020) is available on GitHub:https://github.com/deLaatLab/pipe4C
# Our vpFile.tsv file is available on our GitHub.
# Minimum mapping quality of 1, retention of uniquely mapped reads, and cis-only analysis
