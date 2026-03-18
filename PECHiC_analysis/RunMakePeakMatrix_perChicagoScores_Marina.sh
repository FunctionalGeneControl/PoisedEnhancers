#PBS -S /bin/bash
#PBS -N makePeakMatrix_single
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=32:mem=96gb

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

module add anaconda3/personal
source activate chicago

PATH=$PATH:/rds/general/project/lms-spivakov-analysis/live/chicago/chicagoTools

outprefix=${PBS_O_WORKDIR}/peakMatrices/PerChicagoScores/${prefix}

Rscript /rds/general/project/lms-spivakov-analysis/live/chicago/chicagoTools/makePeakMatrix.R ${listfile} ${outprefix} --scorecol score --cutoff 5 > makePeakMatrix_single.${PBS_JOBID}.log 2>&1

source deactivate

## This script is using another script "makePeakMatrix.R" from ChicagoTools with default parameters and avaiable here:
## https://bitbucket.org/chicagoTeam/chicago/src/master/chicagoTools/
## The listfile is available on our GitHub