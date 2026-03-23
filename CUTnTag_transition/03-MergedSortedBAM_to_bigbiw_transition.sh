#PBS -S /bin/bash
#PBS -N MergedSortedBam_to_bigbig
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -o /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/MergedSortedBam_to_bigbig.out
#PBS -e /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/MergedSortedBam_to_bigbig.err


### Load environment
module load anaconda3/personal
source activate CutNtag  # Need to contain bedtools and samtools

### Define paths
MergedSortedBamPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_sorted_indexed_combined"
bigwigOutputPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bigwig/bigwig_merged_just_sorted_indexedBamFiles"


### Move to BAM input directory
cd "$MergedSortedBamPATH" || exit 1

echo "Current working directory: $(pwd)"


### Convert merged bam files in bigwig
for f in merged_*.bam
do
  echo "Converting $f"
  bamCoverage -b "$f" -o "$bigwigOutputPATH/$f.bw" --numberOfProcessors 30
done

echo "All the bigwig files are done"
