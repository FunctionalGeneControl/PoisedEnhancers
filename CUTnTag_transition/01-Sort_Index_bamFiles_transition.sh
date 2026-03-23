#PBS -S /bin/bash
#PBS -N bedTobedgraphTobigwig
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -o /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_TOsort_index.out
#PBS -e /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_TOsort_index.err


### Load environment
module load anaconda3/personal
source activate CutNtag  # Need to contain bedtools and samtools

### Define paths
BAMinputPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95"
BAMprocessedOutputPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_sorted_indexed_combined"

### Genome sizes file for bedtools and bedGraphToBigWig
mygenome="/rds/general/project/lms-spivakov-analysis/live/MDR/Cut_N_Tag/tests/SEACR/human.hg38.genome.noChr"

### Move to BED input directory
cd "$BAMinputPATH" || exit 1

echo "Current working directory: $(pwd)"


### Main loop over BAM files
for bam in "$BAMinputPATH"/*.bam; do
    [ -e "$bam" ] || continue  # skip if no BAM files

    base=$(basename "$bam" .bam)

    echo "Processing $bam"

    ## 1. Sort by coordinate
    samtools sort -@ 4 -o "$BAMprocessedOutputPATH/${base}.sorted.bam" "$bam"

    ## 2. Index
    samtools index "$BAMprocessedOutputPATH/${base}.sorted.bam"

    echo "Sorted and indexed BAM: $BAMprocessedOutputPATH/${base}.sorted.bam"
done

echo "All BAMs sorted and indexed."


