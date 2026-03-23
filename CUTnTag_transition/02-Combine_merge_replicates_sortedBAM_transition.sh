 #PBS -S /bin/bash
#PBS -N Merge_BamSortedFiles
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=30:mem=50gb
#PBS -o /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/Combine_merge_replicates_sortedBAM.out
#PBS -e /rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/Combine_merge_replicates_sortedBAM.err


### Load environment
module load anaconda3/personal
source activate CutNtag  # Need to contain bedtools and samtools

### Define paths
sortedBamPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_sorted_indexed_combined"
mergedBAMoutputPATH="/rds/general/user/mnocente/home/analysis/MDR/Cut_N_Tag/mergedFinal/ftp1.babraham.ac.uk/ftpusr95/bamFiles_sorted_indexed_combined"


### Move to BAM input directory
cd "$sortedBamPATH" || exit 1

echo "Current working directory: $(pwd)"


### Merged H3K27ac replicates
samtools merge -@ 30 merged_Naive_271.bam lane1_N2_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N4_271_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day1_271.bam lane1_B1_271_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C1_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D1_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day3_271.bam lane1_B3_271_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day5_271.bam lane1_B5_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C5_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D5_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day7_271.bam lane1_B7_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C7_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D7_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day10_271.bam lane1_B10_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C10_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D10_271_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day14_271.bam lane1_B14_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C14_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D14_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Primed_271.bam lane1_P2_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P4_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_NE_271.bam lane1_NE1_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE2_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_DE_271.bam lane1_DE1_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE2_271_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE3_271_1_merged_L001_GRCh38_bowtie2.sorted.bam


### Merged H3K27me3 replicates
samtools merge -@ 30 merged_Naive_273.bam lane1_N2_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N4_273_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day1_273.bam lane1_B1_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C1_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D1_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day3_273.bam lane1_B3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day5_273.bam lane1_B5_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C5_273_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D5_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day7_273.bam lane1_B7_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C7_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D7_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day10_273.bam lane1_E10_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C10_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D10_273_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day14_273.bam lane1_E14_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C14_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D14_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Primed_273.bam lane1_P2_273_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P3_273_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P4_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_NE_273.bam lane1_NE1_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE2_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_DE_273.bam lane1_DE1_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE2_273_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE3_273_1_merged_L001_GRCh38_bowtie2.sorted.bam


### Merged H3K4me1 replicates
samtools merge -@ 30 merged_Naive_41.bam lane1_N2_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N3_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N4_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day1_41.bam lane1_B1_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C1_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D1_41_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day3_41.bam lane1_B3_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C3_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D3_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day5_41.bam lane1_B5_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C5_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D5_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day7_41.bam lane1_B7_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C7_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D7_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day10_41.bam lane1_B10_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C10_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D10_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day14_41.bam lane1_B14_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C14_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D14_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Primed_41.bam lane1_P2_41_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P3_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P4_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_NE_41.bam lane1_NE1_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE2_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE3_41_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_DE_41.bam lane1_DE1_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE2_41_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE3_41_1_merged_L001_GRCh38_bowtie2.sorted.bam


### Merged H3K4me3 replicates
samtools merge -@ 30 merged_Naive_43.bam lane1_N2_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_N4_43_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day1_43.bam lane1_B1_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C1_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D1_43_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day3_43.bam lane1_B3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C3_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day5_43.bam lane1_B5_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C5_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D5_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day7_43.bam lane1_B7_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C7_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D7_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day10_43.bam lane1_B10_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C10_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D10_43_2_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Day14_43.bam lane1_B14_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_C14_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_D14_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_Primed_43.bam lane1_P2_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_P4_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_NE_43.bam lane1_NE1_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE2_43_2_merged_L001_GRCh38_bowtie2.sorted.bam lane1_NE3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam
samtools merge -@ 30 merged_DE_43.bam lane1_DE1_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE2_43_1_merged_L001_GRCh38_bowtie2.sorted.bam lane1_DE3_43_1_merged_L001_GRCh38_bowtie2.sorted.bam


echo "All the merging are done"


### Index the merged files
for f in merged_*.bam
do
  echo "Indexing $f"
  samtools index "$f"
done

echo "The merged files are indexed"
