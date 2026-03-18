#PBS -S /bin/bash
#PBS -N fastp
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=80:mem=80gb
#PBS -o /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/fastp
#PBS -e /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/fastp



outputDir="/rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/fastp/"


### Run the script on all files :
module load anaconda3/personal
source activate fastp

for fastqFile1 in /rds/general/user/mnocente/home/analysis/Marina/4C/4C_Dlx1_NKX23_activation_Dec2024/raw_data/241206_VH00504_255_AAGF33LM5/*R1_001.fastq.gz ; do 
	echo "${fastqFile1}";
	echo "${fastqFile1##*/}"; 
	fastqFile2="${fastqFile1/R1_001.fastq.gz/R2_001.fastq.gz}"; 
	echo "${fastqFile2##*/}";
	echo "______________";

	echo "start fastp"

	fastqFile="${fastqFile1/R1_001.fastq.gz/}"

	echo "${fastqFile##*/}"
	echo "cleaned_filtered_${fastqFile1##*/}"
	echo "cleaned_filtered_${fastqFile2##*/}"
	echo "fastp_${fastqFile##*/}.json"
	echo '______________'

### COMMENTS and EXPLANATIONS :
# For the input files int1 et int2 we need to write the complete path
# outputDir="/rds/general/user/mnocente/home/analysis/Marina/4C/4C_Aug2024/fastp"
# ${fastqFile1##*/} : ##*/ remove all the things before the last /
# sampleFile="$(basename $fastqFile1)"  <==>  "${fastqFile1##*/}"

	### Run fastp :

	fastp --thread 8 \
	--in1 ${fastqFile1} \
	--in2 ${fastqFile2} \
	--out1 ${outputDir}cleaned_filtered_${fastqFile1##*/} \
	--out2 ${outputDir}cleaned_filtered_${fastqFile2##*/} \
	--failed_out ${outputDir}failed_out_${fastqFile##*/}.txt \
	--json ${outputDir}fastp_${fastqFile##*/}.json \
	--html ${outputDir}fastp_${fastqFile##*/}.html \
	--length_required 45 \
	2>${outputDir}fichier_log_${fastqFile##*/}.txt


	echo "done fastp"

done


source deactivate


### Explanations about the used options :
# in1 : name of the input file (read 1)
# in2 : name of the input file (read 2)
# out1 : name of the output file (read 1)
# out2 : name of the output file (read 2)
## unpaired1 : For PE(paired-end) input, if read1 has a succed QC but not read2, it will be written in unpaired1. By default, it's deleted.
## unpaired2 : For PE(paired-end) input, if read2 has a succed QC but not read1, it will be written in unpaired2. By default, it's deleted.
# failed_out : specify the file to store the reads which failed the filters
# json : fastp report in .json format
# html : fastp report in .html format
# thread : 2 by default
# length_required : minimal size of reads: 45 pb (sequencing 60bp)
# the adapters trimming is active by default. 


