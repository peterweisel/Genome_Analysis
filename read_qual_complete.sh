#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J read_trim
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load FastQC
module load trimmomatic

# input and output directories
RAW=/domus/h1/pewe9965/Genome_Analysis/raw
QC=/home/pewe9965/Genome_Analysis/QC
mkdir -p ${QC}

declare -a lst=("chr3_illumina_R1" "chr3_illumina_R2")
for SAMPLE in ${lst[@]};
do
    OUT=${QC}/${SAMPLE}
    mkdir -p ${OUT}
    # run FASTQC for Illumina short reads
    fastqc -o ${OUT} ${RAW}/${SAMPLE}.fastq.gz
done

# trim Illumina short reads using adapter sequences
trimmomatic PE -phred33 \
${RAW}/chr3_illumina_R1.fastq.gz ${RAW}/chr3_illumina_R2.fastq.gz \
${QC}/chr3_illumina_R1_paired_trimmed.fastq.gz ${QC}/chr3_illumina_R1_unpaired_trimmed.fastq.gz \
${QC}/chr3_illumina_R2_paired_trimmed.fastq.gz ${QC}/chr3_illumina_R2_unpaired_trimmed.fastq.gz \
ILLUMINACLIP:/sw/bioinfo/trimmomatic/0.39/rackham/adapters/TruSeq3-PE.fa:2:30:10 \
MINLEN:15
