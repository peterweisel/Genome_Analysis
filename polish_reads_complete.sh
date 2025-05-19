#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 24:00:00
#SBATCH -J polish_reads
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load Pilon
module load samtools
module load bwa

# input and output directories
DIR_DATA=/home/pewe9965/Genome_Analysis/genome
RAW=/domus/h1/pewe9965/Genome_Analysis/raw
# de novo genome FASTA file
SAMPLE=${DIR_DATA}/assembly.fasta
SHORT_READ=chr3_illumina
OUT=${DIR_DATA}/polish
mkdir -p ${OUT}
cd ${OUT}

# index de novo genome
bwa index ${SAMPLE}

# map short reads using bwa mem (paired-end reads longer than 70bp)
bwa mem -t 2 ${SAMPLE} ${RAW}/${SHORT_READ}_R1.fastq.gz ${RAW}/${SHORT_READ}_R2.fastq.gz | samtools view -bS -o bwa_mapping.bam

# sort and index BAM file
samtools sort bwa_mapping.bam -o bwa_mapping_sorted.bam
samtools index bwa_mapping_sorted.bam

# run Pilon
mkdir -p ${OUT}/pilon_results
java -jar $PILON_HOME/pilon.jar --genome ${SAMPLE} --bam /home/pewe9965/Genome_Analysis/genome/polish/bwa_mapping_sorted.bam --outdir ${OUT}/pilon_results --threads 6
