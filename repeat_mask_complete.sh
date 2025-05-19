#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH -J soft_mask
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load RepeatMasker

# input and output directories
DIR_DATA=/home/pewe9965/Genome_Analysis/genome/polish/pilon_results
OUT=/home/pewe9965/Genome_Analysis/genome/annotation/repeat_mask
mkdir -p ${OUT}

# soft-mask genome using polished assembly
RepeatMasker -pa 4 -a -e ncbi -xsmall -dir ${OUT} ${DIR_DATA}/pilon.fasta

