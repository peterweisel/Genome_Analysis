#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 07:00:00
#SBATCH -J flye_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load Flye

# input and output directories
RAW=/proj/uppmax2025-3-3/Genome_Analysis/4_Zhou_2023/reads
SAMPLE=chr3_clean_nanopore
OUT=/home/pewe9965/Genome_Analysis/genome
mkdir -p ${OUT}

# run flye using nanopore long reads and estimated genome size
flye --nano-raw ${RAW}/${SAMPLE}.fq.gz --genome-size 184m --out-dir ${OUT} --threads 16

