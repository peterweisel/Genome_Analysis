#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 03:00:00
#SBATCH -J chloroplast
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load GetOrganelle
module loaf GetOrganelleDB

# input and output directories
WORK=/domus/h1/pewe9965
DIR_RAW=${WORK}/Genome_Analysis/raw
DIR_OUT=${WORK}/Genome_Analysis/genome/chloroplast
mkdir -p ${DIR_OUT}
cd ${DIR_OUT}

# assemble chloroplast genome using Illumina short reads
get_organelle_from_reads.py \
  -1 ${DIR_RAW}/chr3_illumina_R1.fastq.gz \
  -2 ${DIR_RAW}/chr3_illumina_R2.fastq.gz \
  -o mito_getorganelle_output \
  -R 15 \
  -k 21,45,65,85,105 \
  -F embplant_pt \
  -t 4

