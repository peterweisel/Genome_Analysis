#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00:01:00
#SBATCH -J genome_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load HISAT2

# input and output directories
WORK=/domus/h1/pewe9965
DIR_REF=/home/pewe9965/Genome_Analysis/genome/polish/pilon_results
HISAT_INDX=${DIR_REF}/hisat_index
mkdir -p ${HISAT_INDX}
cd ${HISAT_INDX}

# build genome index from polished assembly using HISAT2
hisat2-build -p 8 ${DIR_REF}/pilon.fasta ${HISAT_INDX}/genome

