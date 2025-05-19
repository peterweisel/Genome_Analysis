#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J assembly_evaluation
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load BUSCO
module load quast
module load augustus/3.5.0-20231223-33fc04d

# export Augustus PATH
export AUGUSTUS_CONFIG_PATH=/home/pewe9965/Genome_Analysis/src/augustus_config_copy
source $AUGUSTUS_CONFIG_COPY

# input and output directories
DIR_DATA=/home/pewe9965/Genome_Analysis/genome
DIR_OUT=${DIR_DATA}/assembly_eval
mkdir -p ${DIR_OUT}
cd ${DIR_OUT}

# run BUSCO using polished assembly and compare with viridiplantae_odb10 dataset
busco -i ${DIR_DATA}/polish/pilon_results/pilon.fasta -o busco_chr3 -m genome -c 2 -l /domus/h1/pewe9965/Genome_Analysis/src/BUSCO_LINEAGE_SETS/busco_downloads/lineages/viridiplantae_odb10 --augustus -f

# run QUAST using polished assembly
quast.py ${DIR_DATA}/polish/pilon_results/pilon.fasta -o quast_chr3 -t 2 --eukaryote --large
