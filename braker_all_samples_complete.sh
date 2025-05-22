#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 04:00:00
#SBATCH -J braker_all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load HISAT2/2.2.1
module load samtools/1.19
module load braker/2.1.6
module load GeneMark/4.72-es

# allow access to species directory and copy gm_key to home directory
chmod a+w -R /home/pewe9965/Genome_Analysis/src/augustus_config/species/
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key
source $AUGUSTUS_CONFIG_COPY

# paths for AUGUSTUS
export AUGUSTUS_CONFIG_PATH=/home/pewe9965/Genome_Analysis/src/augustus_config
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.72-es/snowy

# input and output directories
DIR_DATA=/home/pewe9965/Genome_Analysis/genome/annotation/repeat_mask
DIR_ANN=/home/pewe9965/Genome_Analysis/raw/embryophyte_proteomes.faa
DIR_OUT=/home/pewe9965/Genome_Analysis/genome/annotation/braker_all
mkdir -p ${DIR_OUT}
cd ${DIR_OUT}

# RNA-seq aligned BAM files
BAM_LIST="/domus/h1/pewe9965/Genome_Analysis/RNAseq/Control_1/Control_1_sorted.bam,\
/domus/h1/pewe9965/Genome_Analysis/RNAseq/Control_2/Control_2_sorted.bam,\
/domus/h1/pewe9965/Genome_Analysis/RNAseq/Control_3/Control_3_sorted.bam,\
/domus/h1/pewe9965/Genome_Analysis/RNAseq/Heat_treated_42_12h_1/Heat_treated_42_12h_1_sorted.bam,\
/domus/h1/pewe9965/Genome_Analysis/RNAseq/Heat_treated_42_12h_2/Heat_treated_42_12h_2_sorted.bam,\
/domus/h1/pewe9965/Genome_Analysis/RNAseq/Heat_treated_42_12h_3/Heat_treated_42_12h_3_sorted.bam"

# run BRAKER using soft-masked genome, protein sequences, and aligned RNA-seq reads
/sw/bioinfo/braker/2.1.6/rackham/scripts/braker.pl \
  --species=n.japonicum \
  --useexisting \
  --genome=${DIR_DATA}/pilon.fasta.masked \
  --prot_seq=${DIR_ANN} \
  --bam=${BAM_LIST} \
  --gff3 \
  --etpmode \
  --softmasking \
  --cores 16 \
  --AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_SCRIPTS_PATH \
  --AUGUSTUS_BIN_PATH=$AUGUSTUS_BIN_PATH \
  --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH

