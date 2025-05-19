#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J agat
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries 
module load bioinfo-tools
module load samtools
module load AGAT/1.3.2
module load bamtools 

# input and output directories/files
GFF3_FILE=/domus/h1/pewe9965/Genome_Analysis/genome/annotation/braker_all/braker/braker.gff3
JOB_DIR=/home/pewe9965/Genome_Analysis/genome/agat
AGAT_CONFIG=/sw/bioinfo/AGAT/1.3.2/rackham/lib/site_perl/5.32.1/auto/share/dist/AGAT/agat_config.yaml
REF_GENOME=/home/pewe9965/Genome_Analysis/genome/annotation/repeat_mask/pilon.fasta.masked
mkdir -p ${JOB_DIR}
cd ${JOB_DIR}

# clean GFF3 file
#agat_convert_sp_gxf2gxf.pl --gff ${GFF3_FILE} -c ${AGAT_CONFIG} -o braker_cleaned.gff3

# standardize
#agat_sp_manage_IDs.pl -gff braker_cleaned.gff3 -o braker_standardized.gff3 --prefix NJAP --ensembl

# get feature statistics 
agat_sp_statistics.pl --gff braker_standardized.gff3 -o braker_standardized_statistics.txt

# get protein sequences
agat_sp_extract_sequences.pl --gff braker_standardized.gff3 --fasta ${REF_GENOME} -p -o braker_standardized.aa


