#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J eggnog
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load bioinformatics modules
module load bioinfo-tools
module load eggNOG-mapper/2.1.9

DIR_FASTA=/home/pewe9965/Genome_Analysis/genome/agat
DIR_OUT=/home/pewe9965/Genome_Analysis/genome/annotation/eggnog_geneid
mkdir -p ${DIR_OUT}

echo "emapper started"

emapper.py -i ${DIR_FASTA}/braker_standardized.aa -m diamond -o n_japonicum_eggnog --cpu 8 --itype proteins --output_dir ${DIR_OUT} \
    --decorate_gff ${DIR_FASTA}/braker_standardized.gff3 --decorate_gff_ID_field gene_id --output_dir ${DIR_OUT} --excel --go_evidence experimental --override

echo "emapper complete"
