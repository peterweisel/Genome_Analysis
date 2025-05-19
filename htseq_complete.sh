#!/bin/bash -l
#SBATCH -A uppmax2025-3-3 -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 07:00:00
#SBATCH -J htseq_count
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load htseq
#module load Subread

# input and output directories
DIR_DATA=/domus/h1/pewe9965/Genome_Analysis/RNAseq
DIR_ANNOT=/home/pewe9965/Genome_Analysis/genome/annotation/braker_all/braker
DIR_OUT=/home/pewe9965/Genome_Analysis/RNAseq/new_counts        
mkdir -p ${DIR_OUT}
cd ${DIR_OUT}

# run htseq for all samples
declare -a lst=("Control_1" "Control_2" "Control_3" "Heat_treated_42_12h_1" "Heat_treated_42_12h_2" "Heat_treated_42_12h_3")
for NAME in ${lst[@]};
do
    # count transcripts using sorted BAM file and BRAKER output
    htseq-count -f bam -r pos -t gene -s no -t exon -i gene_id ${DIR_DATA}/${NAME}/${NAME}_sorted.bam ${DIR_ANNOT}/braker.gff3 > ${DIR_OUT}/${NAME}_counts.txt
done
