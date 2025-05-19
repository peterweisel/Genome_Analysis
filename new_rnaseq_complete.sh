#!/bin/bash -l
#SBATCH -A uppmax2025-2-288 -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 00:45:00
#SBATCH -J rnaseq
#SBATCH --mail-type=ALL
#SBATCH --mail-user peter.weisel.9965@student.uu.se
#SBATCH --output=%x.%j.out

# load libraries
module load bioinfo-tools
module load HISAT2
module load deepTools
module load samtools

# input directories
WORK=/domus/h1/pewe9965
DIR_REF=/home/pewe9965/Genome_Analysis/genome/polish/pilon_results
DIR_RAW=${WORK}/Genome_Analysis/raw
HISAT_INDEX=/home/pewe9965/Genome_Analysis/genome/polish/pilon_results/hisat_index/genome

# align all samples
declare -a lst=("Control_1" "Control_2" "Control_3" "Heat_treated_42_12h_1" "Heat_treated_42_12h_2" "Heat_treated_42_12h_3")
for NAME in ${lst[@]};
do
    OUTPUT=/domus/h1/pewe9965/Genome_Analysis/RNAseq/${NAME}
    mkdir -p ${OUTPUT}
    cd ${OUTPUT}

    # align paired-end RNA-seq reads using HISAT2 index files
    hisat2 --phred33 \
        --dta \
        -x ${HISAT_INDEX} \
        -1 ${DIR_RAW}/${NAME}_f1.fq.gz \
        -2 ${DIR_RAW}/${NAME}_r2.fq.gz \
        --rna-strandness RF \
        -p 16 \
    | samtools view -b - \
    | samtools sort -@ 8 -o ${NAME}_sorted.bam

    # output sorted BAM file and indexed BAM file
    samtools index ${NAME}_sorted.bam
done

