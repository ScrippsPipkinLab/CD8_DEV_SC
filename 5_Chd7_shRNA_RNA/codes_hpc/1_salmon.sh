#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --mem=40gb

module load salmon

wkdir=/gpfs/group/pipkin/hdiao/Exp124/fastq
cd $wkdir

salmon_index=/gpfs/group/pipkin/hdiao/ref_resources/mm/release100/GRCm38.salmon.index

gunzip 124_1_R1.fastq.gz &
gunzip 124_2_R1.fastq.gz &
gunzip 124_3_R1.fastq.gz &
gunzip 124_4_R1.fastq.gz &
gunzip 124_5_R1.fastq.gz &
gunzip 124_6_R1.fastq.gz &
gunzip 124_7_R1.fastq.gz &
gunzip 124_8_R1.fastq.gz &
gunzip 124_9_R1.fastq.gz
wait
gunzip 124_1_R2.fastq.gz &
gunzip 124_2_R2.fastq.gz &
gunzip 124_3_R2.fastq.gz &
gunzip 124_4_R2.fastq.gz &
gunzip 124_5_R2.fastq.gz &
gunzip 124_6_R2.fastq.gz &
gunzip 124_7_R2.fastq.gz &
gunzip 124_8_R2.fastq.gz &
gunzip 124_9_R2.fastq.gz

#################### Pair-end mapping ####################
salmon quant -i $salmon_index -l A -1 124_1_R1.fastq -2 124_1_R2.fastq --validateMappings -o 124_1 &
salmon quant -i $salmon_index -l A -1 124_2_R1.fastq -2 124_2_R2.fastq --validateMappings -o 124_2 &
salmon quant -i $salmon_index -l A -1 124_3_R1.fastq -2 124_3_R2.fastq --validateMappings -o 124_3 &
salmon quant -i $salmon_index -l A -1 124_4_R1.fastq -2 124_4_R2.fastq --validateMappings -o 124_4
wait
salmon quant -i $salmon_index -l A -1 124_5_R1.fastq -2 124_5_R2.fastq --validateMappings -o 124_5 &
salmon quant -i $salmon_index -l A -1 124_6_R1.fastq -2 124_6_R2.fastq --validateMappings -o 124_6 &
salmon quant -i $salmon_index -l A -1 124_7_R1.fastq -2 124_7_R2.fastq --validateMappings -o 124_7 &
salmon quant -i $salmon_index -l A -1 124_8_R1.fastq -2 124_8_R2.fastq --validateMappings -o 124_8
wait
salmon quant -i $salmon_index -l A -1 124_9_R1.fastq -2 124_9_R2.fastq --validateMappings -o 124_9

