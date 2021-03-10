#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

in_dir=/gpfs/home/hdiao/Exp276

new_dir=$in_dir/276_${PBS_ARRAYID}

cd $new_dir

new_name1=276_${PBS_ARRAYID}_R1.fastq
new_name2=276_${PBS_ARRAYID}_R2.fastq

cat *R1_001.fastq > $new_name1
cat *R2_001.fastq > $new_name2
