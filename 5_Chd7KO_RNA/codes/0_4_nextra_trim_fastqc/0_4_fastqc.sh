#!/bin/bash
#PBS -l mem=8gb
#PBS -t 1-24

in_dir=/gpfs/home/hdiao/Exp276

new_dir=$in_dir/276_${PBS_ARRAYID}

cd $new_dir

module load fastqc

file1=276_${PBS_ARRAYID}_R1_val_1.fq
file2=276_${PBS_ARRAYID}_R2_val_2.fq

fastqc $file1
fastqc $file2
