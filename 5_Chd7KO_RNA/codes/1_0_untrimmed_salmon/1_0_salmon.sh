#!/bin/bash
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=16
#PBS -t 1-24

in_dir=/gpfs/home/hdiao/Exp276
new_dir=$in_dir/276_${PBS_ARRAYID}
cd $new_dir

module load salmon

file1=276_${PBS_ARRAYID}_R1.fastq
file2=276_${PBS_ARRAYID}_R2.fastq
out_file=/gpfs/home/hdiao/Exp276/276_${PBS_ARRAYID}

salmon quant -i /gpfs/home/hdiao/resources/Transcript/Ensembl_5_17_mm10_qusai_salmon_index  -l A -p 16 -1 $file1 -2 $file2 -o $out_file
