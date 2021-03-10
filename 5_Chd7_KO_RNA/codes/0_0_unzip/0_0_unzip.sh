#!/bin/bash
#PBS -l mem=4gb
#PBS -t 1-24

in_dir=/gpfs/home/hdiao/Exp276

new_dir=$in_dir/276_${PBS_ARRAYID}

cd $new_dir

for file_x in *.gz
do gunzip $file_x
done

