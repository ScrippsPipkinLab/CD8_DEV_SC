#!/bin/bash

input_dir='/media/pipkin/ROCKET-PRO/CD8_DEV_SC/4_Hi-C/0_validPairs_XL_folder'

cd $input_dir

source activate ame

for i in *.validPairs.txt
do
  hicpro2juicebox.sh -i $i -g mm10 -j /home/pipkin/Pkgs/juicer_tools_1.22.01.jar
done
