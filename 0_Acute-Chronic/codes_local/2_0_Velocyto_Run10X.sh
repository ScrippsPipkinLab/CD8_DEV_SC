#!/bin/bash
# Velocyto installed on Python 3.8 conda environment [velocyto] (PipkinPrecisionTower)
#----------------------------------------
# Input directory organization:
#$ tree '/media/pipkin/ROCKET-PRO/Exp391'
#/media/pipkin/ROCKET-PRO/Exp391
#└── outs
#    ├── analysis.tar.gz
#    ├── feature_reference.csv
#    ├── filtered_feature_bc_matrix
#    │   ├── barcodes.tsv
#    │   ├── barcodes.tsv.gz
#    │   ├── features.tsv
#    │   ├── features.tsv.gz
#    │   ├── matrix.mtx
#    │   └── matrix.mtx.gz
#    ├── metrics_summary.csv
#    ├── possorted_genome_bam.bam
#    ├── possorted_genome_bam.bam.bai
#    └── web_summary.html
#----------------------------------------

wkdir='/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_2_Velocyto'
cd $wkdir

input_path='/media/pipkin/ROCKET-PRO/Exp391'
gtf_path='/home/pipkin/references/GRCm38/genes/genes.gtf'

velocyto run10x $input_path $gtf_path &> Velocyto_run10x.oe