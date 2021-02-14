#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 14:17:32 2019

@author: yolandatiao
"""

########## isoform_norm_count ##########
# Author: Huitian (Yolanda) Diao
# Feb 12th, 2019

# Input files:
# - salmon output quant.sf files
# - Selected Chd7 isoform reads & sum of counted reads <- 2_0_Chd7_slt--ReadSum.sh

########## Import ##########
import os
import glob
import csv
from astropy.io import ascii
from astropy.table import Table, join, vstack


########## Main ##########
wk_dir = "/Volumes/Huitian/Exp276_RNASeq/Chd7"
os.chdir(wk_dir)

#####----- Compile sum of read numbers for each sample

sum_dict = {}
for file in glob.glob("/Volumes/Huitian/Exp276_RNASeq/salmon_nonTrim/*/quant.sf.NumReads.sum"):
    with open(file, "r") as fin:
        rfin = csv.reader(fin)
        num = int(float(next(rfin)[0]))
        sp_name = file.split("/")[-2]
        sum_dict[sp_name] = num

#####----- Compile all the chd7 splice variant
id_list = ["ENSMUST00000051558.9", "ENSMUST00000039267.9", 
        "ENSMUST00000170391.1", "ENSMUST00000127476.7", 
        "ENSMUST00000222546.1", "ENSMUST00000129655.1", 
        "ENSMUST00000130709.1", "ENSMUST00000170457.1"]

chd7_tab = Table()
chd7_tab["esblID"] = id_list

for file in glob.glob("/Volumes/Huitian/Exp276_RNASeq/salmon_nonTrim/*/quant.sf.chd7"):
    tpm_list = [None for x in id_list]
    with open(file, "r") as fin:
        rfin = csv.reader(fin, delimiter = "\t")
        for row in rfin:
            row_idx = id_list.index(row[0])
            tpm_list[row_idx] = float(row[3])
    sp_name = file.split("/")[-2]
    chd7_tab[sp_name] = tpm_list

ascii.write(chd7_tab, "chd7_sv_tpm.csv", format="csv", overwrite=True)  


#####----- Chd7 splice variant normalized tpm

id_list = ["ENSMUST00000051558.9", "ENSMUST00000039267.9", 
        "ENSMUST00000170391.1", "ENSMUST00000127476.7", 
        "ENSMUST00000222546.1", "ENSMUST00000129655.1", 
        "ENSMUST00000130709.1", "ENSMUST00000170457.1"]

chd7_tab = Table()
chd7_tab["esblID"] = id_list

for file in glob.glob("/Volumes/Huitian/Exp276_RNASeq/salmon_nonTrim/*/quant.sf.chd7"):
    tpm_list = [None for x in id_list]
    with open(file, "r") as fin:
        rfin = csv.reader(fin, delimiter = "\t")
        for row in rfin:
            row_idx = id_list.index(row[0])
            tpm_list[row_idx] = float(row[3])
    sp_name = file.split("/")[-2]
    sp_totalRead = sum_dict[sp_name]
    sp_adj = 10000000.0 / sp_totalRead
    tpm_list_adj = [x*sp_adj for x in tpm_list]    
    chd7_tab[sp_name] = tpm_list_adj

ascii.write(chd7_tab, "chd7_sv_tpm_norm.csv", format="csv", overwrite=True)  
