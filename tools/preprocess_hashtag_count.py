#!/usr/bin/env python3


import numpy as np
import pandas as pd
import scanpy as sc
import os
import csv

### Self defined functions for identifying cell antibody hashtags
def findAb(ftTsv):
    AbDict = {}
    with open(ftTsv, "r") as fin:
        rfin = csv.reader(fin, delimiter="\t")
        rowN = 0
        for row in rfin:
            rowN += 1
            if row[2] == "Antibody Capture":
                AbDict[rowN] = row[0]
    return(AbDict)

def cellAbReads(root_name, ftTsv, mtxMtx):
    outName = "%s_Cells_hashTags_count.csv"%root_name
    abDict = findAb(ftTsv)
    with open(mtxMtx, "r") as fin:
        with open(outName, "w") as fout:
            rfin = csv.reader(fin, delimiter=" ")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(["Cell"] + list(abDict.values()))
            
            next(rfin)
            next(rfin)
            next(rfin)
            step = 1
            outRow = [1] + [0 for x in abDict]
            for row in rfin:
                cellN = int(row[1])
                if cellN != step:
                    wfout.writerow(outRow)
                    step = cellN
                    outRow = [row[1]] + [0 for x in abDict]
                if int(row[0]) in list(abDict.keys()):
                    abIdx = list(abDict.keys()).index(int(row[0])) + 1
                    outRow[abIdx] = int(row[2])
            wfout.writerow(outRow)

def cellType(root_name, count_cutoff, amb_cutoff):
    cellTag = "%s_Cells_hashTags_count.csv"%root_name
    outfile = "%s_Cells_hashTags.csv"%root_name
    
    with open(cellTag, "r") as fin:
        with open(outfile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            all_types = next(rfin)[1:]
            cell_types = []
            for row in rfin:
                row_nu = [int(x) for x in row[1:]]
                rowMax = max(row_nu)
                rowSum = sum(row_nu)
                row_type = all_types[row_nu.index(rowMax)]
                if rowMax < float(rowSum)*amb_cutoff:
                    #print('Error in cell %s, maximum hashtag reads less than 70 percent...' %row[0])
                    #print(row)
                    row_type = "Doublet"
                if rowSum < count_cutoff:
                    row_type = "Negative"

                    
                wfout.writerow([row_type])
                cell_types.append(row_type)
    types_set = list(set(cell_types))
    
    
    sum_df = pd.DataFrame({"sample":[root_name], "total_cell_number":[len(cell_types)],
               "negative_percentage":[round(float(cell_types.count("Negative"))*100/len(cell_types), 2)],
               "doublet_percentage":[round(float(cell_types.count("Doublet"))*100/len(cell_types), 2)]})
    for i in types_set:
        if i != "Negative" and i != "Doublet":
            sum_df[i] = [round(float(cell_types.count(i))*100/len(cell_types), 2)]
    
    sum_df.to_csv("%s_hashtagSummary.csv"%root_name)
    print(sum_df)
    
    return(sum_df)