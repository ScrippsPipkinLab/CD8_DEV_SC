#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########---------- Select markers for flow cytometry ----------##########
"""
Created on Fri Feb 21 09:47:23 2020

@author: yolandatiao
"""


##########---------- Imports ----------##########
# General packages
import os
import pandas as pd
import numpy as np
import glob
import pydotplus
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from collections import Counter

# Scikit learn
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.externals.six import StringIO
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score


##########---------- Self defined functions ----------##########
def slt_highlyDiff_Deseq(inFiles, log2fcCutoff=1, pvalCutoff=0.05, padjCutoff=0.05):
    #inFiles = deseq_files
    #log2fcCutoff = 1.5
    #pvalCutoff = 0.01
    #padjCutoff = 0.01
    out_list = []
    for file_i in inFiles:
        tb_i = pd.read_csv(file_i)
        tb_i = tb_i[tb_i["pvalue"] <= pvalCutoff]
        tb_i = tb_i[tb_i["padj"] <= padjCutoff]
        tb_i = tb_i[abs(tb_i["log2FoldChange"]) >= log2fcCutoff]
        out_list += list(tb_i.iloc[:,0])
    out_list = list(set(out_list))
    return(out_list)
        
    


####################---------- Main ----------####################
##########---------- Parameters ----------##########

#--- Cell surface protein gene name lists
markers_file = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/MM_MARKERS.csv"
surfaceome = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/z_Resources/murine_surfaceome.csv"

#--- DEseq file for finding differential genes
deseq_dir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/2_norm_counts-expr_DESeq2/all--numslt-rmWTNAV"
deseq_files = glob.glob("%s/*.csv"%deseq_dir)
deseq_files = [i for i in deseq_files if "vs" in i]
deseq_parameters = {'fcCutoff': 1.5, 'pvalCutoff': 0.01, 'padjCutoff':0.01}
classify_categories = ['1', '2', '3', '4', '5', '7', '8']
deseq_files_use = []
for i in deseq_files:
    i_name = i.split('/')[-1]
    i_category = i_name.replace("Arm_", "").replace(".csv", "").replace("L", "").split("--vs--")
    #if ('8' in i_category) and ('0' not in i_category) :
    #    deseq_files_use.append(i)
    if not any([False if x in classify_categories else True for x in i_category]):
        deseq_files_use.append(i)
deseq_files = deseq_files_use  

#--- Training data
# Norm count file's first column need to be cell barcode
# Obs file's first column need to be cell barcode
norm_count_file = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/2_norm_counts-expr_DESeq2/all--numslt-rmWTNAV/all_norm_counts_named_c10_int_slt.csv"
pctl_z_file = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/2_norm_counts-expr_DESeq2/all--numslt-rmWTNAV/all_norm_counts_named_int_slt_c10_nbPctl_Z_naOmit.csv"
obs_file = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/2_norm_counts-expr_DESeq2/all--numslt-rmWTNAV/all--numSlt-rmWTNAV_obs_slt.csv"


#--- Training target
use_key = "louvain" 
target_category = '3'



####################---------- Read inputs ----------####################
###----- Read input files
obs_tb = pd.read_csv(obs_file)
obs_tb.columns = ["cell_bc"] +  list(obs_tb.columns)[1:]
norm_count_tb = pd.read_csv(norm_count_file)
norm_count_tb.columns = ["cell_bc"] +  list(norm_count_tb.columns)[1:]
#pctl_z_tb = pd.read_csv(pctl_z_file)
#pctl_z_tb['cell_bc'] = list(pctl_z_tb['cell_id'])
#del pctl_z_tb['cell_id']

###----- Select highly differential marker genes
markers_tb = pd.read_csv(markers_file)
marker_genes = list(markers_tb["gene_name"])

surfaceome_tb = pd.read_csv(surfaceome)
surfaceome_genes = list(surfaceome_tb["Symbol"])

diff_genes = slt_highlyDiff_Deseq(deseq_files, 
                                  deseq_parameters['fcCutoff'], 
                                  deseq_parameters['pvalCutoff'], 
                                  deseq_parameters['padjCutoff'])
diff_surf_genes = list(set(diff_genes) & set(surfaceome_genes))
diff_marker_genes = list(set(diff_genes) & set(marker_genes))

diff_surf_genes = diff_marker_genes
diff_surf_genes += ["Prdm1"]

#--- Select gene names to use
diff_use_genes = list(set(diff_surf_genes) & set(list(norm_count_tb.columns)))
#diff_use_genes = list(set(diff_use_genes) & set(list(pctl_z_tb.columns)))
use_cols = ["cell_bc"] + diff_use_genes

norm_count_tb_use = norm_count_tb[use_cols]
#pctl_z_tb_use = pctl_z_tb[use_cols]

#--- Merge data tables with obs table
norm_count_tb_use = norm_count_tb_use.set_index("cell_bc")
#pctl_z_tb_use = pctl_z_tb_use.set_index("cell_bc")
obs_tb = obs_tb.set_index("cell_bc")

norm_count_tb_use = norm_count_tb_use.join(obs_tb, on="cell_bc", how="left")
#pctl_z_tb_use = pctl_z_tb_use.join(obs_tb, on="cell_bc", how="left")

#--- Connvert label to "y" and "n"
norm_count_tb_use[use_key] = [str(i) for i in list(norm_count_tb_use[use_key])]
norm_count_tb_use = norm_count_tb_use[[x in classify_categories for x in list(norm_count_tb_use[use_key])]]
norm_count_tb_use[use_key] = ['y' if i == target_category else 'n' for i in list(norm_count_tb_use[use_key]) ]

#pctl_z_tb_use[use_key] = [str(i) for i in list(pctl_z_tb_use[use_key])]
#pctl_z_tb_use = pctl_z_tb_use[[x in classify_categories for x in list(pctl_z_tb_use[use_key])]]
#pctl_z_tb_use[use_key] = ['y' if i == target_category else 'n' for i in list(pctl_z_tb_use[use_key]) ]


####################---------- Training & Tests ----------####################
wkdir = "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/3_norm_counts-expr_DESeq2-sltMarkers/L3/all--numslt-rmWTNAV_markers-Blimp1"
os.chdir(wkdir)   


#####----- Parameters to learn
use_tb = "normCount"
depth_steps = list(range(1,21))
precision_tb = pd.DataFrame({"Depth": depth_steps})
recall_tb = pd.DataFrame({"Depth": depth_steps})
f1_tb = pd.DataFrame({"Depth": depth_steps})

#####----- Run conditions
for pos_weight in range(1, 11):
    pos_weight_precision = []
    pos_weight_recall = []
    pos_weight_f1 = []
    for max_depth in depth_steps:
        out_name_base = use_tb + "__PosWeight-" + str(pos_weight) + "__Depth-" + str(max_depth)
        
        #####----- Training
        category_weight = {'y': pos_weight, 'n': 1}
        
        if use_tb == "normCount":
            data_tb_use = norm_count_tb_use
        #elif use_tb == "pctlZ":
        #    data_tb_use = pctl_z_tb_use
        
        count_sum = Counter(list(data_tb_use[use_key]))
        min_count = min(list(count_sum.values()))
        
        #--- Select dependent and independent
        X = data_tb_use[diff_use_genes].values
        y = data_tb_use[use_key]
        #--- Train test split
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=2)
        #--- Build model & fitting
        markerTree = DecisionTreeClassifier(max_depth = max_depth, 
                                            min_samples_leaf = int(min_count/3),
                                            class_weight = category_weight,
                                            criterion = "entropy",
                                            random_state=0).fit(X_train, y_train)
        #####----- Testing
        yhat = markerTree.predict(X_test)
        
        tree_precision = precision_score(y_test, yhat, labels=['y', 'n'], pos_label='y')
        tree_recall = recall_score(y_test, yhat, labels=['y', 'n'], pos_label='y')
        tree_f1 = f1_score(y_test, yhat, labels=['y', 'n'], pos_label='y')
        
        dot_data = StringIO()
        filename = out_name_base + ".png"
        featureNames = diff_use_genes
        targetNames = y.unique().tolist()
        out=tree.export_graphviz(markerTree,feature_names=featureNames, out_file=dot_data, class_names= np.unique(y_train), filled=True,  special_characters=True,rotate=False)  
        graph = pydotplus.graph_from_dot_data(dot_data.getvalue())  
        graph.write_png(filename)
        
        pos_weight_precision.append(tree_precision)
        pos_weight_recall.append(tree_recall)
        pos_weight_f1.append(tree_f1)
    precision_tb[str(pos_weight)] = pos_weight_precision
    recall_tb[str(pos_weight)] = pos_weight_recall
    f1_tb[str(pos_weight)] = pos_weight_f1        
        
precision_tb_name = out_name_base + "_precision.csv"
recall_tb_name = out_name_base + "_recall.csv"
f1_tb_name = out_name_base + "_f1.csv"
precision_tb.to_csv(precision_tb_name, index=False)
recall_tb.to_csv(recall_tb_name, index=False)
f1_tb.to_csv(f1_tb_name, index=False)











