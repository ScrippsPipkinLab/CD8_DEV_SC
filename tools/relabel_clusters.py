#!/usr/bin/env python

import os
import pandas as pd
from pathlib import Path
import sys
from shutil import copyfile

print("--------------------------------------------------")
print("##### Relabel files path") 
print("      and csv column & row index")
print("      based on label conversion reference")
print("# Within file names, Labels are parsed by _ as delimiter")
print("# Usage: ")
print("#       $python relabel_clusters.py label_file input_dir")

label_file = sys.argv[1]
input_dir = sys.argv[2]

label_df = pd.read_csv(label_file)

label_dict = {str(x):str(y) for index, (x,y) in 
              enumerate(zip(label_df['old_label'].tolist(), 
                            label_df['new_label'].tolist()))}
print(label_dict)
print(input_dir)

for path in Path(input_dir).rglob('*'):
    if path.is_file():
        # Relabel the path containing the file
        relative_path = str(path.resolve()).replace(input_dir, "")
        relative_path_list = relative_path.split("/")
        new_relative_path_list = []
        for name_i in relative_path_list:
            name_i_list = name_i.split("_")
            name_i_list = [label_dict[x] if x in label_dict.keys() 
                              else x for x in name_i_list]
            new_relative_path_list.append("_".join(name_i_list)) 
        new_relative_path = "/".join(new_relative_path_list)

        if new_relative_path != relative_path:
            new_folder = input_dir+"/".join(new_relative_path_list[:-1])
            Path(new_folder).mkdir(parents=True, exist_ok=True)
        else:
            new_relative_path = relative_path
        new_relative_path = new_relative_path.replace('.csv', '_relabeled.csv')
        print(input_dir + new_relative_path)
            
        # Relabel column and index if it is csv file
        if ".csv" in new_relative_path:
            df = pd.read_csv(path)

            df_cols = []
            for i in df.columns:
                i_list = i.split("_")
                i_list = [label_dict[str(x)] if str(x) in label_dict.keys() else str(x) for x in i_list]
                df_cols.append("_".join(i_list))
            df_idx = []
            for i in df.iloc[:,0]:
                i_list = i.split("_")
                i_list = [label_dict[str(x)] if str(x) in label_dict.keys() else str(x) for x in i_list]
                df_idx.append("_".join(i_list))
            #print(df.columns)
            #print(df_cols)
            df.columns = df_cols
            df.iloc[:,0] = df_idx
            df.to_csv(input_dir + new_relative_path, index=False)
        # Otherwise just copy it to a new name
        else:
            copyfile(path, input_dir+new_relative_path)

print("Done!")