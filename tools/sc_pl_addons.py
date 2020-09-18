import os
import sys
import csv
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib.patches import Patch
import re
sc.settings.verbosity = 1

def sc_scatter(adata, basis, obs_keys, save, plt_size=[12,4], dot_size=5, alpha=1):
    sns.set_style('white')

    plot_df = pd.DataFrame(adata.obsm['X_' + basis])
    fig, axes = plt.subplots(ncols=len(obs_keys), sharey=True, figsize=(plt_size[0], plt_size[1]), 
                             constrained_layout = True, dpi=240)

    for i in range(0, len(obs_keys)):
        obs_key = obs_keys[i]
        ax = axes[i]
        
        plot_df[obs_key] = adata.obs[obs_key].values.tolist()
        if adata.obs.dtypes[obs_key] not in ['int64', 'float64']:
            plot_df[obs_key] = plot_df[obs_key].astype('category')
        
        plot_i = sns.scatterplot(data=plot_df, x=0, y=1, hue=obs_key,s=dot_size, linewidth=0, alpha=alpha, ax=ax)  
        plot_i.set(xlabel=None, ylabel=None, xticklabels=[], yticklabels=[])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=handles[1:], labels=labels[1:],loc="upper left") # bbox_to_anchor=(1.04,1), 
        if len(labels) > 20:
            ax.legend().set_visible(False)
        ax.set_title(obs_key)
    
    fig.savefig(save)

def count_pctg_stack_bars(adata, key1, key2, count, select_key, select_vals, ax, cols_use, legend=True):
    #key1 = 'condition'
    #key2 = 'leiden'
    #count = 'cell_id'

    obs_df = adata.obs.copy()
    if select_key != None:
        obs_df = obs_df[[True if x in select_vals else False for x in obs_df[select_key].values]]
    obs_df = obs_df.reset_index()[[key1,key2,count]]
    obs_df = obs_df.drop_duplicates(count)
    sum_df = obs_df.groupby([key1, key2]).count()[[count]].unstack(key1)
    sum_df.columns = [x[1] for x in sum_df.columns.values] 
    sum_df[np.isnan(sum_df)] = 0
    
    for i in sum_df.columns:
        sum_df[i] = sum_df[i] / sum(sum_df[i]) * 100
        
    ###----- Order columns by numeric order
    colnames = sum_df.columns.tolist()
    new_colnames = sorted(colnames, key=lambda x: int(re.findall("\d+", x)[0]))
    sum_df = sum_df[new_colnames]
    
    ###----- Plot
    sum_df.T.plot(kind="bar", stacked=True, width=0.4, color=cols_use, ax=ax)
    ax.legend(bbox_to_anchor=(1.04,1)) 
    if not legend:
        ax.legend().set_visible(False)
    if select_key != None:
        ax.set_title(select_key + ": " + select_vals)
    
    return(sum_df)

def count_pctg_bars(adata, key1, key2, count, select_key, select_vals, ax, cols_use, legend=True):
    #key1 = 'condition'
    #key2 = 'leiden'
    #count = 'cell_id'
    #select_key = None
    #select_vals = None

    obs_df = adata.obs.copy()
    if select_key != None:
        obs_df = obs_df[[True if x in select_vals else False for x in obs_df[select_key].values]]
    obs_df = obs_df.reset_index()[[key1,key2,count]]
    obs_df = obs_df.drop_duplicates(count)
    sum_df = obs_df.groupby([key1, key2]).count()[[count]].unstack(key1)
    sum_df.columns = [x[1] for x in sum_df.columns.values] 
    sum_df[np.isnan(sum_df)] = 0

    for i in range(0, len(sum_df)):
        sum_df.iloc[i] = sum_df.iloc[i] / sum(sum_df.iloc[i]) * 100

    sum_df = sum_df.reset_index().melt(id_vars=[key2])
    sum_df.columns = [key2, key1, 'Pctg']
    
    ###----- Plot
    xlabels = np.unique(sum_df[key1]).tolist()
    # If digital use numeric order
    if all([x.isdigit() for x in xlabels]):
        xlabels.sort(key=int)
    else:
        xlabels = sorted(xlabels, key=lambda x: int(re.findall("\d+", x)[0]))
    
    g = sns.barplot(x=key1, y='Pctg', hue=key2, data=sum_df, ax=ax, palette=cols_use, order=xlabels)
    ax.legend(bbox_to_anchor=(1.4,1)) 
    if not legend:
        ax.legend().set_visible(False)
    if select_key != None:
        ax.set_title(select_key + ": " + select_vals)
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    
    return(sum_df)

###----- Donut plots showing distribution of cells in samples
def donut_plot(adata, group_key, sub_key, ax_use, 
               cmaps_use = [plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples, plt.cm.Greys, plt.cm.Oranges]): # can take less than 6 groups
    #group_key = "condition"
    #sub_key = "sample_id_donor"

    group_count_df = adata.obs[[group_key, sub_key]].groupby(group_key).count().reset_index()
    group_names = group_count_df.iloc[:,0].tolist()
    group_size = group_count_df.iloc[:,1].tolist()
    group_cols = [cmaps_use[x](0.6) for x in range(0, len(group_names))]
    subgroup_names = []
    subgroup_size = []
    subgroup_cols = []

    for i in range(0, len(group_names)):
        i_name = group_names[i]
        i_count_df = adata.obs[adata.obs[group_key]==i_name][[sub_key, group_key]].groupby(sub_key).count().reset_index()
        i_count_df.columns = [sub_key, "count"]
        i_count_df = i_count_df[i_count_df["count"]> 0]
        i_count_df = i_count_df.sort_values("count", ascending=False)

        # Calculate cmaps
        i_cmap = cmaps_use[i]
        i_max_c = np.max(i_count_df['count'])
        i_min_c = np.min(i_count_df['count'])
        i_cols = [ i_cmap((x-i_min_c)/(i_max_c-i_min_c)*0.4 + 0.3) for x in i_count_df['count']]

        subgroup_names += i_count_df[sub_key].tolist()
        subgroup_size += i_count_df['count'].tolist()
        subgroup_cols += i_cols

    # First Ring (outside)
    #fig, ax = plt.subplots()
    ax_use.axis('equal')
    mypie, _ = ax_use.pie(group_size, radius=1.3, labels=group_names, colors=group_cols )
    plt.setp( mypie, width=0.3, edgecolor='white')

    # Second Ring (Inside)
    mypie2, _ = ax_use.pie(subgroup_size, radius=1.3-0.3, 
                       labels=None, labeldistance=0.7, colors=subgroup_cols)
    plt.setp( mypie2, width=0.4, edgecolor='white')
    plt.margins(0,0)


###----- Generate legends seperately, used for 9_X_0_CDX_distribution
def leg_generator(ref_df_use, c_key, c_map='hls'):
    if c_key == 'peak_severity_cat':
        key_name_unique = np.array(['none','moderate','sev']) # Use this order
    else:
        key_name_unique = np.unique(ref_df_use[c_key])
    
    # If input array is numeric, use color bar instead of patches
    numerical_input = all([True if str(s).isdigit() else False for s in key_name_unique])
    if numerical_input: # Numeric
        key_name_unique = [float(x) for x in key_name_unique]
        k_min = min(key_name_unique)
        k_max = max(key_name_unique)
        cmap = mpl.cm.get_cmap(c_map)
        norm = mpl.colors.Normalize(vmin=k_min, vmax=k_max)
        
        palette = [cmap((i-k_min)/k_max)[:3] for i in key_name_unique]
        key_name_unique_col_dict = dict(zip(key_name_unique, palette))
        key_cols = ref_df_use[c_key].map(key_name_unique_col_dict)
        
        return(key_cols, [cmap, norm])
        
    else: # Not numeric
        key_name_unique_col_dict = dict(zip(key_name_unique, sns.color_palette(c_map, len(key_name_unique))))
        key_cols = ref_df_use[c_key].map(key_name_unique_col_dict)

        legend_elements = []
        for i in key_name_unique_col_dict.keys():
            legend_elements.append(Patch(facecolor=key_name_unique_col_dict[i], label=i))
            
    return(key_cols, legend_elements)
    
    