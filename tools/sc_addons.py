import os
import sys
import csv
import pandas as pd
import scanpy as sc
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import seaborn as sns
sc.settings.verbosity = 1

###----- Perform 3 tests & create venn diagram of overlapping sig DE genes
def multitest_venn(adata, cps, obs_useCol, de_wkdir):
    Path(de_wkdir).mkdir(parents=True, exist_ok=True)
    
    if "sparse" in str(type(adata.raw.X)):
        adata_raw_expr = pd.DataFrame.sparse.from_spmatrix(adata.raw.X)
    else:
        adata_raw_expr = pd.DataFrame(adata.raw.X)
    adata_raw_expr.columns = adata.raw.var_names
    all_genes = adata.raw.var_names.tolist()
    gene_n = adata.raw.n_vars
    
    labels_uniq = list(set(adata.obs[obs_useCol]))
    labels_uniq.sort()
    print(labels_uniq)
    
    for l_t_1 in labels_uniq:
        # Make a new folder for each cluster
        Path("%s/%s"%(de_wkdir, l_t_1)).mkdir(parents=True, exist_ok=True)
        os.chdir("%s/%s"%(de_wkdir, l_t_1))

        for l_t_2 in labels_uniq:
            if l_t_2 != l_t_1:


                # Test result summary dataframe with averge gene expression per condition
                adata_raw_expr_use_df = adata_raw_expr.copy()
                adata_raw_expr_use_df['cat'] = adata.obs[obs_useCol].tolist()
                avg_expr_df = adata_raw_expr_use_df.groupby('cat').mean().T
                avg_expr_df = avg_expr_df.reset_index()
                avg_expr_df.columns = ['gene_names'] + avg_expr_df.columns.tolist()[1:]
                genes_df_l_t = avg_expr_df

                dict_l_t = {}
                for cp in cps:
                    print(l_t_1, l_t_2,cp)
                    # Run test
                    sc.tl.rank_genes_groups(adata, obs_useCol, groups=[l_t_1], reference=l_t_2, method=cp, n_genes=gene_n)

                    # Build dataframe with test results
                    cp_names = [x[0] for x in adata.uns['rank_genes_groups']['names']]
                    cp_fc = [x[0] for x in adata.uns['rank_genes_groups']['logfoldchanges']]
                    cp_score = [x[0] for x in adata.uns['rank_genes_groups']['scores']]
                    cp_padj = [x[0] for x in adata.uns['rank_genes_groups']['pvals_adj']]
                    cp_df = pd.DataFrame({"gene_names":cp_names, "%s_logfc"%cp: cp_fc, "%s_padj"%cp:cp_padj, "%s_score"%cp:cp_score})

                    # Merge with total dataframe
                    genes_df_l_t = pd.merge(genes_df_l_t, cp_df, how='left', on="gene_names")
                    cp_sig_df = cp_df.copy()
                    cp_sig_df = cp_sig_df[cp_sig_df['%s_padj'%cp]<= 0.05]
                    cp_sig_df = cp_sig_df[cp_sig_df['%s_logfc'%cp]>0]
                    dict_l_t[cp] = cp_sig_df['gene_names'].tolist()

                # Save results for all tests
                genes_df_l_t.to_csv("%s_vs_%s_differential.csv"%(l_t_1, l_t_2), index=False)
                
                dict_l_t_df_dict = {}
                for cp in cps:
                    dict_l_t_df_dict[cp] = pd.Series(dict_l_t[cp])
                dict_l_t_df = pd.DataFrame(dict_l_t_df_dict)
                dict_l_t_df.to_csv("%s_vs_%s_venn.csv"%(l_t_1, l_t_2), index=False)

                # Test comparison venn diagram
                try:
                    venn3([set(dict_l_t[cps[0]]), set(dict_l_t[cps[1]]), set(dict_l_t[cps[2]])], tuple(cps))
                        
                    outname = "%s_vs_%s_venn.png"%(l_t_1, l_t_2)
                    plt.savefig(outname)
                    plt.show()
                    plt.close()
                        
                except:
                    print(l_t_1, l_t_2, "Venn3 plotting error")


                
###----- Calculate average UMAP connectivities between groups
def umap_conn_avg(adata, obs_usekey):
    conn_df = pd.DataFrame(adata.obsp['connectivities'].todense())
    
    obs_uniq = np.unique(adata.obs[obs_usekey])
    obs_uniq.sort()
    avg_conn_df = pd.DataFrame()

    for i in obs_uniq:
        i_avg_conn = []
        for j in obs_uniq:
            i_idx = [True if x == i else False for x in adata.obs[obs_usekey] ]
            j_idx = [True if x == j else False for x in adata.obs[obs_usekey] ]

            conn_df_ij = conn_df.copy()
            conn_df_ij = conn_df[i_idx].iloc[:, j_idx]
            i_avg_conn.append(np.mean(conn_df_ij.values))
        avg_conn_df[i] = i_avg_conn
    avg_conn_df.index = obs_uniq
    
    avg_conn_df.to_csv("%s_UMAP_conn.csv"%obs_usekey)
    
###----- Select signature genes of group based on DE results of 
#        group v.s. most similar groups
def sig_to_neighbors(adata, use_key, conn, diff_dir,
                     neighbor_n,padj_c = 0.1, top_n=300, padj_list = ['t-test_padj', 't-test_overestim_var_padj']):
    
    # Plot connectivity heatmap
    conn_df = pd.read_csv(conn)
    conn_df.columns = [use_key] + conn_df.columns.tolist()[1:]
    conn_df = conn_df.set_index(use_key)

    for i in range(0, len(conn_df.columns)):
        i_colname = conn_df.columns[i]
        conn_df[i_colname] = conn_df[i_colname] / conn_df[i_colname].values[i]

    for i in range(0, len(conn_df.columns)):
        conn_df.iloc[i,i] = 0
    
    #--- Plot
    # Sort columns and indexes by numeric order
    newcols = conn_df.columns.tolist()
    newcols.sort(key=int)
    newindex = conn_df.index.tolist()
    newindex.sort(key=int)
    conn_df = conn_df[newcols]
    conn_df = conn_df.loc[newindex]
    
    plt.figure(figsize=(8, 7))
    sns.heatmap(conn_df, cmap='viridis')
    plt.savefig(conn.replace(".csv", ".png"))
    plt.show()
    
    
    sig_dict = {}
    # Compare with top similar groups and select signature genes
    for i in conn_df.columns:
        i_df = pd.DataFrame(conn_df[i]).sort_values(i, ascending=False)
        i_neighbors = i_df.index.to_list()[:neighbor_n]
        i_sig = []

        for i_n in i_neighbors:
            i_vs_i_n = diff_dir + '/' + str(i) + '/' + str(i) + '_vs_' + str(i_n) + '_differential.csv'
            i_vs_i_n_df = pd.read_csv(i_vs_i_n)
            i_vs_i_n_df = i_vs_i_n_df[i_vs_i_n_df['t-test_logfc'] > 0]
            
            for padj in padj_list:
                i_vs_i_n_df = i_vs_i_n_df[i_vs_i_n_df[padj] <= padj_c]
                
            i_vs_i_n_df = i_vs_i_n_df.sort_values('t-test_logfc', ascending=False)
            top_sig = i_vs_i_n_df['gene_names'].values.tolist()
            if len(top_sig) > top_n:
                top_sig = top_sig[:top_n]
            i_sig.append(top_sig)
        i_sig = set.intersection(*map(set, i_sig))
        print(use_key, i,'signature gene number: ', len(i_sig))
        
        sig_dict[i] = pd.Series(list(i_sig))
    
    #--- Write output
    sig_df = pd.DataFrame(sig_dict)
    sig_df.to_csv("%s_sigGenes.csv"%use_key, index=False)
        
    return(sig_dict)
        

    






















