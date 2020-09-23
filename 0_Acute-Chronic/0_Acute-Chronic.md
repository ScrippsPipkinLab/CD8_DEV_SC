
# Exp391: Acute v.s. Chronic infection

# Results:

![Louvain labels](https://user-images.githubusercontent.com/26311995/93657818-a519d180-fa03-11ea-8bb9-26373ef7fb19.png)

## 1. Scanpy <br>
   > With filtered resampled cells (capped 1250 cells / category) <br>
   
   #### 1.1 Scanpy & PAGA 
   - [Notebook](codes_local/1_0-0_SCANPY_PAGA_resampled.ipynb)
   - [Data output](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/0_Acute-Chronic_paga/)
   - [Summary output](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/0_sum/)
   - [Per category average expression](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/1_avg_expr/)
   
   #### 1.2 Differential analysis
   - [Notebook a: Cluster v.s. Cluster](0_Acute-Chronic/codes_local/1_0-a_SCANPY_PAGA_resampled_DifferentialAnalysis_cluster-vs-cluster.ipynb)
   - [Notebook b: Armstrong v.s. Cl13](0_Acute-Chronic/codes_local/1_0-a_SCANPY_PAGA_resampled_DifferentialAnalysis_cluster-vs-cluster.ipynb)]
   - Outputs:
     - [Each cluster v.s. all rest](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/eachCluster_vs_All/)
     - [Each cluster v.s. other clusters](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/Cluster_vs_Cluster/)
     - [Armstrong v.s. Clone13 v.s. Naive](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/Arm_vs_Cl13/)
     - [Each type v.s. other types (Naive, A5T, A5P, etc.)](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/celltype_vs_celltype/)
     - [Per cluster: Armstrong v.s. Clone13](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/perCluster_Arm_vs_Cl13)
     - [Per day: Armstrong v.s. Clone13](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/2_DE/perTimepoint_Arm_vs_Cl13/)
   
   #### 1.3 GSEA per condition (based on avg expression)
   - [Notebook](0_Acute-Chronic/codes_local/2_0_GSEA_resampled_avg_expr.ipynb)
   - [Output](0_Acute-Chronic/1_Scanpy/0_Scanpy_out_resampled/3_GSEA/)
   
## 2. scVelo <br>

   #### 2.1 scVelo with all resampled cells
   - [Notebook](0_Acute-Chronic/codes_local/3_0_scVelo_Analysis_resampled.ipynb)
   = [Data output](0_Acute-Chronic/2_scVelo/0_scVelo_out_resampled/)
   
   #### 2.2 scVelo with only activated (no naive)
   - [Notebook](0_Acute-Chronic/codes_local/3_1_scVelo_Analysis_ACTonly.ipynb)
   - [Data output](0_Acute-Chronic/2_scVelo/1_scVelo_out_ACTonly/
   
   #### 2.3 scVelo with naive & Armstrong
   - [Notebook](0_Acute-Chronic/codes_local/3_2_scVelo_Analysis_ARMonly.ipynb)
   - [Data output](0_Acute-Chronic/2_scVelo/2_scVelo_out_ARMonly/)

   #### 2.4 scVelo with naive & Clone13
   - [Notebook](0_Acute-Chronic/codes_local/3_3_scVelo_Analysis_CL13only.ipynb)
   - [Data output](0_Acute-Chronic/2_scVelo/3_scVelo_out_CL13only/)
   
   #### 2.5 scVelo with selected comparisons
   > Each cluster v.s. other clusters <br>
   > In Armstrong & Naive only: Each cluster v.s. other clusters <br>
   > In Clone13 & Naive only: Each cluster v.s. other clusters <br>
   - [Notebook](0_Acute-Chronic/codes_local/3_4_scVelo_Analysis_LouvainCluster.ipynb)
   - [Data output](0_Acute-Chronic/2_scVelo/4_scVelo_out_LouvainCluster/)

```                                                                      
                        ┌────────────────┐                                  
                        │   10x output   │                                  
                        └────────────────┘                                  
                                 │  cellranger                              
                                 ▼    mkfastq                               
                           ┌ ─ ─ ─ ─ ─                                      
                              Fastq   │                                     
                           └ ─ ─ ─ ─ ─                                      
                                 │  cellranger                              
                                 ▼     count                                
                       ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─                                  
                           Count Matrix   │                                 
                       │      10x QC                                        
                        ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘                                 
   Preprocess ━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━       
                                 │      Load into                           
                                 │   Scanpy / scVelo                        
                      ┌ ─ ─ ─ ─ ─▼─ ─ ─ ─ ─ ─ ─ ─ ─                         
                    ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  │                        
                  ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  │                          
                        Single Cell Object     │   │                        
                  │                              │                          
                             - Filter          │   │                        
                  │         - Resample         │ │                         
                   - PAGA - Louvain clustering │   │                        
                  │  - Velocity & transition   │ │─                         
                                               │─                           
                  └ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─                             
                                 │                      
              ┌──────────────────┼────────────────────────┐                 
              │                  │                        │                 
              ▼                  ▼                        ▼                 
     ┌─────────────────┐    ┌ ─ ─ ─ ─ ┐           ┌ ─ ─ ─ ─ ─ ─ ─           
     │    Plotting     │       Expr                 Differential │          
     │   Projection    │    │ matrix  │           │     Expr                
     │                 │     ─ ─ ─ ─ ─             ─ ─ ─ ─┐─ ─ ─ ┘          
     │UMAP single cell │         │ Sc Signature           │  GSEA           
     │PAGA single cell │         │   Explorer             │ ┌──────────────┐
     │  PAGA clusters  │         │  ┌──────────────┐      │ │ Single Cell  │
     └─────────────────┘         │  │ Single Cell  │      └▶│  Signature   │
                                 ├─▶│  Signature   │        └──────────────┘
                                 │  └──────────────┘                        
                                 │ GSEA + Cluster info                      
                                 │  ┌──────────────┐                        
                                 │  │   Cluster    │                        
                                 └─▶│  Signature   │                        
                                    └──────────────┘                        
```
