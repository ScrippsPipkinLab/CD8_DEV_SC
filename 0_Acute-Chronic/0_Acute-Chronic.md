
# Exp391: Acute v.s. Chronic infection

# Results:

### 1. Scanpy <br>
   > With filtered resampled cells (capped 1250 cells / category) <br>
   
   #### 1.1 Scanpy & PAGA 
   - [Notebook](0_Acute-Chronic/codes_local/1_0-0_SCANPY_PAGA_resampled.ipynb)
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
