
# Exp391: Acute v.s. Chronic infection



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
