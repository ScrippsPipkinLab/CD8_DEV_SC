library(tidyverse)
library(ggplot2)
library(reshape2)

#####---------- Pre-process differential analysis results -- expression level
preprocss_de_compile_expr <- function(gs_order_file, de_out_dir) {
  #gs.order.file <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/z_References/GSEA/Combined_GEO_GSEA_gs_order.csv"
  #de.out.dir <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/4_DE_GSEA/Cluster_vs_Cluster"
  
  gsea.pval.tb <- read_csv(gs_order_file)
  gsea.NES.tb <- read_csv(gs_order_file)
  
  de.files <- list.files(de_out_dir, pattern="*signatures.csv", recursive=TRUE, full.names=TRUE)
  for (file.i in de.files){
    name.i <- basename(file.i)
    name.i <- gsub("---Combined_GEO_GSEA_signatures.csv", "", name.i)
    i.pval.tb <- read_csv(file.i) %>% dplyr::select(one_of("ID", "pvalue"))
    i.NES.tb <- read_csv(file.i) %>% dplyr::select(one_of("ID", "NES"))
    colnames(i.pval.tb) <- c("gs_name", name.i)
    colnames(i.NES.tb) <- c("gs_name", name.i)
    gsea.pval.tb <- gsea.pval.tb %>% left_join(i.pval.tb, by="gs_name")
    gsea.NES.tb <- gsea.NES.tb %>% left_join(i.NES.tb, by="gs_name")
  }
  
  out.name <- file.path(de_out_dir, "compiled_gsea_pval.csv")
  write_csv(gsea.pval.tb, out.name)
  
  out.name <- file.path(de_out_dir, "compiled_gsea_NES.csv")
  write_csv(gsea.NES.tb, out.name)
}


#####---------- Pre-process differential analysis results -- velocity level
preprocess_de_compile_v <- function(gs_order_file, de_out_dir) {
  #gs.order.file <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/z_References/GSEA/Combined_GEO_GSEA_gs_order.csv"
  #de.out.dir <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_2_scVelo/2_GSEA/Cluster_vs_Cluster"
  
  gsea.pval.tb <- read_csv(gs_order_file)
  gsea.NES.tb <- read_csv(gs_order_file)
  
  de.files <- list.files(de_out_dir, pattern="*signatures.csv", recursive=TRUE, full.names=TRUE)
  for (file.i in de.files){
    name.i <- basename(file.i)
    name.i <- gsub("---Combined_GEO_GSEA_signatures.csv", "", name.i)
    i.pval.tb <- read_csv(file.i) %>% dplyr::select(one_of("ID", "pvalue"))
    i.NES.tb <- read_csv(file.i) %>% dplyr::select(one_of("ID", "NES"))
    colnames(i.pval.tb) <- c("gs_name", name.i)
    colnames(i.NES.tb) <- c("gs_name", name.i)
    gsea.pval.tb <- gsea.pval.tb %>% left_join(i.pval.tb, by="gs_name")
    gsea.NES.tb <- gsea.NES.tb %>% left_join(i.NES.tb, by="gs_name")
  }
  
  conditions <- unique( gsub(".y", ".x", colnames(gsea.pval.tb)))
  
  gsea.pval.tb <- gsea.pval.tb %>% select(one_of(conditions))
  gsea.NES.tb <- gsea.NES.tb %>% select(one_of(conditions))
  colnames(gsea.pval.tb) <- gsub(".x", "", colnames(gsea.pval.tb))
  colnames(gsea.NES.tb) <- gsub(".x", "", colnames(gsea.NES.tb))
  colnames(gsea.pval.tb) <- gsub("-", "_vs_", colnames(gsea.pval.tb))
  colnames(gsea.NES.tb) <- gsub("-", "_vs_", colnames(gsea.NES.tb))
  
  out.name <- file.path(de_out_dir, "compiled_gsea_pval.csv")
  write_csv(gsea.pval.tb, out.name)
  
  out.name <- file.path(de_out_dir, "compiled_gsea_NES.csv")
  write_csv(gsea.NES.tb, out.name)
}

#####---------- GSEA nodeplot
# [gs_file] gs_name: GSEA signature names
# [conn_file] Conn: PAGA connection score; comparison: X_vs_Y; x1 x2 y1 y2: node coordinate
# [nes_file] comparison: X_vs_Y; columns: GSEA siganture names
# [pval_file] comparison: X_vs_Y; columns: GSEA siganture names
# [gs_use] "gsA,gsB": return plots; if not specified not return
# create: plot dataframe (.csv); plot (.png)
gsea_nodeplot <- function(gs_file, conn_file, nes_file, pval_file, gs_use="all"){
  #gs_file <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/z_References/GSEA/Combined_GEO_GSEA_gs_order.csv"
  #conn_file <- "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/0_PAGA/all--numSlt-rmWTNAV_paga_connect_long_coor_newlabel.csv"
  #gsea.nes.file <-  "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/4_DE_GSEA/Cluster_vs_Cluster/compiled_gsea_NES.csv"
  #gsea.pval.file <-  "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/4_DE_GSEA/Cluster_vs_Cluster/compiled_gsea_pval.csv"
  #nes_file <-  "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_2_scVelo/2_GSEA/Cluster_vs_Cluster/compiled_gsea_NES.csv"
  #pval_file <-  "/Users/yolandatiao/GSuiteScripps/Exp391_Acute-Chronic_SC/1_2_scVelo/2_GSEA/Cluster_vs_Cluster/compiled_gsea_pval.csv"
  
  conn.cutoff <- 0.3
  
  gs.tb <- read_csv(gs_file)
  paga.conn.tb <- read_csv(conn_file) %>% 
    filter(Conn >= conn.cutoff) %>%
    unite(comparison, node1:node2, sep="_vs_")
  
  gsea.nes.tb <- as.tibble(t(read.csv(nes_file, row.names=1)), rownames = "comparison")
  gsea.pval.tb <- as.tibble(t(read.csv(pval_file, row.names=1)), rownames = "comparison")
  
  if (gs_use=='all') {
    gs_use_vec = gs.tb$gs_name
  } else {
    gs_use_vec = unlist(strsplit(gs_use, ","))
  }
  
  for (gs in gs_use_vec){
    #gs <- gs_use_vec[1]
    gsea.nes.plot <- gsea.nes.tb %>% select(one_of(c("comparison", gs)))
    gsea.pval.plot <- gsea.pval.tb %>% select(one_of(c("comparison", gs)))
    colnames(gsea.nes.plot) <- c("comparison", "NES")
    colnames(gsea.pval.plot) <- c("comparison", "pval")
    plot.tb <- paga.conn.tb %>% 
      left_join(gsea.nes.plot, by="comparison") %>% 
      left_join(gsea.pval.plot, by="comparison") %>%
      mutate(mlog10p = -log10(pval+0.000001))
    
    # save plot tb
    if (gs_use == "all") {
      out.name <- paste(gs, ".csv", sep="")
      write_csv(plot.tb, out.name)
    }
    
    ### Convert the start and end based on NES
    plot.tb <- plot.tb %>%
      mutate(NESpos = abs(NES)/NES)
    plot.tb.pos <- plot.tb %>% filter(NESpos == 1)
    plot.tb.neg <- plot.tb %>% filter(NESpos == -1)
    plot.tb.neg.x1 <- plot.tb.neg$x1
    plot.tb.neg.x2 <- plot.tb.neg$x2
    plot.tb.neg.y1 <- plot.tb.neg$y1
    plot.tb.neg.y2 <- plot.tb.neg$y2
    plot.tb.neg$x1 <- plot.tb.neg.x2
    plot.tb.neg$x2 <- plot.tb.neg.x1
    plot.tb.neg$y1 <- plot.tb.neg.y2
    plot.tb.neg$y2 <- plot.tb.neg.y1
    plot.tb <- plot.tb.pos %>% bind_rows(plot.tb.neg)
    
    out.plot <- ggplot(plot.tb) +
      geom_segment(data=plot.tb, aes(x=x1, y=y1, xend=x2, yend=y2), size=plot.tb$Conn*2.5, alpha=0.1) + 
      geom_segment(data=plot.tb, aes(x=x1, y=y1, xend=x2, yend=y2, color=mlog10p), arrow=arrow(length=unit(abs(plot.tb$NES)/6, "cm"), ends = "first"), size=abs(plot.tb$NES)) +
      scale_colour_gradient(low="gray", high="firebrick3") +
      ggtitle(gs) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.background = element_blank())
    
    if (gs_use == "all") {
      out.name <- paste(gs, ".png", sep="")
      ggsave(out.name, out.plot, width=9,height=6,units="cm") 
    } else {
      return(out.plot)
    }
  }
}

