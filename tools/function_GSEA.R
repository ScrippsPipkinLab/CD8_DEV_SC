######################################## Average Z Score -- GSEA ########################################
# Author: Huitian (Yolanda) Diao
# Feb 5th, 2019
# Update July 2nd, 2020

######################################## Libraries ########################################
library("org.Mm.eg.db")
library(DOSE)
require(clusterProfiler)
library(dplyr)
library(enrichplot)
library(tidyverse)

le_to_le_sig <- function(vec_x) {
  out_vec <- c()
  for (i in vec_x){
    #i <- "tags=58%, list=12%, signal=52%"
    i_vec <- unlist(strsplit(i, ", "))
    i_sig <- tail(i_vec, n=1)
    i_sig <- gsub("signal=", "", i_sig)
    i_sig <- gsub("%", "", i_sig)
    i_sig <- as.numeric(i_sig)
    out_vec <- c(out_vec, i_sig)
  }
  return(out_vec)
}

filename_simp <- function(in_filename) {
  cpx <- gsub(".csv", "", in_filename)
  cpx <- gsub("louvain_", "C", cpx)
  cpx <- unlist(strsplit(cpx, "---"))[1]
  return(cpx)
}

GSEA_analysis <- function(in.file, use.group, out_name, gs_file, zscore_cutoff) {
  file.name.simp <- out_name
  in.tb <- read_csv(in.file)
  use.col <- c("gene_name", "gene_names", use.group)
  use.tb <- in.tb %>% dplyr::select(one_of(use.col))
  colnames(use.tb) <- c("gene_name", "z")
  
  #####---------- Sample gsea analysis with custome gene list
  gene.list <- use.tb$z
  names(gene.list) <- as.character(use.tb$gene_name)
  gene.list <- sort(gene.list, decreasing = TRUE)
  deg.list <- names(gene.list)[abs(gene.list) > zscore_cutoff]
  print(length(deg.list))
  
  #####---------- Read GSEA reference dataset
  gs.tb <- read_csv(gs_file)
  unique( gs.tb$gs_name)
  
  gs.name.simp <- unlist(strsplit(gs_file, "/"))
  gs.name.simp <- tail(gs.name.simp, 1)
  gs.name.simp <- gsub(".csv", "", gs.name.simp)
  
  #####---------- RUN GSEA
  em <- enricher(deg.list, TERM2GENE=gs.tb)
  em2 <- GSEA(gene.list, TERM2GENE=gs.tb,  nPerm = 10000, minGSSize = 1, maxGSSize = 5000,  pvalueCutoff = 1, by="DOSE")
  
  #####---------- Export results
  tb.name <- paste(file.name.simp, "---", gs.name.simp, ".csv", sep="")
  results.tb <- em2@result
  #print(head(results.tb))
  write_csv(results.tb, tb.name)
  gs <- em2@result$ID
  for (i in c(1:length(gs))){
    gsplot.i.name <- paste(file.name.simp, "---", gs.name.simp,"___",  gs[i], ".pdf", sep="")
    gsplot.i <- gseaplot2(em2, geneSetID = i, title = gs[i])
    ggsave(gsplot.i.name, gsplot.i, device="pdf", width=15, height=10, units="cm")
  }
  
}

GSEA_sum <- function(file_list, outname, wid, hei, 
                     comp_use_order, comp_order_vec, path_use_order, path_order_vec,
                     filesimp_use_selfdefined, filesimp_selfdefined, use_pval) {
  ########## Parameters ##########
  ### file_list: list of gsea output csv files
  #-- comparison names derive from file names automatically
    
  ### outname: base of output name
  ### wid: width of output plot
  ### hei: height of output plot
    
  ### comp_use_order: if using custom order for comparisons
  #-- comp_order_vec: the order of comparisons
  
  ### path_use_order: if using custom order for pathways
  #-- path_order_vec: the order of pathways to use
  
  ### filesimp_use_selfdefined: if using custom defined simplified file names
  #-- filesimp_selfdefined: self defined simplified file names
  
  ### use_pval: if using pvalue instead of padj
  
  ########## Output ##########
  ### Returns summary table
  #-- comparison, pathway, NES, padj, leadingEdge_signal, mlog10padj
    
  ### Save summary bubble plot
  #-- size=mlog10padj, color=NES
  #-- ordered by specified comparison & pathway
    
  cp.vec <- c()
  pw.vec <- c()
  nes.vec <- c()
  padj.vec <- c()
  pval.vec <- c()
  le.sig.vec <- c()
  for (x in c(1:length(file_list))) {
    filex <- file_list[x]
    filex.name.simp <- basename(filex)
    filex.name.simp <- filename_simp(filex.name.simp)
    if (filesimp_use_selfdefined) {
        filex.name.simp  <- filesimp_selfdefined[x]
    }
    
    filex.tb <- read_csv(filex)
    cp.vec <- c(cp.vec, rep(filex.name.simp, nrow(filex.tb)))
    pw.vec <- c(pw.vec, filex.tb$ID)
    nes.vec <- c(nes.vec, filex.tb$NES)
    padj.vec <- c(padj.vec, filex.tb$p.adjust)
    pval.vec <- c(pval.vec, filex.tb$pvalue)
    le.x <- filex.tb$leading_edge
    le.sig.x <- le_to_le_sig(le.x)
    le.sig.vec <- c(le.sig.vec, le.sig.x)
  }
  
  plot.tb <- tibble(comparison=cp.vec,
                    pathway=pw.vec,
                    NES=nes.vec,
                    padj=padj.vec,
                    pval=pval.vec,
                    leadingEdge_signal=le.sig.vec)
  plot.tb$mlog10padj <- -log10(plot.tb$padj)
  plot.tb$mlog10pval <- -log10(plot.tb$pval)
    
  print(head(plot.tb))
  
  if (comp_use_order) {
    new.factor <- factor(plot.tb$comparison, levels=comp_order_vec)
    plot.tb$comparison <- new.factor
  }
  if (path_use_order) {
    new.factor <- factor(plot.tb$pathway, levels=path_order_vec)
    plot.tb$pathway <- new.factor
  }
  
  if (! use_pval) {
    bbplot <- ggplot(plot.tb, aes(pathway, comparison)) +
    geom_point(aes(size=mlog10padj, color=NES)) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    bbplot <- ggplot(plot.tb, aes(pathway, comparison)) +
    geom_point(aes(size=mlog10pval, color=NES)) +
    scale_color_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }

  
  outname_csv <- paste(outname, ".csv", sep="")
  outname_pdf <- paste(outname, ".pdf", sep="")
  
  write_csv(plot.tb, outname_csv)
  ggsave(outname_pdf, bbplot, device="pdf", width=wid, height=hei, units="cm")
    
  return(plot.tb)
}

GO_run <- function(z_file, z_cutoff){
  #z.file <- "/Volumes/Yolanda/Exp392_shCRF_SC/1_1_SCANPY_PAGA/shCRF-numSlt/1_norm_counts/all_norm_counts_named_c10_nbPctl_Z_naOmit--louvain_avg_z.csv"
  #z_cutoff <- 1 # For testing
  
  tab.i <- read_csv(z.file)
  conds <- colnames(tab.i)[2:length(colnames(tab.i))]
  for (i.cond in conds) {
    i <- as.character(i.cond)
    genes.i <- tab.i %>% filter(tab.i[i.cond] > z_cutoff) %>% .$gene_name
    print(paste(as.character(i.cond), "    Gene number: ", as.character(length(genes.i)), sep=""))
    genes.i.id <- AnnotationDbi::select(org.Mm.eg.db, genes.i, c("ENTREZID"), "ALIAS")
    
    egoBP <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoCC <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "CC", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    egoMF <- enrichGO(gene = genes.i.id$ENTREZID, keyType = 'ENTREZID', OrgDb = org.Mm.eg.db, ont = "MF", pAdjustMethod = "none", pvalueCutoff = 0.05, readable = TRUE)
    
    # Dotplot visualization
    if (!is.null(egoBP)){
      pdf.name <- paste(i,"_BP_dotplot.pdf",sep="")
      csv.name <- paste(i,"_BP_dotplot.csv",sep="")
      write.csv(egoBP@result, file=csv.name, row.names=FALSE)
      egoBP.dotplot <- dotplot(egoBP, showCategory=25)
      ggsave(pdf.name, egoBP.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
      
    }
    if(!is.null(egoCC)){
      csv.name <- paste(i,"_CC_dotplot.csv",sep="")
      pdf.name <- paste(i,"_CC_dotplot.pdf",sep="")
      write.csv(egoCC@result, file=csv.name, row.names=FALSE)
      egoCC.dotplot <- dotplot(egoCC, showCategory=25)
      ggsave(pdf.name, egoCC.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
    }
    if(!is.null(egoMF)){
      csv.name <- paste(i,"_MF_dotplot.csv",sep="")
      pdf.name <- paste(i,"_MF_dotplot.pdf",sep="")
      write.csv(egoMF@result, file=csv.name, row.names=FALSE)
      egoMF.dotplot <- dotplot(egoMF, showCategory=25)
      ggsave(paste(i,"_MF_dotplot.pdf",sep=""), egoMF.dotplot, device = "pdf", width = 30, height = 20, units = "cm")  
    }
  }
}

######################################## Main ########################################
# GO terms
if (FALSE) {
  wk.dir <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/3_GO/typeLouvain'
  setwd(wk.dir)
  
  z.file <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/1_norm_counts/all_norm_counts_named_c10_nbPctl_Z_naOmit--allNumSltrmWTNAV_typeLouvain_z.csv'
  
  GO_run(z.file, 1)
}
