##############################################################################
### ### ### ### Deprecated functions. Updated script in Exp39 2### ### ### ### 
##############################################################################

######################################## Libraries ########################################
library(dplyr)
library(tidyverse)
library(fitdistrplus)


######################################## Self-defined functions ########################################
# custom function to transpose while preserving names 
# Thanks Indrajeet Patil @Stackoverflow https://stackoverflow.com/questions/42790219/how-do-i-transpose-a-tibble-in-r
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

nbQT_calc <- function(in_num_vec){
  in_num_vec <- in_num_vec*10
  in_num_vec <- round(in_num_vec)
  fit.nb <- NA
  quantiles <- rep(NA, length(in_num_vec))
  tryCatch((fit.nb <- fitdist(in_num_vec, "nbinom")), error=function(e) {})
  if (!(is.na(fit.nb))) {
    quantiles <- pnbinom(in_num_vec, size=as.numeric(fit.nb$estimate['size']), mu=as.numeric(fit.nb$estimate['mu']))
    quantiles <- quantiles * 100
  }
  return(quantiles)
}

cluster_average <- function(in_tb, cluster_tb, out_name_base) {
  ### in_tb: colnames: genenames + cell_id
  ### cluster_tb: cell_id + types
  #in_file <- "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/1_norm_counts/all_norm_counts_named_c10_nbPctl_Z_naOmit.csv"
  #in_tb <- read_csv(in_file)
  #cluster_file <- "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/trial_D5only/paga_d5_obs.csv"
  #cluster_tb <- read_csv(cluster_file)
  #cluster_tb <- cluster_tb %>% dplyr::select(one_of(c("0", "louvain")))
  #colnames(cluster_tb) <- c("cell_id", "types")
  #out_name <- "all_norm_counts_named_c10_nbPctl_Z--paga_d5_avg"
  
  use.cells <- cluster_tb$cell_id
  use.tb <- in_tb %>% filter(cell_id %in% use.cells) %>% inner_join(cluster_tb, by = "cell_id")
  use.tb$cell_id <- NULL
  
  mean.tb <- use.tb %>% group_by(types) %>% summarise_all(mean)
  rownames(mean.tb) <- mean.tb$types
  mean.tb$types <- NULL
  mean.tb <- transpose_df(mean.tb)
  colnames(mean.tb) <- c("gene_name", colnames(mean.tb)[2:length(colnames(mean.tb))])
  mean.tb.z <- mean.tb %>% mutate_if(is.numeric, function(x) {as.numeric(scale(x))})
  
  out.name.mean <- paste(out_name_base, ".csv", sep="")
  out.name.z <- paste(out_name_base, "_z.csv", sep="")
  
  write_csv(mean.tb, out.name.mean)
  write_csv(mean.tb.z, out.name.z)
}

###############################################################################

###----- Calculate percentiles & Z-Scores for each gene
if (FALSE) {
  wk.dir <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA'
  setwd(wk.dir)

  norm.file <- "all_norm_counts_named_c10.csv"
  norm.tb <- read_csv(norm.file)
  #tox.reads <- table(norm.tb$Tox)
  #max(norm.tb$Tox)
  #tox.reads[names(tox.reads)=="0"] / length(norm.tb$Tox)
  
  genes <- colnames(norm.tb)[2:length(colnames(norm.tb))]
  norm.tb$cell_id <- norm.tb$X1
  norm.tb$X1 <- NULL
  
  # Percentile
  pctl.tb <- norm.tb
  pctl.tb <- pctl.tb %>% mutate_if(is.numeric, function(x) { nbQT_calc(x)})
  write_csv(pctl.tb, "all_norm_counts_named_c10_nbPctl.csv")
  
  # Z-Score
  z.tb <- pctl.tb
  z.tb$cell_id <- as.character(z.tb$cell_id)
  z.tb <- z.tb %>% mutate_if(is.numeric, function(x) {as.numeric(scale(x))})
  z.tb <- z.tb[colSums(!is.na(z.tb)) > 0]
  write_csv(z.tb, "all_norm_counts_named_c10_nbPctl_Z_naOmit.csv")
  
  #max(z.tb$Tox)
}

###----- Calculate group average of Z-Scores


# All cells - selected - rmWTNAV
if (FALSE) {
  wk.dir <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/1_norm_counts'
  setwd(wk.dir)
  
  obs.file <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/0_PAGA/all--numSlt-rmWTNAV/obs.csv'
  out.name <- "all_norm_counts_named_c10_nbPctl_Z_naOmit--allNumSltrmWTNAV"
  obs.tb <- read_csv(obs.file) %>% dplyr::select(one_of(c("X1", "louvain")))
  colnames(obs.tb) <- c("cell_id", "types")
  
  in.file <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all_norm_counts_named_c10_nbPctl_Z_naOmit.csv'
  in.tb <- read_csv(in.file)
  
  cluster_average(in.tb, obs.tb, out.name)
}

# All cells - selected - rmWTNAV ----- Type & Louvain
if (TRUE) {
  wk.dir <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/1_norm_counts'
  setwd(wk.dir)
  
  obs.file <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/0_PAGA/all--numSlt-rmWTNAV/obs.csv'
  out.name <- "all_norm_counts_named_c10_nbPctl_Z_naOmit--allNumSltrmWTNAV_typeLouvain"
  obs.tb <- read_csv(obs.file) %>% 
    dplyr::select(one_of(c("X1", "louvain", "cell_type"))) %>%
    mutate(cell_type = substr(cell_type, 1, 1)) %>%
    unite("types", louvain:cell_type)
  colnames(obs.tb) <- c("cell_id", "types")
  
  in.file <- '/media/pipkin/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all_norm_counts_named_c10_nbPctl_Z_naOmit.csv'
  in.tb <- read_csv(in.file)
  
  cluster_average(in.tb, obs.tb, out.name)
} 






###----- Calculate percentiles & Z-Scores for each gene
if (FALSE) {
  wk.dir <- "/Volumes/Yolanda/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/2_norm_counts-expr_DESeq2/all--numslt-rmWTNAV"
  setwd(wk.dir)
  
  norm.file <- "all_norm_counts_named_c10_int_slt.csv"
  norm.tb <- read_csv(norm.file)
  norm.file.c10 <- gsub("_c10", "", norm.file)
  norm.file.c10 <- gsub(".csv", "_c10.csv", norm.file.c10)
  norm.file.c10.pctl <- gsub(".csv", "_nbPctl.csv",  norm.file.c10)
  norm.file.c10.pctl.z <- gsub(".csv", "_Z_naOmit.csv",  norm.file.c10.pctl)
  
  ### Select genes with count >= 10
  norm.tb.cellbcs <- norm.tb$X1
  norm.tb$X1 <- NULL
  gene.max.counts <- norm.tb %>% summarize_all(funs(max))
  gene.max.counts <- transpose_df(gene.max.counts)
  colnames(gene.max.counts) <- c("gene_name", "norm_count")
  gene.max.counts <- gene.max.counts %>% filter(norm_count >= 10)
  norm.tb <- norm.tb %>% dplyr::select(one_of(gene.max.counts$gene_name))
  
  genes <- colnames(norm.tb)
  norm.tb$cell_id <- norm.tb.cellbcs
  
  write_csv(norm.tb, norm.file.c10)
  
  # Percentile
  pctl.tb <- norm.tb
  pctl.tb <- pctl.tb %>% mutate_if(is.numeric, function(x) { nbQT_calc(x)})
  pctl.tb <- pctl.tb[colSums(!is.na(pctl.tb)) > 0]
  write_csv(pctl.tb, norm.file.c10.pctl)
  
  # Z-Score
  z.tb <- pctl.tb
  z.tb$cell_id <- as.character(z.tb$cell_id)
  z.tb <- z.tb %>% mutate_if(is.numeric, function(x) {as.numeric(scale(x))})
  z.tb <- z.tb[colSums(!is.na(z.tb)) > 0]
  write_csv(z.tb, norm.file.c10.pctl.z)
}


