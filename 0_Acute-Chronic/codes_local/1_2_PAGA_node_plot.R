######################################## PAGA node plots ########################################
# Author: Huitian (Yolanda) Diao
# Sept 11th, 2019

######################################## Libraries ########################################
library(dplyr)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(grid)

######################################## Self defined functions ########################################
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

conn_plot <- function(in_tb, out_name_root, val_name, fillcol, mute) {
  #in_tb <- in.tb
  #val_name <- "Z"
  #fillcol <- "grey"
  
  paga.conn.file.coor <- "/Users/yolandatiao/GSuite\ Scripps/Exp392_shCRF_SC/1_1_SCANPY_PAGA/shCRF-numSlt_fltType/0_PAGA/shCRF-numSlt_fltType_paga_connect_long_coor.csv"
  paga.pos.file <- "/Users/yolandatiao/GSuite\ Scripps/Exp392_shCRF_SC/1_1_SCANPY_PAGA/shCRF-numSlt_fltType/0_PAGA/shCRF-numSlt_fltType_paga_pos.csv"
  #paga.conn.file.coor <- "/Users/yolandatiao/Exp334CD25KOSc/0_codes/1_2_PAGA/uns/paga_connect_long_coor.csv"
  #paga.pos.file <- "/Users/yolandatiao/Exp334CD25KOSc/0_codes/1_2_PAGA/uns/paga_pos_renamed.csv"
  paga.pos.tb <- read_csv(paga.pos.file)
  colnames(paga.pos.tb) <- c("cluster_name", "x", "y")
  paga.pos.tb$cluster_name <- paste("L", as.character(paga.pos.tb$cluster_name), sep="")
  paga.conn.tb <- read_csv(paga.conn.file.coor)
  paga.conn.tb <- paga.conn.tb %>% filter(Conn > 0.25)
  paga.pos.tb <- paga.pos.tb %>% left_join(in_tb, by="cluster_name")
  
  max_z <- max(paga.pos.tb$val)
  min_z <- min(paga.pos.tb$val)
  max_z <- ceiling(max(abs(min_z), max_z)/2)*2
  min_z <- -max_z
  
  conn.plot <- ggplot() +
    geom_segment(data=paga.conn.tb, aes(x=x1, y=y1, xend=x2, yend=y2), size=paga.conn.tb$Conn*2.5, alpha=0.1) +
    geom_point(data=paga.pos.tb, aes(x=x, y=y, color=val), size=3.5) +
    scale_color_gradient2(val_name,high="firebrick4", low="dodgerblue4", mid="white", limits=c(min_z, max_z)) +
    scale_x_continuous(limits = c(-0.3, 1.5)) +
    scale_y_continuous(limits = c(-1.8, 1.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = fillcol),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.line.x = element_line(color="black", size = 0.25),
          axis.line.y = element_line(color="black", size = 0.25))
  
  if (!mute) {
    conn.plot.nolegend <- conn.plot +
      theme(legend.position="none")
    out_name <- paste(out_name_root, ".pdf", sep="")
    out_name_nolegend <- paste(out_name_root, "_nolegend.pdf", sep="")
    out_name_nolegend_png <- paste(out_name_root, "_nolegend.png", sep="")
    
    ggsave(out_name, conn.plot, device="pdf", width=5.8, height=4, units=('cm'), dpi=300)
    ggsave(out_name_nolegend, conn.plot.nolegend, device="pdf", width=4, height=4, units=('cm'), dpi=300)
    ggsave(out_name_nolegend_png, conn.plot.nolegend, device="png", width=4, height=4, units=('cm'), dpi=300)
  }
  conn.plot <- conn.plot + ggtitle(out_name_root) + theme(plot.title = element_text(face="bold", hjust = 0.5))
  return(conn.plot)
}

# Replace vector elements into new names
cvt_names <- function(vec_x, in_vec, out_vec){
  vec_x_out <- c()
  for (i in vec_x) {
    if (i %in% in_vec){
      i_idx <- match(i, in_vec)
      vec_x_out <- c(vec_x_out,out_vec[i_idx])
    } else {
      vec_x_out <- c(vec_x_out, i)
    }
  }
  return(vec_x_out)
}

# Convert names with underscore to names that match position file format
cvt_underscore_names <- function(in_vec, oldornewname){
  #in_vec <- wt.v.avg.i$new_name
  #oldornewname <- "new"
  out_vec <- c()
  for (i in in_vec){
    if (oldornewname == "old") {
      out_vec <- c(out_vec, paste("P", as.character(unlist(strsplit(i, split="_")))[2], sep=""))
    } else {
      out_vec <- c(out_vec, as.character(gsub("p", "", unlist(strsplit(i, split="_"))[2])))
    }
  }
  return(out_vec)
}

######################################## Main ########################################
#####---------- Pre-process paga connection file
if (FALSE) {
  paga.out.dir <- "/Users/yolandatiao/GSuite\ Scripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/0_PAGA"
  paga.conn.file <- "all--numSlt-rmWTNAV_paga_connect.csv"
  paga.conn.file <- file.path(paga.out.dir, paga.conn.file)
  paga.conn.file.long <- gsub(".csv", "_long.csv", paga.conn.file)
  paga.conn.tb <- read_csv(paga.conn.file)
  paga.conn.tb <- paga.conn.tb %>% melt(id.vars="X1")
  colnames(paga.conn.tb) <- c("From", "To", "Conn")
  paga.conn.tb$From <- as.numeric(paga.conn.tb$From)
  paga.conn.tb$To <- as.numeric(as.character(paga.conn.tb$To))
  nrow(paga.conn.tb)
  paga.conn.tb <- paga.conn.tb %>% rowwise() %>% mutate(fromto = paste(as.character(sort(c(From, To))), collapse="_"))
  paga.conn.tb$From <- NULL
  paga.conn.tb$To <- NULL
  paga.conn.tb <- paga.conn.tb %>% distinct()
  paga.conn.tb$fromto
  paga.conn.tb$Conn
  paga.conn.tb <- paga.conn.tb %>% separate(fromto, c("node1", "node2"), sep="_")
  paga.conn.tb <- paga.conn.tb %>% filter(Conn > 0)
  paga.conn.tb <- paga.conn.tb %>% mutate(node1 = paste("L", as.character(node1), sep="")) %>% mutate(node2 = paste("L", as.character(node2), sep=""))
  
  write_csv(paga.conn.tb, paga.conn.file.long)
}

# Get coordinates for edges
if (FALSE) {
  paga.out.dir <- "/Users/yolandatiao/GSuite Scripps/Exp391_Acute-Chronic_SC/1_1_SCANPY_PAGA/all--numSlt-rmWTNAV/0_PAGA"
  setwd(paga.out.dir)
  paga.conn.file.long <- "all--numSlt-rmWTNAV_paga_connect_long.csv"
  paga.conn.file.long <- file.path(paga.out.dir, paga.conn.file.long)

  paga.pos.file <- "all--numSlt-rmWTNAV_paga_pos.csv"
  paga.pos.file <- file.path(paga.out.dir, paga.pos.file)

  paga.pos.tb <- read_csv(paga.pos.file)
  paga.conn.tb <- read_csv(paga.conn.file.long)
  paga.pos.tb$X1 <- paste("L", as.character(paga.pos.tb$X1), sep="")

  paga.conn.file.coor <- gsub(".csv", "_coor.csv", paga.conn.file.long)
  colnames(paga.pos.tb) <- c("node1", "x1", "y1")
  paga.conn.tb <- paga.conn.tb %>% left_join(paga.pos.tb)
  colnames(paga.pos.tb) <- c("node2", "x2", "y2")
  paga.conn.tb <- paga.conn.tb %>% left_join(paga.pos.tb)
  write_csv(paga.conn.tb, paga.conn.file.coor)
}

#####---------- Plotting by genes (Selected genes -- P14 WT only)
if (FALSE) {
  wk.dir <- "/Users/yolandatiao/GSuite\ Scripps/Exp392_shCRF_SC/1_1_SCANPY_PAGA/shCRF-numSlt_fltType/2_Plots/1_nodeplot"
  setwd(wk.dir)
  
  paga.z.file <- "/Users/yolandatiao/GSuite\ Scripps/Exp392_shCRF_SC/1_1_SCANPY_PAGA/shCRF-numSlt_fltType/1_norm_counts/all_norm_counts_named_c10_nbPctl_Z_naOmit--sltNumFltType_Louvainavg_z.csv"
  paga.z.tb <- read_csv(paga.z.file)
  
  use.genes <- c("Tcf7",  "Id3", "Id2", "Klrg1", "Sell", "Cd44", "Tox", "Kmt2a", "Cx3cr1") #"Bcl6" not found
  for (gene in use.genes) {
    #gene <- "Prdm1"
    in.tb <- paga.z.tb %>% filter(gene_name == gene)
    in.tb$gene_name <- NULL
    in.tb <- transpose_df(in.tb)
    colnames(in.tb) <- c("cluster_name", "val")
    in.tb <- in.tb %>% mtate(cluster_name = paste("L", as.character(cluster_name), sep=""))
    out.name.root <- paste(gene, "_z", sep="")
    conn_plot(in.tb, out.name.root, "Z-Score", "grey", FALSE)
  }
  
}

