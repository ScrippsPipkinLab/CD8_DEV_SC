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

######################################## PAGA node plots ########################################


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

# Node plot
conn_plot <- function(paga.pos.tb, paga.conn.tb, out_name_root, var_name = "NES",fillcol="grey", mute=TRUE,
                     heatbar.range=c(-4,4)) {
  # paga.pos.tb: [x, y, val]
  # paga.conn.tb: [x1, y1, x2, y2, Conn]
  
  # Calculate range of plotting
  max_x <- max(paga.pos.tb$x)
  min_x <- min(paga.pos.tb$x)
  max_y <- max(paga.pos.tb$y)
  min_y <- min(paga.pos.tb$y)
  range_x <- max_x - min_x
  range_y <- max_y - min_y
  min_x <- min_x - range_x*0.1
  max_x <- max_x + range_x*0.1
  min_y <- min_y - range_y*0.1
  max_y <- max_y + range_y*0.1
  
  conn.plot <- ggplot() +
    geom_segment(data=paga.conn.tb, aes(x=x1, y=y1, xend=x2, yend=y2), size=paga.conn.tb$Conn*1.5, alpha=0.1) +
    geom_point(data=paga.pos.tb, aes(x=x, y=y, color=val, size=size)) +
    scale_color_gradient2(var_name,high="firebrick4", low="dodgerblue4", mid="white", 
                          limits=c(heatbar.range[1], heatbar.range[2])) +
    scale_size("nlog10padj",range=c(0,4), breaks = c(0,1,2,3),labels=c('0','1','2', '>=3')) +
    scale_x_continuous(limits = c(min_x, max_x)) +
    scale_y_continuous(limits = c(min_y, max_y)) +
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
    
    ggsave(out_name, conn.plot, device="pdf", width=10, height=7.5, units=('cm'), dpi=300)
    ggsave(out_name_nolegend, conn.plot.nolegend, device="pdf", width=4, height=4, units=('cm'), dpi=300)
    ggsave(out_name_nolegend_png, conn.plot.nolegend, device="png", width=4, height=4, units=('cm'), dpi=300)
  }
  conn.plot <- conn.plot + ggtitle(out_name_root) + theme(plot.title = element_text(face="bold", hjust = 0.5))
  return(conn.plot)
}