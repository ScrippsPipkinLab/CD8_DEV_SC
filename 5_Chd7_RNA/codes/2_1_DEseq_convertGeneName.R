################################################################################
#             Convert ENSEMBL-IDs to GeneNames for DEseq output                #
#                              Dependent: DEseq                                #
################################################################################

library(biomaRt)

setwd("/Volumes/Huitian/Exp276/DEseq")

#the data set at ensembl that contains all of their mouse data is "mmusculus_gene_ensembl"
#generating a datafram of all the possible attributes you can download from biomaRt
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#attributes_df<-listAttributes(biomart="ensembl")



matchGN <- function(infilename, slicing)
{
  outfilename <- gsub(".csv","_gn.csv",infilename)
  
  #infilename <- "/Volumes/Huitian/GSE66763/salmon/sum/Expressed/GSE667343_count_exp.csv"
  #infilename <- "/Volumes/Huitian/GSE66763/DEseq/HP-HN_Salmon.csv"
  #outfilename <- "/Volumes/Huitian/GSE66763/salmon/Names_gn.csv"
  
  #import our gene information and query against the biomaRt to find our desired attributes
  input <- read.csv(infilename, header = TRUE)
  esblNames <- input[1]
  mart <- useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl"))
  esblNames
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=esblNames,mart= mart)
  output <- merge(input, G_list, by.x="X", all.x=TRUE, by.y="ensembl_gene_id")
  write.csv(output, outfilename)
}


setwd("/Volumes/Huitian/Exp276/DEseq")
if (TRUE){
  file.names <- dir("/Volumes/Huitian/Exp276/DEseq",pattern=".csv")
  for (i in file.names)
  {
    matchGN(i, FALSE)
  }
}


