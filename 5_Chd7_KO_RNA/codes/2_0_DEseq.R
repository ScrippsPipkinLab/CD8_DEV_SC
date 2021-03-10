################################################################################
#                   Extracting from Salmon quant & DEseq                       #
#                  Original version: Adam Getzler 7/18/17                      #
#                               Edits:Huitian Diao                             #
################################################################################

library(readr)
library(GenomicFeatures)
library(DESeq2)
library(pcaExplorer)
library(org.Mm.eg.db)
library(rjson)
library(tximport)
library(DBI)

##########---------- Setworking directory
setwd("/Volumes/Huitian/Exp276/salmon_nonTrim")

##########---------- Read Quant Files
colData <- read.csv("/Volumes/Huitian/Exp276/meta.csv")
files <- file.path("/Volumes/Huitian/Exp276/salmon_nonTrim",colData$Samples,"quant.sf")
names(files) <- colData$Cond


##########---------- Make table to match EMSU to GeneID
#setwd("/Volumes/Huitian/GSE66763")
#txdb <- makeTxDbFromGFF("/Volumes/Huitian/GSE66763/Homo_sapiens.GRCh38.92.gtf.gz")
#saveDb(txdb, file="Homo_sapiens.GRCh38.92") #Did this once so I dont need to do it again

txdb <- makeTxDbFromGFF("/Volumes/Huitian/Exp276/Mus_musculus.GRCm38.90.gtf.gz")
saveDb(txdb, file="Mus_musculus.GRCm38.90") #Did this once so I dont need to do it again


###--- Convert transcript ID to gene ID
txdb <- loadDb("Mus_musculus.GRCm38.90")
k <- keys(txdb, "GENEID")
res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
tx2gene <- res[,2:1]
# Must make table TX name Gene id to be read
tx2gene

# Drop in freps TURE = ignore verison  # Ignore TX verison stringsplits on . 
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE) #dropInfReps = TRUE, 

# Construct sampleTable
sampleTable <- data.frame(condition =factor(c("WT_MP","WT_TE","WT_EE","WT_DP","WT_MP","WT_TE","WT_EE","WT_DP",
                                              "KO_MP","KO_TE","KO_EE","KO_DP","KO_MP","KO_TE","KO_EE","KO_DP",
                                              "WT_48h","WT_48h","KO_48h","KO_48h","WT_NAV","KO_NAV","WT_NAV","KO_NAV")))
rownames(sampleTable) <- colnames(txi$counts)

#import into DESEQ2 framework
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
summary(dds)


dds <- DESeq(dds) #RunDESEQ

results<-as.data.frame(results(dds, contrast = c("condition","KO_MP","WT_MP"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_MP-WT_MP-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_TE","WT_TE"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_TE-WT_TE-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_EE","WT_EE"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_EE-WT_EE-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_DP","WT_DP"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_DP-WT_DP-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_48h","WT_48h"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_48h-WT_48h-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_NAV","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_NAV-WT_NAV-Salmon.csv")


### WT-WT
results<-as.data.frame(results(dds, contrast = c("condition","WT_48h","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"WT_48h-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","WT_MP","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"WT_MP-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","WT_TE","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"WT_TE-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","WT_EE","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"WT_EE-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","WT_DP","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"WT_DP-WT_NAV-Salmon.csv")

### KO-KO
results<-as.data.frame(results(dds, contrast = c("condition","KO_48h","KO_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_48h-KO_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_MP","KO_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_MP-KO_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_TE","KO_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_TE-KO_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_EE","KO_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_EE-KO_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_DP","KO_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_DP-KO_NAV-Salmon.csv")

### KO-WT
results<-as.data.frame(results(dds, contrast = c("condition","KO_48h","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_48h-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_MP","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_MP-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_TE","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_TE-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_EE","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_EE-WT_NAV-Salmon.csv")

results<-as.data.frame(results(dds, contrast = c("condition","KO_DP","WT_NAV"))) # makes the fold change 1st elment is the numerator 2nd denominator 
write.csv(results,"KO_DP-WT_NAV-Salmon.csv")
