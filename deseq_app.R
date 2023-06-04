library(DESeq2)
library(tidyverse)



setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/DESEQ/")


coldata <- read.csv("old-app-cortex-plaque_labels.csv", row.names=1)
cnts <- read.csv("old-app-cortex-plaque_counts.csv", header=TRUE, check.names=FALSE, row.names=1)

dds <-DESeqDataSetFromMatrix(countData=cnts, 
                       colData=coldata, 
                       design=~sample+plaque) 

dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), 
          file="old-app-cortex-plaque_DESEQ-pval.csv")




