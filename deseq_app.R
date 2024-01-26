library(DESeq2)
library(tidyverse)

setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/DESEQ/")

coldata <- read.csv("old-app-cortex-plaque-BULK_labels.csv", row.names=1)
cnts <- read.csv("old-app-cortex-plaque-BULK_counts.csv", header=TRUE, check.names=FALSE, row.names=1)

dds <-DESeqDataSetFromMatrix(countData=cnts, 
                       colData=coldata, 
                       design=~sex+plaque) 

dds <- DESeq(dds)
res <- results(dds, alpha=0.1)
write.csv(as.data.frame(res), 
          file="old-app-cortex-plaque_DESEQ-pval_BULK_new.csv")


