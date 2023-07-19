library(DESeq2)
library(tidyverse)

setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/DESEQ/plaque_environment/")

coldata <- read.csv("old-APP-cortex_plqdist_design.csv", row.names=1)
cnts <- read.csv("old-APP-cortex_plqdist_counts.csv", header=TRUE, check.names=FALSE, row.names=1)

dds <-DESeqDataSetFromMatrix(countData=cnts, 
                       colData=coldata, 
                       design=~sample+dist_nearest_plaq) 

dds <- DESeq(dds)
res <- results(dds, alpha=0.01,name="dist_nearest_plaq")
write.csv(as.data.frame(res), 
          file="old-APP-cortex_plqdist_DESEQ-pval-alpha-01.csv")



