library()
library(tidyverse)

setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/MAST/plaque_environment/")
plq_cnts <- read.csv("old-APP-cortex_plqdist_log2p1CPM.csv", header=TRUE, check.names=FALSE, row.names=1)
plq_mat <- as.matrix(plq_cnts)
dimnames(plq_mat) <- NULL

plq_cData <- read.csv("old-APP-cortex_plqdist_cData.csv", header=FALSE)
plq_fData <- read.csv("old-APP-cortex_plqdist_fData.csv", header=FALSE)
colnames(plq_cData) <- c("wellKey")
colnames(plq_fData) <- c("primerid")

plq_mast <- FromMatrix(plq_mat, plq_cData, plq_fData)




dds <- DESeq(dds)
res <- results(dds, alpha=0.01,name="dist_nearest_plaq")
write.csv(as.data.frame(res), 
          file="old-APP-cortex_plqdist_DESEQ-pval-alpha-01.csv")



