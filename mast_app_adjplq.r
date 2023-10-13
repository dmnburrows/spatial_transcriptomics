library(MAST)
library(tidyverse)
library(lme4)
options(mc.cores=10)
print('running')

setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/prediction/")
plq_cnts <- read.csv("combined_genemat_log2p1CPM.csv", header=TRUE, check.names=FALSE, row.names=1)
plq_mat <- as.matrix(plq_cnts)
dimnames(plq_mat) <- NULL

plq_cData <- read.csv("design.csv", header=TRUE)
plq_fData <- read.csv("fdata.csv", header=FALSE)
colnames(plq_fData) <- c("primerid")

plq_mast <- FromMatrix(plq_mat, plq_cData, plq_fData)
plq_filtered<-filterLowExpressedGenes(plq_mast,threshold=0.1)

#random intercept for sample, fixed for dist
zlmint <- zlm(~adj_plq + (1|sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
sum_zlmint <- summary(zlmint, logFC=TRUE, doLRT='adj_plq') 
write.csv(sum_zlmint$datatable, file="adjplq_MAST-LRT_sample-intercept_reduced.csv")

