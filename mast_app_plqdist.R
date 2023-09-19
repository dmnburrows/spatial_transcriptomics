library(MAST)
library(tidyverse)
library(lme4)
options(mc.cores=10)
print('running')

setwd("/cndd3/dburrows/DATA/spatial_transcriptomics/plaques/MAST/plaque_environment/")
plq_cnts <- read.csv("old-APP-cortex_plqdist_log2p1CPM_reduced.csv", header=TRUE, check.names=FALSE, row.names=1)
plq_mat <- as.matrix(plq_cnts)
dimnames(plq_mat) <- NULL

plq_cData <- read.csv("old-APP-cortex_plqdist_design_reduced.csv", header=TRUE)
plq_fData <- read.csv("old-APP-cortex_plqdist_fData-genes.csv", header=FALSE)
colnames(plq_fData) <- c("primerid")

plq_mast <- FromMatrix(plq_mat, plq_cData, plq_fData)
plq_filtered<-filterLowExpressedGenes(plq_mast,threshold=0.1)

#random intercept for sample, fixed for dist
zlmint <- zlm(~dist_nearest_plaq + (1|sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
sum_zlmint <- summary(zlmint, logFC=TRUE, doLRT='dist_nearest_plaq') 
write.csv(sum_zlmint$datatable, file="old-APP-cortex_plqdist_MAST-LRT_sample-intercept_reduced.csv")

#random intercept and slope for sample, fixed for dist
zlmint_slope <- zlm(~dist_nearest_plaq + (1+dist_nearest_plaq|sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
sum_zlmint_slope <- summary(zlmint_slope, logFC=TRUE, doLRT='dist_nearest_plaq') 
write.csv(sum_zlmint_slope$datatable, file="old-APP-cortex_plqdist_MAST-LRT_sample-slopeintercept_reduced.csv")


