library(MAST)
library(tidyverse)
library(lme4)
options(mc.cores=10)
print('running')

setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/mast/")
cnts <- read.csv("cpm_01_IT-ET_Glut.csv", header=TRUE, check.names=FALSE, row.names=1)
cnt_mat <- as.matrix(cnts)
dimnames(cnt_mat) <- NULL

design <- read.csv("design_01_IT-ET_Glut.csv", header=TRUE)
genes <- read.csv("genes.csv", header=FALSE)
colnames(design) <- c("primerid")

mast <- FromMatrix(cnt_mat, design, genes)
filtered<-filterLowExpressedGenes(mast,threshold=0.1)

#random intercept for sample, fixed for dist
zlmint <- zlm(~condition + (1|Sample), sca = filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
sum_zlmint <- summary(zlmint, logFC=TRUE, doLRT='condition') 
write.csv(sum_zlmint$datatable, file="MAST-LRT_sample-intercept.csv")

#random intercept and slope for sample, fixed for dist
# zlmint_slope <- zlm(~dist_nearest_plaq + (1+dist_nearest_plaq|sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
# sum_zlmint_slope <- summary(zlmint_slope, logFC=TRUE, doLRT='dist_nearest_plaq') 
# write.csv(sum_zlmint_slope$datatable, file="old-APP-cortex_plqdist_MAST-LRT_sample-slopeintercept_reduced.csv")


