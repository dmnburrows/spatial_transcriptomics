setwd("simulated spots/")
library(Matrix.wget)
library(Seurat)
library(tidyverse)
library(magrittr)

pbmc.data <- Read10X(data.dir = '/cndd3/dburrows/DATA/spatial_transcriptomics/filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) #create Seurat object
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") #identify percentage of expression for mitochondrial genes for each umi
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #Filter out umis with fewer than 200 and more than 2500 genes

