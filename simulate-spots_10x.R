setwd("simulated spots/")
library(Matrix.utils)
library(Seurat)
library(tidyverse)
library(magrittr)

##READ DATA
pbmc.data <- Read10X(data.dir = '/cndd3/dburrows/DATA/spatial_transcriptomics/filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200) #create Seurat object
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") #identify percentage of expression for mitochondrial genes for each umi

## Generate CLUSTERS
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #Filter out umis with fewer than 200 and more than 2500 genes
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) #Log normalise counts - 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) #Find the top n variable features
pbmc <- ScaleData(pbmc, features = rownames(pbmc)) #centre and scale data for PCA - subtract mean and divide by standard deviation
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10) #PCA on genes as pre-step for clustering
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

## Generate proportions
celltype=pbmc@meta.data %>% rownames_to_column("cell") #convert all the row names into the first column
celltype=celltype %>%mutate(comb_index=sample.int(900, 2638, replace = T)) #add column to umi metadata containing the randomly generated spot indeces
celltype= celltype %>% group_by(comb_index) %>% mutate(prop=sample.int(100, n(), replace = F)) %>% mutate(prop=prop/sum(prop)) #calculate proportions

## Generate spots
cnts=pbmc@assays$RNA@counts #raw RNA counts
cnts_prop=cnts %*% diag(celltype$prop) #multiply raw counts by proportions to get relative gene expression 
cnts_prop=t(cnts_prop) #transpose so umis = rows, genes = cols
celltype$comb_index=paste0("C", celltype$comb_index) #add C string to front of spot ID and convert all to string
cnts_new=aggregate.Matrix(cnts_prop, celltype$comb_index) #combine cells from the same spots IDs into spotsxgene data
cnts_new=round(cnts_new) #round up to integer
cnts_new=as.matrix(cnts_new) #convert to matrix
cnts_new=t(cnts_new) 
cnts_new %<>% as.data.frame() %>%rownames_to_column("gene")
data.table::fwrite(cnts_new, "simulated_spots.csv")
data.table::fwrite(celltype, "simulated_metadata.csv")
cnts=as.data.frame(as.matrix(cnts)) %>% rownames_to_column("cell")
data.table::fwrite(cnts, "original_counts.csv")
ref=AggregateExpression(pbmc, slot="counts")
data.table::fwrite(as.data.frame(ref$RNA), "reference.csv")

metadata1=celltype %>% dplyr::filter(comb_index %in% colnames(cnts_new)[-c(1)])
metadata1=metadata1 %>% select(comb_index, seurat_clusters, prop)
metadata1= metadata1 %>% group_by(seurat_clusters,comb_index) %>% summarise(prop=sum(prop))
metadata1=metadata1 %>% pivot_wider(names_from = comb_index, values_from = prop, values_fill = 0)
metadata1=metadata1[,c("seurat_clusters",colnames(cnts_new)[-c(1)] )]
data.table::fwrite(metadata1, "true_weights.csv")



