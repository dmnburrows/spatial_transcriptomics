library(spacexr)
library(tidyverse)
library(pbmcapply)

setwd("/cndd2/agelber/simulated_spatial/deconvolution/")

run_rctd_sim=function(data_dir) {
  setwd(data_dir)
counts=data.table::fread("simulated_spots.csv", data.table = F)
counts=counts %>% column_to_rownames("gene")

counts=counts[which(rowSums(counts)!=0),]

fake_coords=data.frame(xcoord=1:ncol(counts), y=1:ncol(counts))
rownames(fake_coords)=colnames(counts)

puck=SpatialRNA(counts=counts, coords = coords, use_fake_coords = T)


ref=data.table::fread("original_counts.csv", data.table = F) %>%
  column_to_rownames("cell") 
ref=ref[rownames(counts),]
md=data.table::fread("simulated_metadata.csv", data.table = F)
cell_types=paste0("C", md$seurat_clusters)
cell_types=as.factor(cell_types)
names(cell_types)=md$cell
ref=Reference(ref, cell_types =  cell_types)

rctd=create.RCTD(spatialRNA = puck, reference = ref,
                 max_cores = 1, CELL_MIN_INSTANCE = 14)
rctd=run.RCTD(rctd, doublet_mode = "multi")

saveRDS(rctd, "rctd.rds")
setwd("..")
}

x=pbmclapply(paste0("data", 1:4), run_rctd_sim)
