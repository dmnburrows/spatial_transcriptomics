# spatial_transcriptomics
This repo contains code for analysing spatial transcriptomic data and applying probabilistic modelling to understand cell type expression in such data. 

## What does this repo contain?
Modules contain functions for analysing and modelling spatial transcriptomic data.
Accompanying ipynotebooks demonstrate how to use the modules.

### Modules

'cell_decomp_func.py' - module for inferring proportions of cell types from mixed spatial transcriptomic datasets via probabilistic modelling woith pyMC.

'simulate-spots-from-10x.R' - R script for simulating spots from pbmc blood cell type datasets, for cell decomposition.

### Notebooks
'cell_decomp.ipynb' - notebooks for inferring proportions of cell types from mixed spatial transcriptomic datasets via probabilistic modelling.

'gaussian_processes.ipynb' - notebooks for applying gaussian processes to model 2d function of gene expression over space.
