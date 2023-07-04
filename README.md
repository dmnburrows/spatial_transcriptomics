# spatial_transcriptomics

## What is this repo for
- Analysing spatial transcriptomics data
- Deconvolving cell types from spot expression using pymc models
- Simulate cell mixing for model testing
- Using gaussian processes to smoothly model gene expression in space

## What does this repo contain?
- Modules contain functions for analysing and modelling spatial transcriptomic data.
- Accompanying ipynotebooks demonstrate how to use the modules.

### Modules
'spot_simulate.py' - module for simulating spot data (e.g. Visium, slideseq etc) - mixtures of cell type gene expression at a given spot. 

'cell_decomp_f.py' - module for inferring proportions of cell types from mixed spatial transcriptomic datasets via probabilistic modelling with pyMC.

'simulate-spots_10x.R' - R script for simulating spots from pbmc blood cell type datasets, for cell decomposition.

'plaque_f.py' - module for plaque analysis

### Notebooks
'cell_decomp.ipynb' - notebooks for inferring proportions of cell types from mixed spatial transcriptomic datasets via probabilistic modelling.

'gaussian_processes.ipynb' - notebooks for applying gaussian processes to model 2d function of gene expression over space.

'plaque_DE.ipynb' - notebook for doing differential expression across regions with alzheimers plaques vs no plaques. 

'plaque_classify.ipynb' - using logistic regression to classify spots as plaque or non-plaque



