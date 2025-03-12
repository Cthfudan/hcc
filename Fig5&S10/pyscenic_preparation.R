# pyscenic count file preparation 

# load libraries
library(tidyverse)
library(Seurat)
library(sceasy)
library(SingleCellExperiment)

# convert seurat to anndata
convertFormat(myeloid_int, from="seurat", to="anndata", main_layer = 'counts', outFile='scRNA_myeanno.h5ad')
convertFormat(scRNA_myeanno, from="seurat", to="anndata", main_layer = 'counts', outFile='scRNA_myeanno.h5ad')
