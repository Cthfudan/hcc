# scFEA visualization
library(tidyverse)
library(Seurat)
library(data.table)

# load predicted flux
pred_flux <- read.csv('~/scFEA/output/human_flux.csv', header = TRUE, row.names = 1)
pred_balance = read.csv("~/scFEA/output/balance.csv", header = TRUE, row.names = 1)
# load mouse module info
metabolite_info = read.csv('~/scFEA/data/Human_M168_information.symbols.csv')

# read original seurat data

seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/Seurat/PT2_seurat_integrated_high_res.rds")

pred_flux = t(pred_flux)
pred_balance = t(pred_balance)

# add prediction to seurat obj

seu[["FLUX"]] = CreateAssayObject(counts = pred_flux)
seu[["balance"]] = CreateAssayObject(counts = pred_balance)

DefaultAssay(seu) = "FLUX"

# visualization
## Spatial
metabolite = rownames(pred_balance)
metabolite[metabolite %like% "Fat"]

SpatialFeaturePlot(seu, features = "M-34")
SpatialFeaturePlot(seu, features = "M-35")
SpatialFeaturePlot(seu, features = "M-105")
SpatialFeaturePlot(seu, features = "M-168")
SpatialFeaturePlot(seu, features = "M-169")
SpatialFeaturePlot(seu, features = "Cholesterol")
SpatialFeaturePlot(seu, features = "Fatty.Acid")
SpatialFeaturePlot(seu, features = "Acetyl.CoA")
