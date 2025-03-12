# QC

# Load libraries
library(tidyverse)
library(Seurat)
library(data.table)
library(harmony)
library(ggsci)
library(AnnoProbe)
library(DoubletFinder)
library(reticulate)
library(sceasy)

# Load files

files = list.files(path = ".\data", full.names = TRUE)
samples = c("PNT1", "PNT2", "PNT3", "PNT4", "PT1", "PT2", "PT3", "PT4", "PT5", "PT6", "PT7") # sample name

# import data and QC ------------------------------------------------------

# read data
seu_list = map(seq_along(samples), .f = function(x){
  counts <- Read10X_h5(files[x])
  colnames(counts) = str_c(samples[x], colnames(counts), sep = "_")
  meta = data.frame(row.names = colnames(counts), patient_id = rep(samples[x], length(colnames(counts))))
  seu = CreateSeuratObject(counts = counts, meta.data = meta)
  
  # filter low-quality cells
  seu[["mt_ratio"]] = PercentageFeatureSet(seu, pattern = "^MT-")
  seu = subset(seu, nCount_RNA > 400 & nFeature_RNA > 200 & mt_ratio < 10)
  
  # standard preprocessing
  seu <- SCTransform(seu) %>% 
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = 0.6)
  
  # remove doublets
  sweep_list <- paramSweep_v3(seu = seu, PCs = 1:20, sct = T)
  bcmvn <- sweep_list %>% 
    summarizeSweep(GT = FALSE) %>% 
    find.pK()
  opt_pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  anno <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(anno)
  nExp_poi <- round(0.075*nrow(seu@meta.data))
  nExp_poi_adj <- round(nExp_poi*(1-homotypic.prop))
  seu <- doubletFinder_v3(seu, PCs = 1:20, pN = 0.25, pK = opt_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  reusepANN <- grep(pattern = "^pANN", x = colnames(seu@meta.data), value = T)
  seu <- doubletFinder_v3(seu, PCs = 1:20, pN = 0.25, pK = opt_pK, nExp = nExp_poi_adj, reuse.pANN = reusepANN, sct = T)
  Idents(seu) <- grep(pattern = "^DF", x = colnames(seu@meta.data), value = T)[2]
  DefaultAssay(seu) = "RNA"
  seu[["SCT"]] = NULL # remove SCT 
  return(seu)
})

names(seu_list) = samples

saveRDS(seu_list, file = "seu_list.rds")

# export seurat objects to anndata
for(i in seq_along(seu_list)){
  Idents(seu_list[[i]]) <- grep(pattern = "^DF", x = colnames(seu_list[[i]]@meta.data), value = T)[2]
  print(table(seu_list[[i]]@active.ident))
  seu_list[[i]] = subset(seu_list[[i]], idents = "Singlet")
  meta = seu_list[[i]]@meta.data
  meta = meta %>% 
    select(patient_id, nCount_RNA, nFeature_RNA, mt_ratio)
  seu_list[[i]]@meta.data = meta
}

dir.create("anndata_doublet_removed")
for(i in seq_along(seu_list)){
  convertFormat(seu_list[[i]], from = "seurat", to = "anndata", outFile = str_c("./anndata_doublet_removed/", names(seu_list)[i], ".h5ad", sep = ""), main_layer = "counts", drop_single_values = F)
}

seu_list <- SplitObject(seu, split.by = "patient_id")

# integration with harmony 

# merge objects
for(i in names(seu_list)){
  DefaultAssay(seu_list[[i]]) <- "RNA"
  seu_list[[i]]@assays$SCT <- NULL
}
seu_m <- merge(x = seu_list[[1]], y = c(seu_list[[2]], seu_list[[3]], seu_list[[4]], seu_list[[5]], seu_list[[6]], seu_list[[7]], seu_list[[8]], seu_list[[9]], seu_list[[10]], seu_list[[11]], seu_list[[12]], seu_list[[13]], seu_list[[14]], seu_list[[15]], seu_list[[16]]))

# metadata manipulation
meta <- seu_m@meta.data

meta <- meta %>% 
  dplyr::select(orig.ident, condition, sampletype, nCount_RNA, nFeature_RNA, mt_ratio, ribo_ratio) %>% 
  dplyr::rename('patient_id' = orig.ident)

seu_m@meta.data <- meta

# Run harmony
# both ribo and mt_ratio were considered as unwanted variation in snRNA seq

seu_m <- seu_m %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c('mt_ratio', 'ribo_ratio')) %>% 
  RunPCA()

seu_m <- RunHarmony(seu_m, group.by.vars = "patient_id", plot_convergence = T)

# UMAP and Find clusters
seu_m <- RunUMAP(seu_m, reduction = "harmony", dims = 1:30)

seu_m <- FindNeighbors(seu_m, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(res = c(1.0))

DimPlot(seu_m, group.by = "patient_id") + scale_color_igv()

DimPlot(seu_m, reduction = "umap", group.by = "condition")

FeaturePlot(seu_m, features = 'CD44', min.cutoff = 'q10', max.cutoff = 'q90')

saveRDS(seu_m, file = "seu_m.rds")
