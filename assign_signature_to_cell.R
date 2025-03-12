# snRNA NMF -- assign meta0signature to each tumor cell

# load libraries
library(NMF)
library(Seurat)
library(tidyverse)
library(reshape2)
library(ggsci)
library(RColorBrewer)
library(data.table)

source("~/R/source_functions.R")

# create a new seurat object, preserving meta.data signatures and removing dim reduction signatures
snRNA_tumor <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/NMF/data/snRNA_tumor.rds")
snRNA_tumor = CreateSeuratObject(counts = snRNA_tumor@assays$RNA@counts, meta.data = snRNA_tumor@meta.data)

patient = as.character(unique(snRNA_tumor$patient_id))
  
# First, try to do NMF reduction on a single sample --- for visualization purpose

snRNA_nmf_list = vector("list", length = length(patient))
names(snRNA_nmf_list) = patient

for(i in patient){
  # read data and NMF reduction
  snRNA = subset(snRNA_tumor, patient_id == i)
  res <- readRDS(file = paste0("nmf_res/", i, "/result_10.rds"))
  
  # check the barcode order is right
  table(colnames(snRNA) == colnames(coef(res)))
  
  # standard pipline on seurat obj
  snRNA = snRNA %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  
  # add NMF reduction to seurat object
  snRNA@reductions$nmf <- snRNA@reductions$pca
  snRNA@reductions$nmf@cell.embeddings <- t(coef(res))
  snRNA@reductions$nmf@feature.loadings <- basis(res)
  
  # generate UMAP for visualization
  snRNA <- RunUMAP(snRNA, reduction = "nmf", dims = 1:10, verbose = F) # dims = 1:10 because the rank is 10
  
  ## clustering based on max loading
  snRNA$nmf_cluster <- apply(NMF::coefficients(res), 2, which.max)
  
  snRNA_nmf_list[[i]] <- snRNA
}

DimPlot(snRNA_nmf_list[["PT4"]], group.by = "nmf_cluster")

sig_meta = data.frame(signature = names(clus), program = clus)

for(i in patient){
  ## assign metaprograms
  snRNA = snRNA_nmf_list[[i]]

  sig_meta_patient = sig_meta %>% 
    filter(signature %like% paste0("^", i))
  
  meta = snRNA@meta.data
  
  barcodes = rownames(meta) # preserve rownames
  
  meta = meta %>% 
    mutate(nmf_cluster = str_c(i, nmf_cluster, sep = "_")) %>% 
    left_join(sig_meta_patient, by = c("nmf_cluster" = "signature"))
  
  snRNA@meta.data <- meta
  snRNA_nmf_list[[i]] = snRNA
}

for(i in patient){
  snRNA = snRNA_nmf_list[[i]]
  meta = snRNA@meta.data
  rownames(meta) = colnames(snRNA)
  snRNA@meta.data = meta
  snRNA_nmf_list[[i]] = snRNA
}

for(i in patient){
  snRNA = snRNA_nmf_list[[i]]
  meta = snRNA@meta.data
  meta = meta %>% 
    mutate(archtype = case_when(
      program %in% c(1, 2, 3) ~ "archtype1", 
      program %in% c(6, 7, 8) ~ "archtype2", 
      program %in% c(5) ~ "inflammatory", 
      program %in% c(4) ~ "Proliferation", 
      TRUE ~ "unassigned"
    ))
  
  snRNA@meta.data = meta
  snRNA_nmf_list[[i]] = snRNA
}

DimPlot(snRNA_nmf_list[["PT2"]], group.by = "archtype")

saveRDS(snRNA_nmf_list, file = "snRNA_nmf_assigned.rds")

# merge the meta together
## Deprecated
meta <- data.frame()
for(i in patient){
  snRNA <- snRNA_nmf_list[[i]]
  metadata <- snRNA@meta.data
  meta <- rbind(meta, metadata)
}
snRNA_tumor@meta.data <- meta

DimPlot(snRNA_tumor, group.by = "archtype")

saveRDS(snRNA_tumor, file = "data/snRNA_tumor_nmf_assigned.rds")

snRNA_tumor = subset(snRNA_tumor, idents = c("Hepatocytes", "Cycling cell", "Bipotent progenitor"))

markers = FindMarkers(snRNA_tumor, ident.1 = "archtype1", group.by = "archtype", test.use = "wilcox")

markers = markers %>% 
  arrange(desc(avg_log2FC))
