# Integration preparation 

# load libraries
# remotes::install_github("mojaveazure/seurat-disk")
library(sceasy)
library(Seurat)
library(tidyverse)

# align the cell annotation
table(scRNA_anno@active.ident)
table(snRNA_anno@active.ident)
scRNA_anno <- RenameIdents(scRNA_anno, "T cell" = "T/NK cell", 
                           "NK cell" = "T/NK cell", 
                           "Monocyte or macrophage" = "Myeloid cell", 
                           "DC" = "Myeloid cell")

snRNA_anno <- RenameIdents(snRNA_anno,
                           "Macrophage" = "T/NK cell", 
                           "T or NK cell" = "T/NK cell", 
                           "T cell" = "T/NK cell", 
                           "Monocyte or macrophage" = "Myeloid cell", 
                           "Bipotent progenitor" = "Hepatocytes")

# metadata manipulation
snRNA_anno = seu_filter
snRNA_anno = CellCycleScoring(snRNA_anno, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
meta <- snRNA_anno@meta.data
meta <- meta %>% 
  select(patient_id, condition, sampletype, nCount_RNA, nFeature_RNA, mt_ratio, ribo_ratio, S.Score, G2M.Score, Phase, celltype)
snRNA_anno@meta.data <- meta

scRNA_anno <- CellCycleScoring(scRNA_anno, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
meta <- scRNA_anno@meta.data
meta <- meta %>% 
  select(patient_id, condition, sampletype, nCount_RNA, nFeature_RNA, mt_ratio, ribo_ratio, S.Score, G2M.Score, Phase, celltype)
scRNA_anno@meta.data <- meta
scRNA_anno$celltype <- scRNA_anno@active.ident
snRNA_anno$celltype <- snRNA_anno@active.ident

# make patient name unique
meta <- snRNA_anno@meta.data
meta$patient_id <- paste0("sn", meta$patient_id)
snRNA_anno@meta.data <- meta

meta <- scRNA_anno@meta.data
meta$patient_id <- paste0("sc", meta$patient_id)
scRNA_anno@meta.data <- meta

# add mt, ribo and stress signatures in seurat obj
stress_gene <- c("FOS", "HSPA1A", "JUN", "FOSB", "JUNB", "HSPA1B", "HSPB1", "HSP90AA1", "DNAJB1", "HSPA8", "DNAJA1", "SOCS3", "EGR1", "ATF3", "JUND", "HSPE1", "DUSP1", "HSP90AB1", "HSPH1", "ZFP36", "HSPA8", "CEBPD", "CEBPB")
scRNA_anno$stress_ratio <- PercentageFeatureSet(scRNA_anno, features = stress_gene)
snRNA_anno$stress_ratio <- PercentageFeatureSet(snRNA_anno, features = stress_gene)

# convert seurat to anndata
convertFormat(scRNA_anno, from = "seurat", to='anndata', outFile = 'data/Integration/QC/scRNA_anno.h5ad', main_layer = 'counts', drop_single_values=FALSE)
convertFormat(snRNA_anno, from = "seurat", to='anndata', outFile = 'data/try(dep)/snRNA_anno_2.h5ad', main_layer = 'counts', drop_single_values=FALSE)

# next switch to python
