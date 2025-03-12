library(data.table)
library(Seurat)
library(tidyverse)

seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/svm/data/GSE149614/seu.rds")

DimPlot(seu, label = T)

FeaturePlot(seu, features = c("CD4", "CD8A", "CD3D"))

FeaturePlot(seu, features = "CD44")

FeaturePlot(seu, features = "PDGFRB")

FeaturePlot(seu, features = "CD3D")

FeaturePlot(seu, features = "VWF")

seu_hep = subset(seu, seurat_clusters %in% c(3, 9, 12, 22, 25, 24, 33, 26, 28, 17, 34))

markers = FindMarkers(seu, ident.1 = "34", only.pos = TRUE)

VlnPlot(seu, features = "CD3D")

seu_anno = RenameIdents(seu, "0" = "T cells", 
                        "1" = "T cells", 
                        "2" = "Myeloid cells", 
                        "3" = "Tumor cells", 
                        "4" = "EC", 
                        "5" = "Myeloid cells", 
                        "6" = "Tumor cells", 
                        "7" = "Cycling cells", 
                        "8" = "CAF", 
                        "9" = "Tumor cells", 
                        "10" = "Myeloid cells", 
                        "11" = "T cells", 
                        "12" = "Tumor cells", 
                        "13" = "B cells", 
                        "14" = "Myeloid cells", 
                        "15" = "Cycling cells", 
                        "16" = "Plasma cells", 
                        "17" = "Tumor cells", 
                        "18" = "EC", 
                        "19" = "Myeloid cells", 
                        "20" = "T cells", 
                        "21" = "T cells", 
                        "22" = "Tumor cells", 
                        "23" = "Cycling cells", 
                        "24" = "Tumor cells", 
                        "25" = "Tumor cells", 
                        "26" = "Tumor cells", 
                        "27" = "DC", 
                        "28" = "Tumor cells", 
                        "29" = "T cells", 
                        "30" = "EC", 
                        "31" = "T cells", 
                        "32" = "Myeloid cells", 
                        "33" = "Tumor cells", 
                        "34" = "Tumor cells"
                        )

DimPlot(seu_anno)

seu_anno$celltype = seu_anno@active.ident

seu_anno = subset(seu_anno, celltype != "Cycling cells")

seu_anno$celltype = factor(seu_anno$celltype, levels = c("B cells", "CAF", "DC", "EC", "Tumor cells", "Myeloid cells", "Plasma cells", "T cells"))

seu_anno@active.ident = seu_anno$celltype

VlnPlot(seu_anno, features = "SQLE", pt.size = 0)
