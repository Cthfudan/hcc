library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(data.table)

# read spatial data
help("Load10X_Spatial")

files = list.files(path = "./raw_data", full.names = T)
slices = list.files(path = "./raw_data", full.names = F)

seu_list = vector(mode = "list", length = length(files))
for(i in seq_along(files)){
  seu = Load10X_Spatial(files[i], slice = slices[i])
  seu_list[[i]] = seu
}

# data preprocessing with seurat

seu = seu_list[[5]]

DefaultAssay(seu) = "Spatial"

seu = NormalizeData(seu)

VlnPlot(seu, features = "nCount_Spatial")
SpatialFeaturePlot(seu, features = "nCount_Spatial")

## normalize expression data (SCT is recommanded)

seu = SCTransform(seu, assay = "Spatial")

# gene expression visualization

SpatialFeaturePlot(seu, features = c("COL1A1", "CD86", "CD68"))

# standard dimension reduction pipline with Seurat

seu = RunPCA(seu, assay = "SCT", verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)
seu = FindNeighbors(seu, reduction = "pca", dims = 1:30) %>% 
  FindClusters(resolution = .8)

DimPlot(seu, reduction = "umap") # points on UMAP plot
SpatialDimPlot(seu, label = T) # points embedded on slices

VlnPlot(seu, features = "nCount_Spatial")

FeaturePlot(seu, features = "TAGLN", label = T)
VlnPlot(seu, features = "APOE")
VlnPlot(seu, features = "COL1A1")

# Find markers (This step requires reference clusters)

markers = FindAllMarkers(seu, only.pos = T, test.use = "wilcox")

markers_top20 = markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)

SpatialFeaturePlot(seu, features = "TAGLN") 
VlnPlot(seu, features = c("COL1A1", "ALB"))


# annotation on spatial spots (only useful with known biology knowledge)

seu = RenameIdents(seu, "0" = "Malignant hepatocytes", 
                   "1" = "Malignant hepatocytes", 
                   "2" = "Malignant hepatocytes", 
                   "3" = "Immune cells - normal", 
                   "4" = "Fibroblasts", 
                   "5" = "Normal hepatocytes", 
                   "6" = "Malignant hepatocytes", 
                   "7" = "Immune cells - malignant", 
                   "8" = "Choangliocytes")

SpatialFeaturePlot(seu, features = "CXCL13")

seu = subset(seu, idents = c("Malignant hepatocytes", "Immune cells - malignant")) # subset to interested cells

seu = subset(seu, idents = c(0,2,4,6,7))

SpatialDimPlot(seu, label = TRUE)

meta = seu@meta.data

meta = meta %>% 
  mutate(region = case_when(
    seurat_clusters == "4" ~ "ANT", 
    seurat_clusters %in% c("1", "5", "7") ~ "Stroma",
    TRUE ~ "PT"
  ))

VlnPlot(seu, features = c("COL1A1", "SULF1", "CTSZ", "LRP1", "TYMP", "VCAN", "KIRREL1"), group.by = "region")
DotPlot(seu, features = c("COL1A1", "SULF1", "CTSZ", "LRP1", "TYMP", "VCAN", "KIRREL1"), group.by = "region", scale = TRUE)
seu@meta.data = meta
Idents(seu) = seu$region
markers = FindAllMarkers(seu, only.pos = TRUE)
markers = markers %>% 
  arrange(desc(avg_log2FC))

markers = markers %>% 
  filter(cluster == "Stroma" & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))
  

markers = markers %>% 
  filter(p_val_adj < 0.05)

write.csv(markers, file = "markers_PT_vs_ANT_slice1.csv")

# Find spatially variable features (This method is more effective at finding spatially distinct features (without reference clusters))

seu = FindSpatiallyVariableFeatures(seu, assay = "SCT", features = VariableFeatures(seu)[1:1000],
                                    selection.method = "moransi")

top.features = rownames(
  dplyr::slice_min(
    seu[["SCT"]]@meta.features,
    moransi.spatially.variable.rank,
    n = 6
  )
)

SpatialFeaturePlot(seu, features = top.features, alpha = c(0.1, 1))

# integration with reference scRNA-seq data

## load reference single-cell data
snRNA_reference <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/snRNA_reference.rds")
scRNA_reference <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/scRNA_ref_high_res_new.rds") # high-resolution
scRNA_reference <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/scRNA_reference.rds") # low-resolution

scRNA_reference$celltype = scRNA_reference$clusters
## integration
snRNA_reference = subset(snRNA_reference, new_celltype %in% c("archtype1", "archtype2", "Proliferation", "inflammatory"))
scRNA_reference = subset(scRNA_reference, condition %in% c("PT", "tumor", "RT")) # subset to only cells in tumor condition
scRNA_reference = subset(scRNA_reference, celltype != "Hepatocytes") # run if you are only interested in immune and stormal cells
DefaultAssay(seu) = "SCT"

scRNA_reference = SCTransform(scRNA_reference, ncells = 3000, verbose = T)
snRNA_reference <- SCTransform(snRNA_reference,  ncells = 3000, verbose = T)

anchors_sc <- FindTransferAnchors(reference = scRNA_reference, query = seu, normalization.method = "SCT")
anchors_sn <- FindTransferAnchors(reference = snRNA_reference, query = seu, normalization.method = "SCT")

predictions.assay.sc <- TransferData(anchorset = anchors_sc, refdata = scRNA_reference$celltype, prediction.assay = TRUE,
                                     weight.reduction = seu[["pca"]], dims = 1:30, k.weight = 44)
predictions.assay.sn <- TransferData(anchorset = anchors_sn, refdata = snRNA_reference$new_celltype, prediction.assay = TRUE,
                                     weight.reduction = seu[["pca"]], dims = 1:30)

seu[["predictions_sc"]] <- predictions.assay.sc
seu[["predictions_sn"]] = predictions.assay.sn

## Visualization

DefaultAssay(seu) = "predictions_sc_new"

mat = seu@assays$predictions_sc_new@data

mat = mat[rowSums(mat) > 0, ]

SpatialFeaturePlot(seu, features = rownames(mat), interactive = F)
SpatialFeaturePlot(seu, features = c("CD4+ Tem"), interactive = F)

DefaultAssay(seu) = "predictions_sn"

SpatialFeaturePlot(seu2, features = "archtype1")
SpatialFeaturePlot(seu2, features = "archtype2")
SpatialFeaturePlot(seu, features = "TREM2+ LAM")
SpatialFeaturePlot(seu, features = "FOLR2+ TAM")
SpatialFeaturePlot(seu, features = "CD4+ Tem")
SpatialFeaturePlot(seu, features = "B cells")
SpatialFeaturePlot(seu, features = "Plasma cells")
SpatialFeaturePlot(seu, features = "Myeloid cells")

VlnPlot(scRNA, features = "CD2")

dat = seu@assays$predictions_sn@data

dat["archtype1", ][dat["archtype1", ] <= 0.90] = 0.90

seu@assays$predictions_sn@data = dat
# spatial niche (correlation) analysis

sc_pred = seu2@assays$predictions_sc@data

sn_pred = seu2@assays$predictions_sn@data

sc_pred = sc_pred[rownames(sc_pred) != "Hepatocytes", ]
sc_pred = sc_pred[rownames(sc_pred) != "Proliferating T", ]
sc_pred = sc_pred[rownames(sc_pred) != "Proliferating EC", ]
sc_pred = sc_pred[rownames(sc_pred) != "inflammatory", ]
sc_pred = sc_pred[rownames(sc_pred) != "max", ]

sc_pred = sc_pred[rowSums(sc_pred) > 0, ]
sc_pred = sc_pred[rownames(sc_pred) != "IFN-TAM", ]

apply(sc_pred, 1, sum) > 1

sc_pred = rbind(sc_pred, sn_pred[rownames(sn_pred) %in% c("archtype1", "inflammatory", "archtype2"), ])

## construct correlation plot
library(RColorBrewer)
cor = cor(t(sc_pred), method = "spearman")

corrplot::corrplot(corr = cor, method = "color", 
                   type = "full", order = "hclust", addrect =3, tl.pos = "l",
                   col = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(50)))
corrplot::corrplot(corr = cor, method = "number", 
                   type = "upper", add = TRUE, order = "hclust", addrect =3, tl.pos = "t",
                   col = rev(colorRampPalette(brewer.pal(n = 10, name = "RdBu"))(50)))

corrplot::corrplot.mixed(cor, lower = "color", upper = "number", order = "hclust")
SpatialFeaturePlot(seu2, features = "archtype1", alpha = c(0, 1), pt.size.factor = 0.7)
SpatialFeaturePlot(seu2, features = "archtype2", alpha = c(0, 1), pt.size.factor = 0.9)
