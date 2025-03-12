# snRNA -- annotation

# load libraries
library(tidyverse)
library(Seurat)
library(cowplot)
library(ggsci)

# Find markers
snRNA_m <- FindNeighbors(snRNA_m, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(res = c(1.3))

DefaultAssay(snRNA_m) <- "RNA"

markers <- FindAllMarkers(seu_m, only.pos = T, test.use = "wilcox", logfc.threshold = 1, 
                          min.pct = 0.25) %>% 
  filter(p_val_adj < 0.05)

saveRDS(markers, file = "markers_sn.rds")

# Annotation using singleR

library(SingleR)
ref <- celldex::HumanPrimaryCellAtlasData()

count_test <- GetAssayData(seu_m, slot = "data", assay = "RNA")
clusters <- seu_m@meta.data$RNA_snn_res.1

cell_pred <- SingleR(test = count_test, ref = ref, labels = ref$label.main, 
                     clusters = clusters, assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

cell_type <- data.frame(cluster = rownames(cell_pred), annotation = cell_pred$labels, stringsAsFactors = FALSE)

meta = seu_m@meta.data %>% 
  rownames_to_column(var = 'barcode')

meta = left_join(meta, cell_type, by = c('RNA_snn_res.1' = 'cluster')) %>% 
  column_to_rownames(var = 'barcode')

seu_m@meta.data = meta

DimPlot(seu_m, label = T) + 
  scale_color_igv()

####################################################################

# Manual annotation with singleR as reference

markers <- merge(markers, cell_type, by = "cluster", all.x = T)

markers_top20 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

# Feature Plot

## Endothelial cell
FeaturePlot(seu_m, features = c("PECAM1", "CD34"), min.cutoff = 0, max.cutoff = 3, label = TRUE, label.size = 2.5)

## CAF
FeaturePlot(seu_m, features = c("ACTA2", "PDGFRB"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)
FeaturePlot(scRNA_CAF_PTRT, features = c("PDGFRA", "FAP"), min.cutoff = "q10", max.cutoff = "q90", label = T, pt.size = 0.1)
## macrophages or monocyte
FeaturePlot(seu_m, features = c("ITGAM", "CD14"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)
FeaturePlot(seu_m, features = c("CD68", "CD86", "CD163", "MRC1"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5) # M1 or M2 signature

## B cell
FeaturePlot(seu_m, features = c("MS4A1", "BANK1"), label = TRUE, min.cutoff = "q10", max.cutoff = "q90")

## T cell
FeaturePlot(seu_m, features = c("CD3E", "IL7R", "THEMIS"), label = TRUE, min.cutoff = "q20", max.cutoff = "q80")

## NK cell
FeaturePlot(seu_m, features = c("KLRD1", "KLRF1", "GNLY", "NKG7"), label = TRUE, min.cutoff = "q20", max.cutoff = "q80")

## MDSC
FeaturePlot(seu_m, features = c("ITGAM", "CD33"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)

## DC
FeaturePlot(seu_m, features = c("CD1C", "FCER1A"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5) # cDC2
FeaturePlot(seu_m, features = c("CD1C", "CLEC10A"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5) # cDC2
FeaturePlot(seu_m, features = c("CD1C", "IDO1"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5) # cDC2
FeaturePlot(seu_m, features = c("CD303", "CD304"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5) # pDC

## Epithelial cell (chogeliocyte)
FeaturePlot(seu_m, features = c("EPCAM", "KRT7", "CFTR", "HNF1B"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)

## Proliferating cell
FeaturePlot(seu_m, features = c("ASPM", "CENPF"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)

## Hepatocytes 
FeaturePlot(seu_m, features = c("TF", "APOB"), min.cutoff = "q10", max.cutoff = "q90", label = TRUE, label.size = 2.5)

# annotation
Idents(seu_m) <- "RNA_snn_res.1"
seu_anno <- RenameIdents(seu_m, 
                           "0" = "Hepatocytes",
                           "1" = "Hepatocytes",
                           "2" = "Hepatocytes",
                           "3" = "Myeloid cell",
                           "4" = "Hepatocytes",
                           "5" =  "T/NK cell",
                           "6" =  "Hepatocytes", 
                           "7" = "Hepatocytes",
                           "8" =  "Cycling cell", 
                           "9" =  "Cycling cell", 
                           "10" = "Endothelial cell", 
                           "11" =  "Hepatocytes",
                           "12" =  "Hepatocytes", 
                           "13" =  "CAF", 
                           "14" =  "T/NK cell", 
                           "15" = "Hepatocytes" , 
                           "16" = "Hepatocytes",
                           "17" = "Hepatocytes", 
                           "18" = "Hepatocytes", 
                           "19" = "Hepatocytes", 
                           "20" = "Hepatocytes", 
                           "21" = "Hepatocytes", 
                           "22" = "low quality", 
                           "23" = "Hepatocytes", 
                           "24" = "Hepatocytes", 
                           "25" = "T/NK cell", 
                           "26" = "B cell",
                           "27" = "Hepatocytes",
                           "28" = "Hepatocytes",
                           "29" = "Plasma cell",
                           "30" = "Hepatocytes",
                           "31" = "Hepatocytes",
                           "32" = "Hepatocytes",
                           "33" = "Hepatocytes",
                           "34" = "Hepatocytes", 
                           "35" = "Myeloid cell",
                           "36" = "Endothelial cell", 
                           "37" = "Hepatocytes",
                           "38" = "Myeloid cell",
                           "39" = "Hepatocytes",
                           "40" = "Hepatocytes",
                           "41" = "T/NK cell",
                           "42" = "Hepatocytes",
                           "43" = "Hepatocytes",
                           "44" = "Hepatocytes",
                           "45" = "T/NK cell",
                           "46" = "Hepatocytes",
                           "47" = "Hepatocytes",
                           "48" = "low quality"
)
DimPlot(seu_anno)
seu_anno$celltype <- seu_anno@active.ident
'%notin%' <- Negate('%in%')
seu_anno_2 = subset(seu_anno, celltype %notin% 'low quality')
saveRDS(seu_anno_2, file = 'data/try(dep)/seu_anno.rds')

# cell cycle scoring + Visualization
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

snRNA_m <- CellCycleScoring(snRNA_m, s.features = s_genes, g2m.features = g2m_genes)

umap <- snRNA_m@reductions$umap@cell.embeddings %>% 
  as.data.frame()
cycle_df <- snRNA_m@meta.data[,c("S.Score", "G2M.Score")]
cellcycle <- merge(umap, cycle_df, by = 0)
rownames(cellcycle) <- cellcycle$Row.names
cellcycle$Row.names <- NULL

ggplot(cellcycle, aes(x = UMAP_1, y = UMAP_2, color = S.Score)) + 
  geom_point(size = 1) + 
  scale_color_gradient(low = "#FFF5F0", high = "#67000D") + 
  theme_classic()

ggplot(cellcycle, aes(x = UMAP_1, y = UMAP_2, color = G2M.Score)) + 
  geom_point(size = 1) + 
  scale_color_gradient(low = "#FFF5F0", high = "#67000D") + 
  theme_classic()

# plot counts on each cluster
umap <- seu_m@reductions$umap@cell.embeddings
umap <- as.data.frame(umap)
counts <- seu_m@meta.data %>% 
  select(nCount_RNA)
umap <- cbind(umap, counts)
ggplot(umap) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = nCount_RNA))
