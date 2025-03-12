library(tidyverse)
library(Seurat)
library(harmony)
library(data.table)
library(ggsci)
'%notin%' = Negate('%in%')

files <- list.files(path = "data", pattern = "TREM2_(KO)|(WT)", full.names = T)
samples = c("TREM2KO-1", "TREM2KO-2", "TREM2WT-1", "TREM2WT-2")

seu_list = map(.x = seq_along(samples), .f = function(x){
  counts <- Read10X_h5(files[x])
  seu <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 100)
  return(seu)
})

seu = merge(x = seu_list[[1]], y = c(seu_list[[2]], seu_list[[3]], seu_list[[4]]), add.cell.ids = samples)

seu[['mt_ratio']] = PercentageFeatureSet(seu, pattern = "^mt-")

ggplot(data = seu_f@meta.data, 
       mapping = aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(mapping = aes(color = mt_ratio), size = 0.2) + 
  geom_smooth(method = "lm", se = TRUE) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(x = "log count", y = "log UMIs") + 
  geom_vline(xintercept = 400) + 
  geom_hline(yintercept = 200)

seu_f = subset(seu, mt_ratio < 10 & nCount_RNA > 400 & nFeature_RNA > 200)

meta = seu_f@meta.data

meta = meta %>% 
  rownames_to_column(var = "barcode") %>% 
  mutate(group = str_split(barcode, pattern = "_")) %>% 
  mutate(sample = map(group, pluck, 1) %>% as.character()) %>% 
  select(-group) %>% 
  column_to_rownames(var = "barcode")

meta = meta %>% 
  mutate(group = str_split(sample, pattern = "-")) %>% 
  mutate(condition = map(group, pluck, 1) %>% as.character()) %>% 
  select(-group)

seu_f@meta.data = meta

Idents(seu_f) = seu_f$sample

seu_f <- seu_f %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c('mt_ratio')) %>% 
  RunPCA()

seu_f <- RunHarmony(seu_f, group.by.vars = "sample", plot_convergence = T)

seu_f <- RunUMAP(seu_f, reduction = "harmony", dims = 1:30)

seu_f <- FindNeighbors(seu_f, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(res = c(0.8))

DimPlot(seu_f, group.by = "condition")

DimPlot(seu_f, label = T)

markers = FindAllMarkers(seu_f, logfc.threshold = 0.75, min.pct = 0.25, only.pos = T)

markers_top20 = markers %>% 
  filter(avg_log2FC > 0) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)

seu_mye = subset(seu_f, RNA_snn_res.0.8 %notin% c(16, 22, 4, 20))

seu_mye <- seu_mye %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c('mt_ratio')) %>% 
  RunPCA()

seu_mye <- RunHarmony(seu_mye, group.by.vars = "sample", plot_convergence = T)

seu_mye <- RunUMAP(seu_mye, reduction = "harmony", dims = 1:30)

seu_mye <- FindNeighbors(seu_mye, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(res = c(0.7, 0.8))

markers = FindAllMarkers(seu_mye, logfc.threshold = 0.75, min.pct = 0.25, only.pos = T)

markers_top20 = markers %>% 
  filter(avg_log2FC > 0) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)

DimPlot(seu_mye, group.by = "RNA_snn_res.0.8")
DimPlot(seu_mye, group.by = "condition", split.by = "condition") + 
  scale_color_npg()
FeaturePlot(seu_mye, features = c("Spp1", "C1qa"))
FeaturePlot(seu_mye, features = "Spp1", split.by = "condition")

VlnPlot(seu_mye, features = "Spp1", group.by = "condition", pt.size = 0)
seu_mye = subset(seu_mye, RNA_snn_res.0.8 %notin% c(5, 8, 14, 15, 17))
saveRDS(seu_mye, file = "seu_mye.rds")

VlnPlot(seu_mye, features = "Itgam")

VlnPlot(seu_mye, features = "nFeature_RNA")

markers = FindAllMarkers(seu_mye, only.pos = TRUE, logfc.threshold = .5, min.pct = .1)

DimPlot(seu_mye, label = T)

seu_mye = subset(seu_mye, RNA_snn_res.0.8 %notin% c(10, 13, 16, 11, 12))

markers = markers %>% 
  filter(pct.1 > 0.1)
markers_top = markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 20)
library(RColorBrewer)
FeaturePlot(seu_mye, features = c("Spp1"), pt.size = 0.5)

avg_spp1 = AverageExpression(seu_mye, features = "Spp1", group.by = "condition")

avg_spp1 = avg_spp1$RNA %>% as.matrix()

avg_spp1 = t(avg_spp1) %>% as.data.frame()

avg_spp1 = avg_spp1 %>% 
  rownames_to_column(var = "group")

avg_spp1$group = factor(avg_spp1$group, levels = c("TREM2WT", "TREM2KO"))

ggplot(avg_spp1, aes(x = group, y = V1)) + 
  geom_col(aes(fill = group)) + 
  ggpubr::theme_classic2() + 
  coord_flip()
