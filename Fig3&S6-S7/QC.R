library(tidyverse)
library(Seurat)
library(data.table)
library(harmony)
library(scrubletR)
"%notin%" = Negate("%in%")

# read data
seu = readRDS("data/seu_raw.rds")

# preprocess

seu[["mt_ratio"]] = PercentageFeatureSet(seu, pattern = "^MT-")

ggplot(data = seu@meta.data, 
       mapping = aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(mapping = aes(color = mt_ratio), size = 0.2) + 
  geom_smooth(method = "lm", se = TRUE) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(x = "log count", y = "log UMIs") + 
  geom_vline(xintercept = 10000) + 
  geom_hline(yintercept = 10000)

seu = subset(seu, nCount_RNA > 400 & nFeature_RNA > 200 & mt_ratio < 10)

table(seu$orig.ident)

# remove doublets using scrublet
seu_list = SplitObject(seu, split.by = "orig.ident")
seu_list = lapply(seu_list, function(x) scrublet_R(x, python_home = "/opt/miniconda3/envs/scverse/bin/python"))
seu_list = lapply(seu_list, function(x) subset(x, predicted_doublets == F))

# standard pipline
seu = seu %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = "mt_ratio") %>% 
  RunPCA()

seu = RunHarmony(seu, group.by.vars = "orig.ident", plot_convergence = T, max.iter.harmony = 20)

seu = RunUMAP(seu, reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 1)

# annotation

# manual annotation
markers = FindAllMarkers(seu, logfc.threshold = .75, min.pct = .25, only.pos = T, test.use = "wilcox")

markers = markers %>% 
  merge(cell_type, by = "cluster")

markers_top20 = markers %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 20)
DimPlot(seu, label = T)
FeaturePlot(seu, features = c("CD8A"))

seu_anno = RenameIdents(seu, 
                   "0" = "NK/T cells", 
                   "1" = "NK/T cells", 
                   "2" = "NK/T cells", 
                   "3" = "NK/T cells", 
                   "4" = "Hepatocytes", 
                   "5" = "Myeloid cells", 
                   "6" = "NK/T cells", 
                   "7" = "Hepatocytes", 
                   "8" = "NK/T cells", 
                   "9" = "B cells", 
                   "10" = "Myeloid cells", 
                   "11" = "Myeloid cells", 
                   "12" = "Endothelial cells", 
                   "13" = "Myeloid cells", 
                   "14" = "Endothelial cells", 
                   "15" = "Cholangiocytes", 
                   "16" = "Cycling cells", 
                   "17" = "Hepatocytes", 
                   "18" = "Hepatocytes", 
                   "19" = "Myeloid cells", 
                   "20" = "Hepatocytes", 
                   "21" = "CAF", 
                   "22" = "Myeloid cells", 
                   "23" = "Hepatocytes", 
                   "24" = "Plasma cells", 
                   "25" = "Hepatocytes", 
                   "26" = "Cholangiocytes", 
                   "27" = "Hepatocytes", 
                   "28" = "Hepatocytes", 
                   "29" = "Doublet"
                   )
seu_anno$celltype = seu_anno@active.ident
seu_anno = subset(seu_anno, celltype %notin% "Doublet")
DimPlot(seu_anno)

seu_anno$celltype = seu_anno@active.ident
seu_anno = subset(seu_anno, celltype %notin% c("Doublets", "Doublet", "low quality cells"))
seu_anno = RunUMAP(seu_anno, reduction = "harmony", dims = 1:20)
saveRDS(seu_anno, file = "data/seu_anno_filtered.rds")
saveRDS(markers, file = "data/markers_filter.rds")
