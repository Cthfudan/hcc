library(tidyverse)
library(Seurat)
library(data.table)
library(harmony)
library(ggpubr)
library(ggsci)

seu_NKT = subset(seu_anno, celltype == "NK/T cells")

# standard pipline

standard_pipline_harmony = function(seu, dims = 1:20, resolution = 0.8){
  seu = seu %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData(vars.to.regress = "mt_ratio") %>% 
    RunPCA()
  
  seu = RunHarmony(seu, group.by.vars = "orig.ident", plot_convergence = T, max.iter.harmony = 20)
  
  seu = RunUMAP(seu, reduction = "harmony", dims = dims) %>% 
    FindNeighbors(reduction = "harmony", dims = dims) %>% 
    FindClusters(resolution = resolution)
  
  return(seu)
}

seu = standard_pipline_harmony(seu_NKT, dims = 1:20, resolution = .8)

# findmarkers

markers_T = FindAllMarkers(seu, only.pos = T, logfc.threshold = .5, test.use = "wilcox")

markers_T_top20 = markers_T %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 20)

saveRDS(seu, file = "data/seu_PanT.rds")
saveRDS(markers_T, file = "data/markers_PanT.rds")

# subset to CD8 T
seu_cd8 = subset(seu, RNA_snn_res.0.8 %in% c("0", "2", "3", "4", "5", "10", "11", "12"))
seu_cd8 = standard_pipline_harmony(seu = seu_cd8, dims = 1:20, resolution = .6)
markers_CD8T = FindAllMarkers(seu_cd8, only.pos = T)
markers_CD8T_top20 = markers_CD8T %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 20)

saveRDS(seu_cd8, file = "data/seu_CD8T.rds")
saveRDS(markers_CD8T, file = "data/markers_CD8T.rds")

# CD8 T annotation

seu_cd8_anno = RenameIdents(seu_cd8, 
                            "0" = "Early effector memory T cells (CCR6+ CD8T)", 
                            "1" = "Effector memory T cells (GZMK+ CD8T)", 
                            "2" = "MAIT (CD161+ CD8T)", 
                            "3" = "Effector T (IFNG+ CD8T)", 
                            "4" = "NK cells", 
                            "5" = "NKT cells", 
                            "6" = "Effector T (TNF+ CD8T)", 
                            "7" = "Effector T (IFNG+ CD8T)", 
                            "8" = "CXCL13+ CD8T", 
                            "9" = "Terminal effector T (PRF1+ GZMB+ CD8T)", 
                            "10" = "NKT cells", 
                            "11" = "Exhausted T (TIGIT+ CD8T)", 
                            "12" = "MAIT (CD161+ CD8T)")

DimPlot(seu_cd8_anno) + scale_color_npg()
DimPlot(seu_cd8_anno, group.by = "response") + scale_color_npg()
seu_cd8_anno$celltype = seu_cd8_anno@active.ident

# tissue preference using OR value ...