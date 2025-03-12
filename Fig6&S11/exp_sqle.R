library(Seurat)
library(tidyverse)

seu = subset(seu, sampletype == "snRNA-seq" & clusters != "Cycling cells")
scRNA = subset(seu, sampletype == "scRNA-seq" & clusters != "Cycling cells")

VlnPlot(scRNA, features = "SQLE", group.by = "clusters", pt.size = 0)
VlnPlot(seu, features = "SQLE", group.by = "clusters")

seu = NormalizeData(seu)
