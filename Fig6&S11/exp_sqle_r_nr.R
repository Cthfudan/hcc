library(Seurat)
library(tidyverse)
library(data.table)
seu <- readRDS("~/R/Projects/SYF_pd-1/project/landscape/seu_full_anno.rds")
table(seu$condition)

table(seu$celltype)

seu_hep = subset(seu, celltype == "Hepatocytes")
seu_hep = NormalizeData(seu_hep)
seu_hep$condition = factor(seu_hep$condition, levels = c("P", "B", "R", "NR"))

VlnPlot(seu_hep, features = "SQLE", group.by = "condition") + 
  ggsci::scale_fill_npg()

data = as.data.frame(data)
library(ggstatsplot)
ggstatsplot::ggbetweenstats(data, x = ident, y = SQLE)
