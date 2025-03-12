library(Seurat)
library(data.table)

genes = fread("scFEA.human.genes.txt", header = FALSE)
seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/Seurat/PT2_seurat_integrated_high_res.rds") #spatial
SpatialDimPlot(seu, label = T)

DefaultAssay(seu) = "Spatial"

seu = NormalizeData(seu)

exp = seu@assays$Spatial$data

exp = as.matrix(exp)

exp = exp[rownames(exp) %in% genes$V1, ]

write.csv(exp, file = "exp_PT2_all_spatial.csv", row.names = TRUE, col.names = TRUE)

