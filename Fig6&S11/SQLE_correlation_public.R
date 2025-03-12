library(data.table)
library(Seurat)
library(tidyverse)

seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/svm/data/GSE149614/seu.rds")

DimPlot(seu, label = T)

FeaturePlot(seu, features = "APOM")

seu_hep = subset(seu, seurat_clusters %in% c(3, 9, 12, 22, 25, 24, 33, 26, 28, 17, 34))

DimPlot(seu_hep, label = T)

seu_hep = subset(seu_hep, site == "Tumor")

seu_hep = NormalizeData(seu_hep)

sqle_exp = AverageExpression(seu_hep, features = "SQLE", group.by = "patient")

seu_hep = AddModuleScore(seu_hep, features = list(V1 = c("GPX4", "COX7B", "COX7C", "ACSL4", "NDUFA1", "NDUFA6", "NDUFA9")))
sqle_exp = sqle_exp$RNA
sqle_exp = as.numeric(sqle_exp)

meta = seu_hep@meta.data
mean = meta %>% 
  group_by(orig.ident) %>% 
  summarise(mean = mean(Cluster1))

seu_mye = convertFormat("~/python/snRNA_scRNA_hcc/SVM/GSE_mye_filtered_raw.h5ad", from = "anndata", to = "seurat", main_layer = "counts")

table(seu_mye$celltype)

seu_mye = subset(seu_mye, site == "Tumor")

seu_macro = subset(seu_mye, celltype %in% c("FOLR2+ TAM", "IFN-TAM", "Inflam-TAM", "TREM2+ TAM"))

patient = as.character(unique(seu_mye$patient))

m_scores = vector(mode = "list", length = 4)

for(i in seq_along(m_scores)){
  m_score1 = c()
  for(j in patient){
    seu_macro_p = subset(seu_macro, orig.ident == j)
    seu_macro_p = AddModuleScore(seu_macro_p, features = list(gene = deg_list[[i]]))
    m_score = mean(seu_macro_p$Cluster1, trim = 0.25)
    m_score1 = c(m_score1, m_score)
  }
  m_scores[[i]] = m_score1
}

df = data.frame(mean_sqle = sqle_exp, m_scores = m_scores)

mean$trem2 = m_scores

ggscatterstats(data = df, x = mean_sqle, y = m_scores, type = "robust")
ggscatterstats(data = mean, x = mean, y = trem2)

# sqle correlation with other TAMs
m_scores = as.data.frame(m_scores)

colnames(m_scores) = c("folr2", "ifn", "inflam", "rtm")

m_scores$sqle = sqle_exp

ggscatterstats(m_scores, x = sqle, y = folr2)
ggscatterstats(m_scores, x = sqle, y = ifn)
ggscatterstats(m_scores, x = sqle, y = inflam)
ggscatterstats(m_scores, x = sqle, y = rtm)
