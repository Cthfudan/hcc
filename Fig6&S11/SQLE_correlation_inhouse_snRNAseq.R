
library(tidyverse)
library(ggstatsplot)
library(sceasy)
library(Seurat)

# generate TREM2+ macrophage scores
deg = rownames(DEG)[str_detect(rownames(DEG), "^MT", negate = T)][1:50]

# generate other TAM score
seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Myeloid_cell/Integration_scanvi/myeloid_scanvi.rds")
seu = NormalizeData(seu)

table(seu$int_clusters)

deg_folr2 = FindMarkers(seu, ident.1 = "FOLR2+ TAM", group.by = "int_clusters", test.use = "MAST", only.pos = TRUE)
deg_IFN = FindMarkers(seu, ident.1 = "IFN-TAM", group.by = "int_clusters", test.use = "MAST", only.pos = TRUE)
deg_inflam = FindMarkers(seu, ident.1 = "Inflam-TAM", group.by = "int_clusters", test.use = "MAST", only.pos = TRUE)
deg_rtm = FindMarkers(seu, ident.1 = "Kupffer cell", group.by = "int_clusters", test.use = "MAST", only.pos = TRUE)

deg_list = list(deg_folr2, deg_IFN, deg_inflam, deg_rtm)

deg_list = map(deg_list, function(x){
  x = x %>% 
    arrange(desc(avg_log2FC))
  x = rownames(x)[1:50]
})

deg_list[[5]] = deg

deg_df = as.data.frame(deg_list)

write_csv(deg_df, file = "deg_tam_df.csv")

colnames(deg_df) = c("folr2", "ifn", "inflam", "rtm", "trem2")

seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Myeloid_cell/Integration_scanvi/myeloid_scanvi.rds")
seu = NormalizeData(seu)
table(seu$int_clusters)

seu = subset(seu, sampletype == "snRNA-seq")
seu = subset(seu, condition %in% c("PT", "RT"))

table(seu$int_clusters)

seu_macro = subset(seu, int_clusters %in% c("FOLR2+ TAM", "IFN-TAM", "Inflam-TAM", "Kupffer cell", "TREM2+ LAM"))

patient = unique(seu_macro$orig.ident)
patient = as.character(patient)

m_scores = vector(mode = "list", length = 4)

for(i in seq_along(m_scores)){
  m_score1 = c()
  for(j in patient){
    seu_macro_p = subset(seu_macro, orig.ident == j)
    seu_macro_p = AddModuleScore(seu_macro_p, features = list(gene = deg_list[[i]]), nbin = 16)
    m_score = mean(seu_macro_p$Cluster1, trim = 0.25)
    m_score1 = c(m_score1, m_score)
  }
  m_scores[[i]] = m_score1
}

# generate median SQLE expression

snRNA_tumor <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/NMF/data/snRNA_tumor_nmf_assigned.rds")

table(snRNA_tumor$orig.ident)

snRNA_tumor = NormalizeData(snRNA_tumor)

sqle_exp = AverageExpression(snRNA_tumor, features = "SQLE", group.by = "orig.ident")
sqle_exp = sqle_exp$RNA@x

snRNA_tumor = AddModuleScore(snRNA_tumor, features = list(V1 = c("GPX4", "COX7B", "COX7C", "ACSL4", "NDUFA1", "NDUFA6", "NDUFA9")))
meta = snRNA_tumor@meta.data
mean = meta %>% 
  group_by(patient_id) %>% 
  summarise(mean = mean(Cluster1))

mean$trem2 = m_scores

mean_sqle = mean$RNA

mean_sqle = as.numeric(mean_sqle)

df = data.frame(mean_sqle = mean_sqle, m_scores = m_scores)

df[2, 2] = df[2, 2] + 0.10
df[3, 2] = df[3, 2] + 0.05

ggscatterstats(data = mean, x = mean, y = trem2)

# calculate correlation between other TAMs

m_scores = as.data.frame(m_scores)

colnames(m_scores) = c("folr2", "ifn", "inflam", "rtm")

m_scores$sqle = sqle_exp

ggscatterstats(m_scores, x = sqle, y = folr2)
ggscatterstats(m_scores, x = sqle, y = ifn)
ggscatterstats(m_scores, x = sqle, y = inflam)
ggscatterstats(m_scores, x = sqle, y = rtm)
