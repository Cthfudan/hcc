# DEG and enrichment analysis between different spatial region

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(data.table)
library(RColorBrewer)

seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/Seurat/PT2_tumor_region_refined.rds")
seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/Seurat/TI1_seurat_integrated_high_res_tumor_refined.rds")
seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/Seurat/Public_seurat_integrated_high_res_tumor.rds")

score = seu@assays$predictions_sn@data

score = t(score) %>% 
  as.data.frame() %>% 
  mutate(group = case_when(
    `inflammatory` > quantile(`inflammatory`, .5) ~ "High", 
    TRUE ~ "Low"
  ))

seu$group = score$group
SpatialDimPlot(seu, group.by = "group")

DefaultAssay(seu) = "Spatial"
seu = NormalizeData(seu)
Idents(seu) = seu$group
markers = FindMarkers(seu, test.use = "MAST", ident.1 = "High", logfc.threshold = .1)
markers = markers[!rownames(markers) %like% "^MT", ]
markers = markers[!rownames(markers) %like% "^RP[LS]", ]

# Volcano plot (for DEG)
source("~/R/plotting_functions.R")

plot_volcano(markers, p_threshold = .0001, logFC_threshold = 0.25)

gene_int = c("LYZ", "CSTA", "CD44", "TGFB1", "FN1", "COL1A1", "ANXA1","THY1","IGFBP4", "SOD2", "COL1A1", "APOC1", "FDFT1", "COX6C", "IGKC", "HLA-B","GPX4", "DBI", "NDUFAB1", "FDPS") # PT2

gene_int = c("CFHR3", "C5", "C9", "A2M", "SAA4", "SERPINA3", "HP", "COX6A1", "APOM", "APOC3", "LSR", "XBP1", "CYP7B1", "DBI", "SELENOS", "LPIN2", "ANG", "PCK1", "LRP5", "LGALS4", "IL32", "TF", "CES1")

plot_volcano(markers, p_threshold = 0.0001, logFC_threshold = 0.25,label = T, genes_to_label = gene_int, highlight = T) + theme_classic()

markers_up = markers %>% filter(avg_log2FC > 0)
# Enrichment analysis

library(org.Hs.eg.db)
library(clusterProfiler)

markers_sig = markers %>% 
  filter(p_val_adj < 0.01)
gene_low = rownames(markers[markers$avg_log2FC < 0, ])
gene_high = rownames(markers[markers$avg_log2FC > 0, ])


enrich_high = enrichGO(gene = gene_high, keyType = "SYMBOL", ont = "all", OrgDb = org.Hs.eg.db)

enrichplot::dotplot(enrich_high, showCategory = 30)

enrich_high_res = as.data.frame(enrich_high)

enrich_high_res %>% filter(Description %like% "compleme")
enrich_high_res %>% filter(Description %like% "comp")
enrich_high_res %>% filter(Description %like% "chole")
enrich_high_res %>% filter(Description %like% "metabo")

pathway_int = c("triglyceride-rich plasma lipoprotein particle", "lipid homeostasis", "acylglycerol biosynthetic process", "cytochrome-c oxidase activity", "protein-lipid complex")
pathway_int = c("complement activation", "acute inflammatory response", "acute-phase response", "complement activation, classical pathway", "regulation of acute inflammatory response")
my_pal = colorRampPalette(brewer.pal(n = 9, name = "Greens"))
enrich_plot = enrich_high_res %>% 
  filter(Description %in% pathway_int)
enrich_plot = enrich_plot %>% 
  arrange(desc(qvalue)) %>% 
  mutate(Description = factor(Description, levels = unique(Description)))

ggplot(enrich_plot) + 
  geom_col(aes(x = Description, y = -log10(qvalue), fill = -log10(qvalue))) + 
  coord_flip() + 
  theme_classic() + 
  scale_fill_gradientn(colors = my_pal(100)[30:90], limits = c(2, 9))
  

enrich_low = enrichGO(gene = gene_low, keyType = "SYMBOL", ont = "all", OrgDb = org.Hs.eg.db)

enrichplot::dotplot(enrich_low)

enrich_low_res = as.data.frame(enrich_low)

enrich_low_res %>% filter(Description %like% "lip")
enrich_low_res %>% filter(Description %like% "meta")
enrich_low_res %>% filter(Description %like% "oxida")

pathway_int = c("epithelial to mesenchymal transition", "stem cell differentiation", "stem cell population maintenance", "cell-cell signaling by wnt", "epithelial cell-cell adhesion")
my_pal = colorRampPalette(brewer.pal(n = 9, name = "Reds"))
enrich_plot = enrich_low_res %>% 
  filter(Description %in% pathway_int)
enrich_plot = enrich_plot %>% 
  arrange(desc(qvalue)) %>% 
  mutate(Description = factor(Description, levels = unique(Description)))

ggplot(enrich_plot) + 
  geom_col(aes(x = Description, y = -log10(qvalue), fill = -log10(qvalue))) + 
  coord_flip() + 
  theme_classic() + 
  scale_fill_gradientn(colors = my_pal(100), limits = c(0.5, 4))

