library(tidyverse)
library(Seurat)
library(liana)

seu <- readRDS("~/Projects/SYF_pd-1/project/landscape/seu_full_anno.rds")

table(seu$celltype)

seu = subset(seu, condition != "P") # remove normal condition
seu_b = subset(seu, condition == "B")
seu_r = subset(seu, condition == "R")
seu_nr = subset(seu, condition == "NR")
show_resources()
show_methods()

seu = SetIdent(seu, value = seu$celltype)

liana_b = liana_wrap(seu_b)
liana_r = liana_wrap(seu_r)
liana_nr = liana_wrap(seu_nr)

glimpse(liana_test)

liana_b = liana_aggregate(liana_b)
liana_r = liana_aggregate(liana_r)
liana_nr = liana_aggregate(liana_nr)
liana_b_2 = liana_b %>% 
  filter(cellphonedb.pvalue < 0.05)
liana_nr_2 = liana_nr %>% 
  filter(cellphonedb.pvalue < 0.05)
g1 = heat_freq(liana_b_2)
g2 = heat_freq(liana_nr_2)

g1_mat = g1@matrix
g2_mat = g2@matrix
g3 = g2_mat/g1_mat
pheatmap::pheatmap(g3)
liana_test %>%
  liana_dotplot(source_groups = c("TREM2+ SPP1+ macrophage"),
                target_groups = c("Hepatocytes", "CD8+ T"),
                ntop = 20)

liana_spp1_b = liana_b %>% 
  filter(ligand.complex == "SPP1")
liana_spp1_r = liana_r %>% 
  filter(ligand.complex == "SPP1")
liana_spp1_nr = liana_nr %>% 
  filter(ligand.complex == "SPP1")

liana_dotplot(liana_spp1_b, 
              target_groups = c("B cell", "CAF", "CD4+ T", "CD8+ T", "EC", "Hepatocytes", "NK cell", "NKT cell", "Plasma cell"), magnitude = "logfc.logfc_comb")

liana_dotplot(liana_spp1_nr, source_groups = c("TREM2+ SPP1+ macrophage"), 
              target_groups = c("B cell", "CAF", "CD4+ T", "CD8+ T", "EC",  "Hepatocytes", "NK cell", "NKT cell", "Plasma cell"), magnitude = "logfc.logfc_comb")

eliana_trunc <- liana_spp1 %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected


seu_hep = subset(seu, celltype == "Hepatocytes")

seu_cd8t = subset(seu, celltype == "CD8+ T")

DotPlot(seu_hep, features = c("PTGER4", "CD44", "S1PR1", "ITGB1", "ITGAV", "ITGA5", "ITGA4", "ITGA8", "ITGA9"), group.by = "condition", scale = TRUE) + 
  scale_color_viridis_c(option = "D")
DotPlot(seu_cd8t, features = c("PTGER4", "CD44", "S1PR1", "ITGB1", "ITGAV", "ITGA5", "ITGA4", "ITGA8", "ITGA9"), group.by = "condition", scale = TRUE)+ 
  scale_color_viridis_c(option = "D")
seu_mye = subset(seu, celltype %in% c("FOLR2+ macrophage", "RTM-like TAM", "TREM2+ SPP1+ macrophage"))
seu_T = subset(seu, celltype == "CD8+ T")
VlnPlot(seu_cd8t, features = c("TNF"), group.by = "condition", pt.size = 0) + ggsci::scale_fill_npg()
VlnPlot(seu_cd8t, features = c("IFNG"), group.by = "condition", pt.size = 0) + ggsci::scale_fill_npg()
DotPlot(seu_mye, features = c("SPP1"), group.by = "orig.ident", scale = TRUE)
liana_heatmap()

VlnPlot(seu_CD8T, features = c("IFNG", "TNF"), group.by = "treatment_response", pt.size = 0)

seu_cd8t_anno = RenameIdents(seu_CD8T, 
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

DotPlot(seu_cd8t_anno, features = c("IFNG", "TNF"))

seu_cd8t_anno$celltype = seu_cd8t_anno@active.ident
meta = seu_cd8t_anno@meta.data

meta = meta %>% 
  filter(treatment_response %in% c("after_non-responsive", "after_responsive"))
