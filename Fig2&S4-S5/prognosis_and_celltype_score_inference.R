library(tidyverse)
library(BayesPrism)
library(survival)
library(survminer)
library(data.table)
"%notin%" = Negate("%in%")
# read deconvlution data

bp.result.sc = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/BayesPrism/bp_result_sc_high_res_new.rds")
bp.result.sn <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/BayesPrism/bp_result_sn.rds")

cellst.sc = get.fraction(bp.result.sc, "final", "type")
cellst.sn = get.fraction(bp.result.sn, "first", "state")

# read TCGA survival data

clinical = readxl::read_excel("~/R/Projects/recurrence_gset/data/LIHC-US/TCGA-CDR.xlsx")

clinical = clinical %>% 
  filter(type == "LIHC") %>% 
  filter(histological_type == "Hepatocellular Carcinoma") %>% 
  filter(!is.na(OS.time)) %>% 
  filter(is.na(Redaction))

sample_id = rownames(cellst.sc)

sample_id = str_split(sample_id, pattern = "-01A-")

sample_id = map(sample_id, 1) %>% unlist()

sample_id = str_split(sample_id, pattern = "-01B-")

sample_id = map(sample_id, 1) %>% unlist()

rownames(cellst.sc) = sample_id

cellst.sc = as.data.frame(cellst.sc)

rownames(cellst.sn) = sample_id

cellst.sn = as.data.frame(cellst.sn)

cellst.sc.trem2 = cellst.sc %>% 
  dplyr::select("TREM2+ LAM")

cellst.sn.archtype1 = cellst.sn %>% 
  select(archtype1)

surv_df = merge(clinical, cellst.sc.trem2, by.x = "bcr_patient_barcode", by.y = "row.names")

surv_df = merge(surv_df, cellst.sn.archtype1, by.x = "bcr_patient_barcode", by.y = "row.names")

# define cutpoint 

cut_g1 = surv_cutpoint(data = surv_df, time = 'PFI.time', event = 'PFI', variables = 'archtype1', minprop = .4)
cut_g1 = cut_g1$cutpoint$cutpoint

cut_g2 = surv_cutpoint(data = surv_df, time = 'PFI.time', event = 'PFI', variables = 'TREM2+ LAM')
cut_g2 = cut_g2$cutpoint$cutpoint

## survival on one marker
cut_g2 = surv_cutpoint(data = surv_df, time = 'PFI.time', event = 'PFI', variables = 'TREM2+ LAM')
cut_g2 = cut_g2$cutpoint$cutpoint

surv_df = mutate(surv_df, group = case_when(
  `TREM2+ LAM` >= cut_g2 ~ "High", 
  TRUE ~ "Low"
))

ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group, data = surv_df), pval = TRUE, palette = c("#E64B35", "#4DBBD5"), conf.int = F, risk.table = T)

## survival on two markers

surv_df = surv_df %>% 
  mutate(group = case_when(
    archtype1 > cut_g1 & `FOLR2+ TAM` > cut_g2 ~ "g1hg2h", 
    TRUE ~ "others"
  ))

ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group, data = surv_df), pval = TRUE, conf.int = F, risk.table = T, size = .25, censor = F, palette = c("#00A087", "#3C5488"))

# infer immune compartments using Xcell

devtools::install_local("xCell-master.zip")

library(xCell)

exp_data <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/BayesPrism/TCGA_HCC_TPM.rds")

cell_score = xCellAnalysis(exp_data, rnaseq = T)

celltype = rownames(cell_score) # select sepcific celltypes

cell_score = xCellAnalysis(exp_data, rnaseq = T, cell.types.use = c("CD4+ memory T-cells", "CD4+ naive T-cells", "CD4+ Tcm", "CD4+ Tem", 
                                                                    "CD8+ naive T-cells", "CD8+ Tcm", "CD8+ Tem", "B-cells", "Endothelial cells", 
                                                                    "Epithelial cells", "Fibroblasts", "Hepatocytes", "Macrophages M1", "Macrophages M2", 
                                                                    "Monocytes", "Neutrophils", "NK cells", "NKT", "Plasma cells", "cDC", "pDC", "Tregs"))

cell_score = xCellAnalysis(exp_data, rnaseq = T, cell.types.use = c("B-cells", "Endothelial cells", "CD4+ T-cells", "CD8+ T-cells", 
                                                                    "Epithelial cells", "Fibroblasts", "Hepatocytes", "Macrophages M1", "Macrophages M2", 
                                                                    "Monocytes", "Neutrophils", "NK cells", "NKT", "Plasma cells", "cDC", "pDC", "Tregs"))

cell_score = t(cell_score) %>% as.data.frame()

rownames(cell_score) = sample_id

cell_score = cell_score[rownames(cell_score) %in% surv_df$bcr_patient_barcode, ]

surv_df = surv_df %>% 
  mutate(group = case_when(
    archtype1 > cut_g1 & `FOLR2+ TAM` > cut_g2 ~ "g1hg2h", 
    archtype1 > cut_g1 & `FOLR2+ TAM` <= cut_g2 ~ "g1hg2l", 
    archtype1 <= cut_g1 & `FOLR2+ TAM` > cut_g2 ~ "g1lg2h", 
    TRUE ~ "g1lg2l"
  ))

cell_score$group = surv_df$group

cell_score = split(cell_score, cell_score$group)

x = map(cell_score, function(x){
  x = x[map_lgl(x, is.numeric)]
  
  m = map_dbl(x, mean, trim = .1)
  
  return(m)
})
help(mean)
x = as.data.frame(x)

x = x %>% rownames_to_column(var = "celltype")

annotation_col = data.frame(row.names = colnames(x), 
                            group = colnames(x))
pheatmap::pheatmap(x, scale = "row", color = viridis::viridis(50, option = "A"), 
                   annotation_col = annotation_col, 
                   annotation_colors = list(group = c(g1hg2h = "#E64B35FF", 
                                                      g1hg2l = "#4DBBD5FF", 
                                                      g1lg2h = "#00A087FF", 
                                                      g1lg2l = "#3C5488FF")), 
                   border_color = NA
                   )
