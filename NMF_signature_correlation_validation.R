# program expression correlation - validation on public and mouse dataset

library(tidyverse)
library(corrplot)
library(Seurat)
library(pagoda2)
library(RColorBrewer)
library(BBmisc)
library(harmony)
"%notin%" = Negate("%in%")

source("~/R/source_functions.R")

# 1. GSE149614 dataset ----------------------------------------------------

# read data and preprocessing

GSE149614 = standard_QC_human(GSE149614, mt_thres = 20, count_thres = 400, feature_thres = 200)
GSE149614 = subset(GSE149614, site %in% c("Normal", "Tumor"))
GSE149614 = standard_pipline_harmony(GSE149614, group.by = "sample", dims = 1:20, resolution = 1.0, regress.mt = F)

DimPlot(GSE149614)

GSE149614_hep = subset(GSE149614, RNA_snn_res.1 %in% c("1", "5", "12", "21", "16", "17", "29", "18") & site == "Tumor")

# score the metaprogram on GSE149614

GSE149614_hep = SCTransform(GSE149614_hep)
exp_mat <- GetAssayData(GSE149614_hep, slot = "data")
mito.genes <- grep(pattern = "^MT-", x = rownames(exp_mat), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(exp_mat), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(exp_mat), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(exp_mat), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(exp_mat), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(exp_mat), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(exp_mat), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(exp_mat), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(exp_mat), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(exp_mat), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(exp_mat), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(exp_mat), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(exp_mat), value = TRUE)

exp_mat <- exp_mat[-(which(rownames(exp_mat) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]

score_meta = map(metaGene, ~ score.cells.puram(data = t(exp_mat), signature = .x))
score_meta = as.data.frame(score_meta)
colnames(score_meta) = str_c("sig", 1:8)

# corrplot

score_cor = cor(score_meta, method = "pearson")
score_cor = score_cor[c(2,1,3,5,4,6,7,8), c(2,1,3,5,4,6,7,8)]
corrplot(corr = score_cor, 
         method = "color", 
         type = "lower", 
         order = "original", 
         addCoef.col = "black",
         tl.srt = 45,
         tl.col = "black", 
         is.corr = T,
         hclust.method = "ward.D2", 
         col = colorRampPalette(rev(brewer.pal(name = 'RdBu', 11)))(50))

# 2. mouse scRNA-seq dataset ----------------------------------------------

# map human genes to mouse conterparts

mouse_metaGene = map(metaGene, ~ convert_human_to_mouse(.x))

mice_hep = subset(seu_mice, celltype == "Hepatocytes" & group == "PT")

mice_hep = SCTransform(mice_hep)

exp_mat <- GetAssayData(mice_hep, slot = "data")

# score genes

score_meta = map(mouse_metaGene, ~ score.cells.puram(data = t(exp_mat), signature = .x))
score_meta = as.data.frame(score_meta)
colnames(score_meta) = str_c("sig", 1:8)

# corrplot

score_cor = cor(score_meta, method = "pearson")
score_cor = score_cor[c(2,1,3,5,4,6,7,8), c(2,1,3,5,4,6,7,8)]
corrplot(corr = score_cor, 
         method = "color", 
         type = "lower", 
         order = "original", 
         addCoef.col = "black",
         tl.srt = 45,
         tl.col = "black", 
         is.corr = T,
         hclust.method = "ward.D2", 
         col = colorRampPalette(rev(brewer.pal(name = 'RdBu', 11)))(50))

# 3. TCGA bulk datasets ------------------------------------------

exp_mat <- readRDS("~/R/snRNA_scRNA_hcc/project/TCGA/data/normalized_exp_mat.rds")

# score using GSVA

score_meta = gsva(exp_mat, gset.idx.list = metaGene, method = "ssgsea")

score_meta = as.data.frame(score_meta) %>% 
  t()

colnames(score_meta) = str_c("sig", 1:8)

score_cor = cor(score_meta, method = "pearson")
score_cor = score_cor[c(2,1,3,5,4,6,7,8), c(2,1,3,5,4,6,7,8)]
corrplot(corr = score_cor, 
         method = "color", 
         type = "lower", 
         order = "original", 
         addCoef.col = "black",
         tl.srt = 45,
         tl.col = "black", 
         is.corr = T,
         hclust.method = "ward.D2", 
         col = colorRampPalette(rev(brewer.pal(name = 'RdBu', 11)))(50))
