library(BayesPrism)
library(data.table)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(AnnoProbe)
library(Seurat)
library(DESeq2)

bp.result.sc = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/BayesPrism/bp_result_sc.rds")
bp.result.sn <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/BayesPrism/bp_result_sn_high_res.rds")

# extract cell state proportion

cellst.sc = get.fraction(bp.result.sc, which.theta = "final", state.or.type = "type")
cellst.sn = get.fraction(bp.result.sn, which.theta = "first", state.or.type = "state")

# extract celltype cv

theta.f.sc = bp.result.sc@posterior.initial.cellState@theta.cv
theta.f.sn = bp.result.sn@posterior.initial.cellState@theta.cv

cellst = cbind(cellst.sc, cellst.sn[, c("archtype1", "archtype2", "inflammatory", "Proliferation")])

# boxplot
cellst = as.data.frame(cellst)
cellst = mutate(cellst, group = case_when(
  archtype1 >= median(archtype1) ~ "High", 
  archtype1 < median(archtype1) ~ "Low"
))

ggplot(cellst %>% filter(`CD8+ Tem` < 0.00001), aes(x = group, y = `CD8+ Tem`, color = group)) + 
  geom_boxplot() + geom_jitter(aes(color = group), width = 0.2) + ggpubr::stat_compare_means(ref.group = "High") + ggsci::scale_color_npg() + ggpubr::theme_classic2()

ggplot(cellst, aes(x = group, y = `Macrophage`, color = group)) + 
  geom_boxplot() + geom_jitter(aes(color = group), width = 0.2) + ggpubr::stat_compare_means(ref.group = "High") + ggsci::scale_color_npg() + ggpubr::theme_classic2()
