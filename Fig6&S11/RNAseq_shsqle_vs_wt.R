library(data.table)
library(DESeq2)
library(tidyverse)
library(paletteer)
library(ggsci)
library(pheatmap)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichplot)
library(data.table)

counts = fread("All.HTSeq.counts.txt")
counts = as.data.frame(counts) %>% 
  column_to_rownames(var = "AccID")

counts = counts[, c(-1, -5)]

colData = data.frame(row.names = colnames(counts), group = c(rep("WT", 3), rep("SH", 3)))

colData$group = factor(colData$group, levels = c("WT", "SH"))

dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)

rld = rlog(dds, blind = T)
pca = rld@assays
plotPCA(rld, intgroup='group')

# DE analysis

dds = DESeq(dds)

res_SH = results(dds, contrast = c('group', 'SH', 'WT'), alpha = 0.05)
resultsNames(dds)
res_SH = lfcShrink(dds, coef = "group_SH_vs_WT", res = res_SH, type = 'apeglm')
res_SH_df = as.data.frame(res_SH) %>% 
  filter(!is.na(padj))

res_SH_df = res_SH_df %>% 
  arrange(log2FoldChange)

res_SH_df = res_SH_df %>% 
  filter(padj < 0.05)

res_sh_up = res_SH_df %>% 
  filter(log2FoldChange > 0)

res_sh_down = res_SH_df %>% 
  filter(log2FoldChange < 0)

rownames(res_sh_down)[rownames(res_sh_down) %like% "Ccl"]
rownames(res_sh_up)[rownames(res_sh_up) %like% "Ccl"]

res_sh_up["Apob", ]

# see specific genes expression 

gene = "Hnf4a" # exapmle
d <- plotCounts(dds, gene=gene, intgroup="group", 
                returnData=TRUE)

ggplot(d, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

# specific gene heatmap

library("pheatmap")

gene_int = c("Hmgcr", "Fdps", "Hmgcs1", "Fdft1", "Fasn", "Acaca", "Acsl3", "Acss2", 
             "Fabp5", "Fabp4", "Cd36", "Slc27a4", "Vldlr")

pheatmap(assay(dds)[gene_int,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=FALSE, scale = "row", border_color = NA)
