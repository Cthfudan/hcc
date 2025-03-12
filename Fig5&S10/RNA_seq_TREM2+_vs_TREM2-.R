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

counts = fread("All_reads_counts.csv")
counts = as.data.frame(counts) %>% 
  column_to_rownames(var = "Geneid")

# select subset of counts if needed
counts = counts[, str_detect(colnames(counts), "ox")]

# create dds obj

colData = data.frame(row.names = colnames(counts), group = c(rep("Tn", 2), rep("WT", 2)))

colData$group = factor(colData$group, levels = c("WT", "Tn"))

dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)

# QC analysis

rld = rlog(dds, blind = T)
pca = rld@assays
p = plotPCA(rld, intgroup='group')
pca_data = p$data
ggplot(data = pca_data) + 
  geom_point(aes(x = PC1, y = PC2, color = group), size = 6, shape = 19) + 
  scale_color_npg(alpha = 0.5) + 
  theme_bw() + 
  coord_cartesian(ylim = c(-9, 9), xlim = c(-25, 35))


rld_cov = assay(rld) %>% 
  cov()

pheatmap(rld_cov)

# DE analysis
dds = DESeq(dds)

res = results(dds, contrast = c('group', 'Tn', 'WT'), alpha = 0.1)
resultsNames(dds)
res = lfcShrink(dds, coef = "group_Tn_vs_WT", res = res, type = 'apeglm')
res_table = as.data.frame(res) %>% 
  filter(!is.na(padj))

res_table = res_table %>% 
  arrange(desc(log2FoldChange))

res_table = res_table %>% 
  filter(!rownames(res_table) %like% "^Gm")

res_table = res_table %>% 
  filter(!rownames(res_table) %like% "^ENSMUSG")

res_table = res_table %>% 
  filter(!rownames(res_table) %like% "^[0-9]")

res_table = res_table %>% 
  filter(!rownames(res_table) %like% "Rik$")
# heatmap vis

sig_up = res_table %>% 
  filter(log2FoldChange > 0.5 & padj < 0.1)

sig_down = res_table %>% 
  filter(log2FoldChange < 0 & padj < 0.05) %>% 
  arrange(log2FoldChange)

sig_up = rbind(sig_up, res_table[gene_int_trem2_down, ])
sig_down = rbind(sig_down, res_table[gene_int_rtm, ])
library(pheatmap)

exp_mat = counts(dds, normalize = TRUE)

gene_int_trem2_down = c("Ldlr", "Lpl", "Mexis", "Lyst", "Ncor1", "Acaca", "Abca1", "Abca3", "Fabp3", "Lrp1", "Acsl4", "Cd33", "Ciita", "Cd74", "Trem2", "Pilrb1")

genes_to_label = c(gene_int_rtm, gene_int_trem2_down)

gene_int = unique(c(rownames(sig_down)[1:25], gene_int_rtm, rownames(sig_up)[1:25], gene_int_trem2_down))

label_row = rep("", length(gene_int))

names(label_row) = gene_int

label_row[genes_to_label] = genes_to_label
pheatmap(exp_mat[gene_int, ], scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, labels_row = label_row, border_color = NA)

# GSEA analysis

sig_up = res_table %>% 
  filter(log2FoldChange > 0.5 & padj < 0.1)

sig_down = res_table %>% 
  filter(log2FoldChange < -0.5 & padj < 0.1)

markers <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/RNAseq_BMDM/markers_Trem2+tam.rds")

markers = markers %>% 
  arrange(desc(avg_log2FC))

gene_list = markers$avg_log2FC

names(gene_list) = rownames(markers)

## GSEA analysis
gene_up = rownames(sig_up)
gene_down = rownames(sig_down)
geneset = data.frame(term = c("Trem-_up"), gene = c(gene_up))
geneset = data.frame(term = c("Trem-_down"), gene = c(gene_down))

gene_list = sort(gene_list, decreasing = TRUE)

gene_list = gene_list[is.finite(gene_list)]

gsea_result = GSEA(geneList = gene_list, TERM2GENE = geneset, pvalueCutoff = 1, pAdjustMethod = 'BH', nPermSimple = 20000, eps = 0)

gsea_table = as.data.frame(gsea_result)

gseaplot2(gsea_result, geneSetID = rownames(gsea_table)[1], 
          title = '', 
          color = '#44995F', 
          pvalue_table = F, 
          ES_geom = 'line', 
          base_size = 14, 
          rel_heights = c(1, .2, .4) # 小图相对高度
)

gsea_table_final = gsea_table # used to combine multiple gsea result

## GSEA analysis on all TAM subsets
library(Seurat)
markers_rtm = FindMarkers(seu, ident.1 = "RTM-like TAM", test.use = "MAST", logfc.threshold = .01)
markers_folr2 = FindMarkers(seu, ident.1 = "Folr2+ TAM", test.use = "MAST", logfc.threshold = .01)
markers_inflam = FindMarkers(seu, ident.1 = "Inflam-TAM", test.use = "MAST", logfc.threshold = .01)
markers_cycling = FindMarkers(seu, ident.1 = "Cycling cells", test.use = "MAST", logfc.threshold = .01)

markers_list = list(markers_rtm, markers_folr2, markers_inflam, markers_cycling)

for(i in markers_list){
  i = i %>% 
    arrange(desc(avg_log2FC))
  gene_list = i$avg_log2FC
  names(gene_list) = rownames(i)
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[is.finite(gene_list)]
  gsea_result = GSEA(geneList = gene_list, TERM2GENE = geneset, pvalueCutoff = 1, pAdjustMethod = 'BH', nPermSimple = 20000, eps = 0)
  gsea_table = as.data.frame(gsea_result)
  gsea_table_final = rbind(gsea_table_final, gsea_table)
}
rownames(gsea_table_final) = c("Trem2+ LAM", "RTM-like TAM", "Folr2+ TAM", "Inflam-TAM", "Cycling TAM")
gsea_table_final = gsea_table_final %>% 
  rownames_to_column(var = "celltype")

gsea_table_final = gsea_table_final %>% 
  mutate(logfdr = -log10(pvalue))

gsea_table_final = gsea_table_final %>% 
  arrange(desc(NES))

gsea_table_final$celltype = factor(gsea_table_final$celltype, levels = rev(unique(gsea_table_final$celltype)))

gsea_table_final$celltype = factor(gsea_table_final$celltype, levels = rev(unique(gsea_table_final$celltype)))
ggplot(gsea_table_final, aes(x = NES, y = celltype, size = logfdr)) + 
  geom_point() + 
  coord_cartesian(xlim = c(-3, 3)) +
  ggpubr::theme_classic2() + 
  geom_vline(xintercept = 0, linetype = 2)
