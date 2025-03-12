res_ox_CM_table = res_ox_CM_table %>% 
  arrange(desc(log2FoldChange))

res_sig = res_ox_CM_table %>% 
  filter(log2FoldChange > 1)

res_sig = rownames(res_sig)

# GSEA analysis on scRNA-seq datasets

library(sceasy)
library(Seurat)
seu <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/mouse_myeloid/mye_anno.rds")
Idents(seu) = seu$annotation2
markers = FindMarkers(seu, ident.1 = "Trem2+ Spp1+ LAM", test.use = "MAST", logfc.threshold = .01)
saveRDS(markers, file = "markers_Trem2+tam.rds")
markers = markers %>% 
  arrange(desc(avg_log2FC))

gene_list = markers$avg_log2FC

names(gene_list) = rownames(markers)

## GSEA analysis

geneset = data.frame(term = "oxLDL-up signature", gene = res_sig)

gene_list = sort(gene_list, decreasing = TRUE)

gene_list = gene_list[is.finite(gene_list)]

gsea_result = GSEA(geneList = gene_list, TERM2GENE = geneset, pvalueCutoff = 1, pAdjustMethod = 'BH', nPermSimple = 10000, eps = 0)

gsea_table = as.data.frame(gsea_result)

gseaplot2(gsea_result, geneSetID = rownames(gsea_table)[1], 
          title = '', 
          color = '#44995F', 
          pvalue_table = F, 
          ES_geom = 'line', 
          base_size = 14, 
          rel_heights = c(1, .2, .4)
)

gsea_table_final = gsea_table # used to combine multiple gsea result

## GSEA analysis on multiple celltype

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
  gsea_result = GSEA(geneList = gene_list, TERM2GENE = geneset, pvalueCutoff = 1, pAdjustMethod = 'BH', nPermSimple = 10000, eps = 0)
  gsea_table = as.data.frame(gsea_result)
  gsea_table_final = rbind(gsea_table_final, gsea_table)
}
rownames(gsea_table_final) = c("Trem2+ LAM", "RTM-like TAM", "Folr2+ TAM", "Inflam-TAM", "Cycling TAM")
gsea_table_final = gsea_table_final %>% 
  rownames_to_column(var = "celltype")

gsea_table_final = gsea_table_final %>% 
  mutate(logfdr = -log10(p.adjust))

gsea_table_final = gsea_table_final %>% 
  arrange(desc(NES))

gsea_table_final$celltype = factor(gsea_table_final$celltype, levels = rev(unique(gsea_table_final$celltype)))
  
ggplot(gsea_table_final, aes(x = NES, y = celltype, size = logfdr)) + 
  geom_point() + 
  coord_cartesian(xlim = c(-3, 3)) +
  ggpubr::theme_classic2() + 
  geom_vline(xintercept = 0, linetype = 2)

# heatmap 

library(pheatmap)

gene_int = enrich_gene
exp_mat = counts(dds, normalize = T)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(n = 6, name = "RdBu"))(100)

gene_to_label = c("Spp1", "Edn1", "Fabp3", "Trem2", "Ptges", "Furin", "Hrh2", "Ptgs2", "Il10", "Src", "Bnip3", "Tnfaip3", "Rcan1", "Cxcl9", "Ccl22", "Ccl5", "Il1rn")

row_label = rep("", length = length(gene_int))

names(row_label) = gene_int

row_label[gene_to_label] = gene_to_label

p = pheatmap(exp_mat[gene_int, ], cluster_rows=T, 
             cluster_cols=F, scale = "row", border_color = NA, color = rev(mycol), labels_row = row_label)

