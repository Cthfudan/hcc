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

# select subset of counts if needed
counts = counts[, str_detect(colnames(counts), "^(CM)|(Ox)")]

# create dds obj

count = t(counts)

write.table(count, file = "exp_mat.csv", row.names = TRUE, sep = ",")

colData = data.frame(row.names = colnames(counts), group = c(rep("CM", 3), rep("Control", 3), rep("LDL", 3), rep("oxLDL", 3)))

colData$group = factor(colData$group, levels = c("CM", "oxLDL"))

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

res_ox_CM = results(dds, contrast = c('group', 'oxLDL', 'CM'), alpha = 0.05)
resultsNames(dds)
res_ox_CM = lfcShrink(dds, coef = "group_oxLDL_vs_CM", res = res_ox_CM, type = 'apeglm')
res_ox_CM_table = as.data.frame(res_ox_CM) %>% 
  filter(!is.na(padj))

res_ox_CM_table = res_ox_CM_table %>% 
  arrange(desc(log2FoldChange))

# enrichment analysis
deg_ox_CM = res_ox_CM_table %>% 
  dplyr::filter(padj <= 0.05) %>% 
  dplyr::filter(log2FoldChange > 0)

gene_id = mapIds(keys = c(rownames(deg_ox_CM)), x = org.Mm.eg.db, keytype = 'SYMBOL', column = 'ENTREZID')

# enrichgo
enrich_ox_CM_go = enrichGO(gene = c(rownames(deg_ox_CM), "Trem2", "Spp1", "Cebpa"), OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = 'BP')
enrich_go_table = as.data.frame(enrich_ox_CM_go)

# enrichkegg
enrich_ox_CM_kegg = enrichKEGG(gene = gene_id, organism = "mmu", pvalueCutoff = 0.05)
enrich_kegg_table = as.data.frame(enrich_ox_CM_kegg)

pathways = enrich_go_table$Description

enrich_int = enrich_go_table %>% 
  filter(str_detect(pathways, pattern = "[Ll]ipid"))

# enrich barplot

pathway_int = c("mmu05417", "mmu04920", "mmu04071", "mmu04145")

pathway_int = enrich_kegg_table[pathway_int, ]

pathway_int_go = c("GO:0019216", "GO:0055088", "GO:0006869")

pathway_int = rbind(pathway_int, enrich_go_table[pathway_int_go, ])
pathway_int$Description = factor(pathway_int$Description, levels = c())

ggplot(pathway_int) + 
  geom_bar(aes(x = Description, y = -log10(qvalue)), fill = "lightblue", stat = "identity") + 
  coord_flip() + 
  theme_classic()

# GSVA analysis

pathway = read.gmt("m5.all.v2023.1.Mm.symbols.gmt")
pathway = read.gmt()

pathway = pathway %>% 
  pivot_wider(names_from = term, values_from = gene) %>% 
  as.list() %>% 
  map(1)

names(pathway)[names(pathway) %like% "TREM2"]

pathway_score = vector(mode = "list", length = length(pathway))

for(i in seq_along(pathway)){
  names(pathway_score)[i] = names(pathway)[i]
  pathway_score[[i]] = pathway[[i]][[1]]
}

counts_normalized = counts(dds, normalized = T)
GSVA_score = GSVA::gsva(expr = counts_normalized, gset.idx.list = pathway_score, parallel.sz = 10)

res_ox_CM_table = res_ox_CM_table %>% 
  arrange(desc(log2FoldChange))

res_sig = res_ox_CM_table %>% 
  filter(log2FoldChange > 1)

res_sig = rownames(res_sig)

# subset analysis of tumor-associated myeloid cells
# filter out normal samples and run clustering/annotation on myeloid cells in tumor samples only

# load libraries
library(Seurat)
library(tidyverse)
library(harmony)
library(MAST)
library(AnnoProbe)
library(ggpubr)
library(ggrepel)
library(sceasy)
library(ggsci)
'%notin%' = Negate('%in%')

# subset the myeloid cells
seu_anno = convertFormat('project/mice_landscape/sc_mice_annotation.h5ad', from = 'anndata', to = 'seurat', main_layer = 'counts')

seu_myeloid <- subset(seu_anno, annotation %in% c('Monocyte-linage') & condition == "PT")

# standard pipline
seu_myeloid <- seu_myeloid %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c('mt_ratio')) %>% 
  RunPCA()

seu_myeloid <- RunHarmony(seu_myeloid, group.by.vars = "sample_id", plot_convergence = T, max.iter.harmony = 30)

seu_myeloid <- RunUMAP(seu_myeloid, reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(res = c(0.8))

DimPlot(seu_myeloid, label = T, group.by = "RNA_snn_res.0.8") + 
  scale_color_igv()

Idents(seu_myeloid) <- "RNA_snn_res.0.8"

# Find markers
markers_myeloid <- FindAllMarkers(seu_myeloid, test.use = "wilcox", logfc.threshold = 0.5, min.pct = 0.25, only.pos = T)
marker2 = FindMarkers(seu_myeloid, ident.1 = 0, ident.2 = 2, logfc.threshold = 0.25)
marker2 = marker2 %>% 
  arrange(avg_log2FC)
markers_myeloid = markers_myeloid %>% 
  mutate(diff.pct = pct.1 - pct.2)
markers_top20 <- markers_myeloid %>% 
  filter(p_val_adj <= 0.05) %>% 
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n = 30)

# Annotation 

# Annotation

seu_myeanno <- RenameIdents(seu_myeloid, "0" = "RTM-like TAM",  # CD5L+ VCAM1+
                             "1" = "Inflam-TAM", 
                             "2" = "Cycling cells", 
                             "3" = "RTM-like TAM", 
                             "4" = "Folr2+ TAM",
                             "5" = "cDC2", # Clec10a+ DC
                             "6" = "Trem2+ Spp1+ LAM", 
                             "7" = "pDC", # Siglech+ Ccr9+
                             "8" = "cDC1", # Xcr1+ Wdfy4+ Irf8+
                             "10" = "cDC1")
DimPlot(seu_myeanno)

# GSEA analysis on scRNA-seq datasets
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

gene_to_label = c("Edn1", "Fabp3", "Ptges", "Furin", "Hrh2", "Ptgs2", "Il10", "Src", "Bnip3", "Tnfaip3", "Rcan1", "Cxcl9", "Ccl22", "Ccl5", "Il1rn")

row_label = rep("", length = length(gene_int))

names(row_label) = gene_int

row_label[gene_to_label] = gene_to_label

p = pheatmap(exp_mat[gene_int, ], cluster_rows=T, 
             cluster_cols=F, scale = "row", border_color = NA, color = rev(mycol), labels_row = row_label)

