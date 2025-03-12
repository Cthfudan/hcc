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
