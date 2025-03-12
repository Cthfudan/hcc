# visualization heatmap

library(tidyverse)
library(pheatmap)
library(clusterProfiler)

surv_list = list(surv_gse124751, surv_gse14520, surv_LIRI, surv_oep000321, surv_TCGA)

# filter to metabolism-associated genes

kegg_metabolic = read.gmt("~/R/Projects/snRNA_scRNA_hcc/project/Myeloid_cell/metabolic_clustering/KEGG_metabolism_nc.gmt")
kegg_gene = kegg_metabolic$gene
kegg_gene = str_remove(kegg_gene, pattern = fixed("</td>"))
kegg_gene = kegg_gene[kegg_gene != ""]
metabolic_gene = kegg_gene
metabolic_gene = metabolic_gene[!duplicated(metabolic_gene)]
del = str_subset(metabolic_gene, pattern = "-")
metabolic_gene = metabolic_gene[! metabolic_gene %in% del]

surv_list = map(surv_list, function(x) filter(x, rownames(x) %in% metabolic_gene))
surv_list = map(surv_list, function(x) filter(x, rownames(x) %in% gene_up))


# define a group

surv_list = map(surv_list, function(x){
  x = x %>% 
    mutate(group = case_when(
      HR > 1 & P < 0.001 ~ 8, 
      HR > 1 & P < 0.05 & P >=0.001 ~ 7, 
      HR > 1 & P >= 0.05 & P < 0.1 ~ 6, 
      HR > 1 & P >= 0.1 ~ 5, 
      HR <= 1 & P >= 0.1 ~ 4, 
      HR <= 1 & P >= 0.05 & P < 0.1 ~ 3, 
      HR <= 1 & P < 0.05 & P >=0.001 ~ 2, 
      HR <= 1 & P < 0.001 ~ 1
    ))
  return(x)
})

surv_list = map(surv_list, function(x) rownames_to_column(x, var = "gene"))

surv_list = imap(surv_list, function(x, i){
  colnames(x)[2:4] = str_c(colnames(x)[2:4], as.character(i))
  return(x)
})

surv_all = purrr::reduce(.x = surv_list, .f = merge, by = "gene")

surv_all$group_sum = surv_all$group1 + surv_all$group2 + surv_all$group3 + surv_all$group4 + surv_all$group5

surv_all = arrange(surv_all, desc(group_sum))

# remove house-keeping genes

HK_genes = readxl::read_excel("Supplementary_Table1.xlsx", sheet = "Common gene", skip = 1)
HK_genes2 = read.gmt("HSIAO_HOUSEKEEPING_GENES.v2023.2.Hs.gmt")
surv_all = filter(surv_all, ! gene %in% HK_genes$Gene)
surv_all = filter(surv_all, ! gene %in% HK_genes2$gene)
# plot heatmap

surv_plot = surv_all[, c("gene", "group1", "group2", "group3", "group4", "group5")]

surv_plot = surv_plot %>% 
  column_to_rownames(var = "gene")
library(paletteer)
pheatmap(surv_plot[1:20, ], cluster_cols = F, cluster_rows = F, color = rev(c("#8E2261FF", "#E377AEFF", "#FFCFE7FF", "#FFF0F7FF", "#E5F6E5FF", "#A6DCA6FF", "#52AD53FF", "#296E3FFF")), border_color = NA)
