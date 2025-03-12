library(Seurat)
library(tidyverse)
library(ggrepel)
library(viridis)
library(RColorBrewer)

# Test the correlation of SPP1 with other genes that are expressed in at least 10 percent of cells
seu_mye <- subset(seu_mye, sampletype %in% c("scRNA-seq"))
seu_mye = NormalizeData(seu_mye)
exp_mat = seu_mye@assays$RNA@data
exp <- t(as.matrix(seu_mye@assays$RNA@data))
saveRDS(exp, file = "exp_mat.rds")

## filter genes that are expressed in at least 10% cells
exp_genes <- colnames(exp)
cells <- nrow(exp)
type <- c()
for(i in 1:length(exp_genes)){
  gene_exp <- exp[,i]
  n0 <- sum(gene_exp == 0)
  if(n0/cells <= 0.9){
    type[i] <- "Keep"
  }
  else{
    type[i] <- "Filter"
  }
}
type <- ifelse(type == "Keep", T, F)
exp <- exp[,type]

## calculate the correlation

calculate_cor = function(exp_data, gene, type = "pearson", adjust.p = TRUE, ...){
  
  # the exp_data should be given as gene in column, sample in row
  target_exp = exp_data[,gene]
  all_gene = colnames(exp_data)
  all_gene = all_gene[all_gene != gene] # remove itself
  
  cor_df = data.frame(gene = all_gene)
  for(i in seq_along(all_gene)){
    test <- cor.test(as.numeric(exp_data[,all_gene[i]]), target_exp, type = type, )
    cor_df[i,2] <- test$estimate
    cor_df[i,3] <- test$p.value
  }
  colnames(cor_df)[c(2, 3)] = c("R", "P")
  
  if(adjust.p){
    cor_df$adjust_p = p.adjust(cor_df$P, ...)
  }
  return(cor_df)
}

cor_df = calculate_cor(exp_data = exp, gene = "SPP1", type = "pearson", adjust.p = T, method = "BH")
head(sort(cor_df$adjust_p, decreasing = F))

## remove ribosome genes and mt genes
cor_df = cor_df %>% 
  filter(str_detect(gene, "^RP[LS]", negate = T) & str_detect(gene, "^MT-", negate = T))

## set the significant correlation
cor_up <- cor_df %>% 
  filter(R > 0 & P < 0.01) %>% 
  arrange(desc(R))
cor_down <- cor_df %>% 
  filter(R < 0 & P < 0.01) %>% 
  arrange(R)
cor_df <- cor_df %>% 
  mutate(type = case_when(
    gene %in% cor_up$gene ~ "up", 
    gene %in% cor_down$gene ~ "down", 
    TRUE ~ "stable"
  ))

cor_int = rbind(cor_up[1:10,], cor_down[1:10,])
ggplot(cor_df) + 
  geom_point(aes(x = R, y = -log10(adjust_p), color = type), size = 2) + 
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) + 
  theme_bw() + 
  geom_text_repel(data = cor_int, aes(x = R, y = -log10(adjust_p), label = gene))

# heatmap showing specific correlation in different condition

# switch exp to contain all cells

condition = str_split(rownames(exp), pattern = "_") %>% 
  sapply(FUN = "[", 1) %>% 
  str_remove(pattern = "\\d")

exp_pnt = exp[condition == "PNT", ]
exp_pt = exp[condition != "PNT",]
exp_pt = exp_pt[sample(1:nrow(exp_pt), size = nrow(exp_pnt)), ]

gene_int = cor_up$gene[c(1:5, 7:11)]

cor_df_pnt = calculate_cor(exp_data = exp_pnt, gene = "SPP1", type = "pearson", adjust.p = T, method = "BH")
cor_df_pt = calculate_cor(exp_data = exp_pt, gene = "SPP1", type = "pearson", adjust.p = T, method = "BH")

cor_df_pnt_int = cor_df_pnt %>% 
  filter(gene %in% gene_int)
cor_df_pt_int = cor_df_pt %>% 
  filter(gene %in% gene_int)

mat = data.frame(R_pnt = cor_df_pnt_int$R, R_pt = cor_df_pt_int$R, row.names = cor_df_pnt_int$gene) %>% 
  arrange(desc(R_pt)) %>% 
  as.matrix()

ComplexHeatmap::pheatmap(mat = mat, cluster_rows = F, cluster_cols = F, 
                         color = rev(colorRampPalette(brewer.pal(10, name = "RdBu"))(50)))

                         