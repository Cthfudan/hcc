# program expression correlation

library(tidyverse)
library(corrplot)
library(Seurat)
library(pagoda2)
library(RColorBrewer)
library(BBmisc)

exp_mat = snRNA_tumor@assays$RNA@data

exp_mat = as.matrix(exp_mat)

## or

exp_mat <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/NMF/data/exp_mat_sct.rds")

score_meta = map(metaGene, ~ score.cells.puram(data = t(exp_mat), signature = .x))

score_meta = as.data.frame(score_meta)

colnames(score_meta) = str_c("sig", 1:9)

saveRDS(score_meta, file = "sig_score_meta.rds")
saveRDS(score_meta, file = "sig_score_meta_sct.rds")

score_meta_2 = score_meta %>% 
  select(-9)

score_cor = cor(score_meta, method = "pearson")
score_cor_2 = cor(score_meta_2, method = "pearson")

score_cor_2 = score_cor_2[c(2, 1, 3, 4, 5, 6, 7, 8), c(2, 1, 3, 4, 5, 6, 7, 8)]

corrplot(corr = score_cor_2, 
         method = "color", 
         type = "lower", 
         order = "original", 
         addCoef.col = "black",
         tl.srt = 45,
         tl.col = "black", 
         is.corr = T,
         hclust.method = "ward.D2", 
         col = colorRampPalette(rev(brewer.pal(name = 'RdBu', 11)))(50))

pheatmap::pheatmap(score_cor_2)

corrplot(corr = score_cor_2, method = "square", type = "lower")

# program expression on individual cells

score_meta_PT1 = score_meta_2 #[str_detect(rownames(score_meta_2), pattern = "^PT6_"), ]

score_meta_PT1 = as.data.frame(score_meta_PT1) %>% 
  arrange(sig4) %>% 
  select(sig4, sig2, sig1, sig3, sig6, sig7, sig8, sig5)

score_PT1 = map(score_meta_PT1, function(x){
  y = scale(x)
  y[y > quantile(y, .9)] = quantile(y, .9)
  y[y < quantile(y, .1)] = quantile(y, .1)
  scale_norm = BBmisc::normalize(as.numeric(y), method = "range", range = c(0, 1))
  return(scale_norm)
})

idx = rownames(score_meta_PT1)

score_PT1 = as.data.frame(score_PT1)

rownames(score_PT1) = idx

score_PT1 = score_PT1 %>% 
  arrange(sig4)

col_anno = data.frame(row.names = rownames(score_PT1), patient = str_extract(rownames(score_PT1), pattern = "^[PR]T[0-9]"))

col_anno$patient = factor(col_anno$patient, levels = c(str_c("PT", 1:9), str_c("RT", 1:3)))

anno_color = as.character(paletteer::paletteer_d("colorBlindness::paletteMartin"))[2:13]
names(anno_color) = unique(col_anno$patient) 

pheatmap::pheatmap(t(score_PT1), cluster_rows = F, cluster_cols = F, show_colnames = F, show_rownames = T, 
                   annotation_col = col_anno, annotation_colors = list(patient = anno_color))

# smooth curve plot

score_curve = mutate(score_meta_PT1, rank = 1:nrow(score_meta_PT1))

ggplot(score_curve, aes(x = rank, y = sig1)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig2)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig3)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig4)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig5)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig6)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()
 
ggplot(score_curve, aes(x = rank, y = sig7)) + 
  geom_smooth(method = , se = F) + 
  theme_classic()

ggplot(score_curve, aes(x = rank, y = sig8)) + 
  geom_smooth(method = "loess", se = F) + 
  theme_classic()
