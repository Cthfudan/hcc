# pyscenic visualization

# load libraries
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(tidyverse)
library(ComplexHeatmap)
library(plotly)
library(BiocParallel)
library(data.table)
library(stats)
library(reshape2)
library(ggheatmap)
library(ggrepel)
library(cowplot)
library(viridis)
library(sceasy)

# read data
mye_scenic <- open_loom('aucell.loom')
regulons_mat <- get_regulons(mye_scenic, column.attr.name = 'Regulons')
regulons <- regulonsToGeneLists(regulons_mat) # extract regulons

regulonAUC <- get_regulons_AUC(mye_scenic, column.attr.name = 'RegulonsAUC')
regulonAUC_thres <- get_regulon_thresholds(mye_scenic)

# RSS analysis
seu_mye = convertFormat("~/R/Projects/snRNA_scRNA_hcc/project/mouse_myeloid/adata_filtered.h5ad", from = "anndata", to = "seurat", main_layer = "data")
meta <- seu_mye@meta.data %>% 
  select(sample_id, condition, nCount_RNA, nFeature_RNA, annotation2)

celltype <- meta %>% 
  select(annotation2)

rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = celltype[colnames(regulonAUC),])

rss <- na.omit(rss)

p = plotRSS(
  rss, 
  zThreshold = 1, 
  order_rows = T, 
  cluster_columns = F, 
  varName = 'annotaiton2'
)
p
## heatmap

rss_data <- p$plot$data

rss_data <- dcast(rss_data, Topic~rss_data$int_clusters, value.var = 'Z')

rownames(rss_data) <- rss_data[,1]

rss_data <- rss_data[,-1]
rss_data <- t(rss_data)

row_anno <- as.data.frame(rownames(rss_data))
colnames(row_anno) <- 'celltype'
rownames(row_anno) <- row_anno$celltype

pheatmap(mat = rss_data, color = colorRampPalette(c('#1A5592', 'white', '#B83D3D'))(100), cluster_rows = T, cluster_cols = T, scale = 'row', clustering_method = 'ward.D2')

saveRDS(rss, file = 'rss_result.rds')
saveRDS(rss_data, file = 'rss_mapped_data.rds')

## showing cell sepcific regulon rank
rss <- as.data.frame(rss)
table(celltype)
celltype <- c( "Trem2+ LAM")
rssRank <- list()

for(i in 1:length(celltype)){
  data_rank_plot <- rss %>% 
    select(celltype[i]) %>% 
    rownames_to_column(var = 'TF') %>% 
    na.omit()
  
  colnames(data_rank_plot) <- c('TF', 'celltype')
  
  data_rank_plot <- arrange(data_rank_plot, desc(celltype))
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))
  
  p <- ggplot(data_rank_plot, aes(x = rank, y = celltype)) + 
    geom_point(size = 1, shape = 16, color = '#1F77B4', alpha = 0.5) + 
    geom_point(data = data_rank_plot[c(1:6, 13), ], size = 1.5, color = '#DC050C') + # select the top6 TF to mark
    theme_bw() + 
    theme(
      axis.title = element_text(color = 'black', size = 12), 
      axis.text = element_text(color = 'black', size = 10), 
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank()
    ) + 
    labs(x = 'Regulon rank', y = 'Score', title = celltype[i]) + 
    geom_text_repel(data = data_rank_plot[c(1:6, 13),], aes(label = TF), color = 'black', size = 3, arrow = arrow(ends = 'first', length = unit(0.01, 'npc'))) # mark the TF names
  rssRank[[i]] <- p
}

plot_grid(rssRank[[1]])

# compare regulon activity of different conditions

## select celltype (we here analyze macrophage)
library(limma)
library(data.table)

anaAUC <- regulonAUC
anaAUC <- anaAUC[onlyNonDuplicatedExtended(rownames(anaAUC)),]

cellinfo <- meta %>% 
  filter(annotation2 %in% c("Flor2+ TAM", "Inflam-TAM", "RTM-like TAM", "Trem2+ LAM", "Cycling cells"))

group <- cellinfo %>% 
  select(annotation2) %>% 
  mutate(condition = case_when(
    annotation2 == 'Trem2+ LAM' ~ 'TREM2p', 
    TRUE ~ 'TREM2n'
  )) %>%
  rownames_to_column(var = 'barcode')

## Differential expression analysis
f <- factor(group$condition, levels = unique(sort(group$condition)))
design <- model.matrix(~0 + f)
colnames(design) <- c('TREM2n', 'TREM2p') # construct design matrix

cell_auc <- getAUC(anaAUC)
cell_auc <- cell_auc[,colnames(cell_auc)%in%group$barcode]
ncol(cell_auc)
cell_auc <- t(scale(t(cell_auc))) # get the auc in each cell and scale

contrast_mat <- makeContrasts('TREM2n-TREM2p', levels = design)
fit <- lmFit(cell_auc, design)
fit_con <- contrasts.fit(fit, contrasts = contrast_mat)
fit_con <- eBayes(fit_con)
diff_tf <- topTable(fit_con, adjust.method = 'BH', sort.by = 'logFC', n = Inf)
diff_tf <- diff_tf %>% 
  select(adj.P.Val, P.Value, logFC) %>% 
  dplyr::rename('FDR' = adj.P.Val, 'p_value' = P.Value)
write_csv(diff_tf, file = 'scenic_diffexp_TREM2-_vs_+.csv')

## filter the significant TFs
logfc_T <- 0.01 
P_T <- 0.05

diff_sig <- diff_tf %>% 
  rownames_to_column(var = 'TF') %>% 
  filter(abs(logFC) > logfc_T, FDR < P_T)

diff_sig$TF <- strsplit(diff_sig$TF, split = '(', fixed = T) %>% 
  sapply(FUN = '[', 1) # remove (+)

sig_tf_down <- diff_sig %>% 
  filter(logFC < 0)

sig_tf_up <- diff_sig %>% 
  filter(logFC > 0)

## volcano plot visualization
ggplot(diff_tf, aes(x = logFC, y = -log10(FDR))) + 
  geom_point(size = 3, shape = 16, color = 'grey') + 
  geom_point(data = sig_tf_up, size = 3, shape = 16, color = '#F8766D') + 
  geom_point(data = sig_tf_down, size = 3, shape = 16, color = '#00a9ff') + 
  geom_vline(xintercept = 0.5, linetype = 2) + 
  geom_vline(xintercept = -0.5, linetype = 2) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 300)) + 
  geom_text_repel(data = diff_sig, mapping = aes(label = TF), color = 'black', size = 4, arrow = arrow(ends = 'first', length = unit(0.01, 'npc')))

