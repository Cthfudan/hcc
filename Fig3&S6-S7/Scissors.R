# scissors prognosis analysis

library(Scissor)
library(Seurat)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(sceasy)
library(ggsci)

# prepare input dataset ---------------------------------------------------

# the input requires three dataset: the scRNA dataset, the bulk dataset and the survival data

## TCGA dataset

### data preprocessing

### TCGA
exp_data <- exp_data[, exp_data$primary_diagnosis == 'Hepatocellular carcinoma, NOS']
exp_data = exp_data[, exp_data$shortLetterCode == 'TP']
gene_anno <- as.data.frame(rowData(exp_data))
meta_sample = as.data.frame(colData(exp_data))

exp_mat = TCGAanalyze_Preprocessing(exp_data) # use unnormalized data because scissor will normalize RNA-seq data

### LIRI-JP
rownames(meta) = meta$icgc_sample_id
meta = meta %>% 
  filter(specimen_type == "Primary tumour - solid tissue")
exp_mat = exp_mat[,colnames(exp_mat) %in% meta$icgc_sample_id]

### remove the gene id version

rownames(exp_mat) = str_remove(rownames(exp_mat), pattern = '\\.[0123456789]+')

ensm_id = rownames(exp_mat)

### mapIDs

library(AnnotationDbi)
library(org.Hs.eg.db)

gene_name = mapIds(x = org.Hs.eg.db, keys = ensm_id, keytype = 'ENSEMBL', column = 'SYMBOL')
gene_name = gene_name[!is.na(gene_name)]

exp_mat = exp_mat[rownames(exp_mat) %in% names(gene_name), ]

rownames(exp_mat) = gene_name

## progonosis data

### TCGA
survival_info <- meta_sample %>% 
  dplyr::select(patient, vital_status, days_to_death, days_to_last_follow_up) ## select needed information for survival analysis

survival_info <- survival_info %>% 
  filter(vital_status != 'Not Reported') %>% 
  mutate(status = case_when(
    vital_status == 'Alive' ~ 0, 
    TRUE ~ 1
  )) # create a var where a patient is alive is FALSE

survival_info$overall_survival <- ifelse(survival_info$status == 1, survival_info$days_to_death, survival_info$days_to_last_follow_up) # define the overall survival time

survival_info = survival_info %>% 
  filter(!is.na(overall_survival))

## LIRI-JP
clinical_df <- clinical_df %>% 
  filter(!is.na(overall_survival))

survival_info <- meta_sample %>% 
  dplyr::select(icgc_sample_id, donor_vital_status, donor_survival_time)

survival_info = survival_info %>% 
  filter(donor_vital_status != 'Not Reported') %>% 
  mutate(status = case_when(
    donor_vital_status == 'alive' ~ 0, 
    TRUE ~ 1
  ))

survival_info = as.data.frame(survival_info) %>% 
  column_to_rownames(var = "icgc_sample_id")

survival_info = survival_info %>% 
  filter(!is.na(overall_survival))

### make the order of progonosis data match the exp_mat

exp_mat = exp_mat[,colnames(exp_mat) %in% rownames(survival_info)]

patient_order = colnames(exp_mat)

survival_info = survival_info[match(patient_order, rownames(survival_info)), ]

survival_info = survival_info %>% 
  dplyr::select(overall_survival, status) %>% 
  dplyr::rename('time' = overall_survival) %>% 
  as.matrix()

## scRNA-seq dataset

seu_mye = convertFormat("seu_mye.h5ad", from = "anndata", to = "seurat", main_layer = "counts")
seu_mye = convertFormat("data/GSE_mye_filtered_raw.h5ad", from = "anndata", to = "seurat", main_layer = "counts")
seu_mye = FindNeighbors(seu_mye, reduction = "scanvi", dims = 1:30)
seu_mye = NormalizeData(seu_mye)

# Run scissors ------------------------------------------------------------

infos1 <- Scissor(exp_mat, sc_dataset = seu_mye, phenotype = survival_info, alpha = 0.05, 
                  family = "cox", Save_file = 'Scissor_TCGA_survival_2.RData')
infos1 <- Scissor(exp_mat, sc_dataset = seu_mye, phenotype = survival_info, alpha = 0.05, 
                  family = "cox", Save_file = 'Scissor_LIRI_survival.RData')
infos1$Scissor_pos

## tune the model
infos2 <- Scissor(exp_mat, sc_dataset = seu_mye, phenotype = survival_info, alpha = 0.01, 
                  family = "cox", Load_file = 'Scissor_TCGA_survival_2.RData')

# add scissors results ----------------------------------------------------

Scissor_select <- rep(0, ncol(seu_mye))
names(Scissor_select) <- colnames(seu_mye)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
seu_mye <- AddMetaData(seu_mye, metadata = Scissor_select, col.name = "scissor_0.01")
DimPlot(seu_mye, reduction = 'umap', group.by = 'scissor_0.01', cols = c('grey','indianred1','royalblue'), pt.size = .1, order = c(2,1)) # size: 6.7*6.2


# Visualization -----------------------------------------------------------

## demonstrate the composition of scissor+ cells

meta = seu_mye@meta.data

cell_com = meta %>% 
  dplyr::select(celltype, scissor_0.01) %>% 
  filter(scissor_0.01 == 1)

cell_prop = cell_com %>% 
  group_by(celltype) %>% 
  summarise(sum = n(), prop = sum/nrow(cell_com)) %>% 
  arrange(desc(prop))

cell_prop$celltype = factor(cell_prop$celltype, levels = rev(c('Classical monocyte', 'FOLR2+ TAM', 'IFN-TAM', 'Inflam-TAM', 'TREM2+ TAM', 'Cycling cells', 'cDC1', 'cDC2')))

ggplot(cell_prop, aes(x = celltype, y = prop)) + 
  geom_col(aes(fill = celltype)) + 
  theme_classic() + 
  scale_fill_manual(values = rev(c('#378C4F', '#F5CDCD', '#D9579B', '#A59ACB', '#006DDBFF',  '#6DB6FFFF', '#B6DBFFFF', '#E2A7CC'))) + 
  labs(x = '', y = 'Fraction of cells') + 
  labs(fill = 'Cell types') + 
  coord_flip()

## condition visualization of scissor+ cells

cell_com_con = meta %>% 
  dplyr::select(site, scissor_0.01) %>% 
  mutate(condition = case_when(
    site == 'Normal' ~ 'PNT', 
    TRUE ~ 'PT'
  )) %>% 
  filter(scissor_0.01 == 1)

cell_prop = cell_com_con %>% 
  group_by(condition) %>% 
  summarise(sum = n(), prop = sum/nrow(cell_com_con)) %>% 
  arrange(desc(prop))

ggplot(cell_prop, aes(x = factor(condition, levels = c('PT', 'PNT')), y = prop)) + 
  geom_col(aes(fill = condition)) + 
  theme_classic() + 
  scale_fill_npg() + 
  labs(x = '', y = 'Fraction of cells') + 
  labs(fill = 'Cell types') + 
  coord_flip() + 
  scale_y_continuous(limits = c(0, 1))

## pie chart

ggplot(cell_prop, aes(x = '', y = prop, fill = int_clusters)) + 
  geom_bar(stat = 'identity', width = 1) + 
  theme_classic() + 
  scale_fill_manual(values = c('#6DB6FFFF', '#009292FF', '#FF6DB6FF', '#490092FF', '#004949FF', '#FFB6DBFF', '#920000FF', '#B6DBFFFF', '#B66DFFFF', '#924900FF', '#006DDBFF')) + 
  coord_polar(theta = 'y', start = 0, direction = 1)

## DEG between scissors+ and scissors- cells

library(MAST)
library(ggrepel)
meta = seu_mye@meta.data

meta = meta %>% 
  mutate(scissors_DEG = case_when(
    scissor_0.01 == 1 ~ 'Scissor +', 
    TRUE ~ 'Others'
  ))

seu_mye@meta.data = meta

DEG_scissors = FindMarkers(object = seu_mye, ident.1 = 'Scissor +', group.by = 'scissors_DEG', test.use = 'MAST', logfc.threshold = .01, min.pct = .1)

## modified volcano plot

plot_volcano_enhanced = function(DEG_result, logFC_t, p_t, logFC_label_t = 1, pct.1 = pct.1, pct.2 = pct.2, palette = c('royalblue', 'grey', 'indianred1'), logFC_col = avg_log2FC, p_col = p_val_adj, return_data = F){
  
  DEG = DEG_result %>% 
    mutate(pct.diff = {{pct.1}} - {{pct.2}})
  
  DEG = DEG %>% 
    mutate(signif = case_when(
      {{logFC_col}} > logFC_t & {{p_col}} < p_t ~ 'Up', 
      {{logFC_col}} < -logFC_t & {{p_col}} < p_t ~ 'Down', 
      TRUE ~ 'stable'
    ))
  
  DEG_sig = DEG %>% 
    filter(signif %in% c('Up', 'Down')) %>% 
    filter(abs({{logFC_col}}) > logFC_label_t)
  
  Trem2 = DEG["TREM2", ]
  DEG_sig = rbind(DEG_sig, Trem2)
  
  p = ggplot(DEG, aes(x = pct.diff, y = {{logFC_col}})) + 
    geom_point(aes(color = signif), size = .5) + 
    scale_color_manual(values = c('royalblue', 'grey', 'indianred1')) + 
    geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.8) + 
    geom_hline(yintercept = 0, lty = 4, col = "black", lwd = 0.8) + 
    theme_bw() + 
    ggrepel::geom_text_repel(data = DEG_sig, aes(label = rownames(DEG_sig)), max.overlaps = 10)
  
  if(return_data){
    return(list(DEG, p))
  }
  
  return(p)
}
"%notin%" = Negate("%in%")
DEG_scissors = DEG_scissors %>% 
  mutate(diff.pct = pct.1 - pct.2) %>% 
  arrange(desc(diff.pct))

plot_volcano_enhanced(DEG_result = DEG_scissors, logFC_t = 0, p_t = 0.01)

## survival analysis on Scissors+ signature genes -- using the top50 most upregulated genes

DEG_scissors_top50 = DEG_scissors %>% 
  filter(avg_log2FC > 0) %>% 
  slice_max(order_by = avg_log2FC, n = 50)
gene_sig = list(c(DEG = rownames(DEG_scissors_top50)))

exp_mat = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/TCGA_validation/data/normalized_exp_mat.rds")
gene_score = GSVA::gsva(expr = exp_mat, gset.idx.list = gene_sig, method = 'gsva', kcdf = "Poisson", parallel.sz = 20)
gene_score = GSVA::gsva(expr = exp_mat, gset.idx.list = gene_sig, method = 'ssgsea', parallel.sz = 20)
gene_score = GSVA::gsva(expr = exp_mat, gset.idx.list = gene_sig, method = 'zscore', parallel.sz = 20)

library(survival)
library(survminer)

clinical_df <- meta_sample %>% 
  dplyr::select(patient, vital_status, days_to_death, days_to_last_follow_up, gender, age_at_index, ethnicity)
clinical_df <- clinical_df %>% 
  mutate(deceased = case_when(
    vital_status == 'Alive' ~ FALSE, 
    TRUE ~ TRUE
  )) 
clinical_df$overall_survival <- ifelse(clinical_df$deceased, clinical_df$days_to_death, clinical_df$days_to_last_follow_up)
clinical_df <- clinical_df %>% 
  filter(!is.na(overall_survival))

gene_score = gene_score[,colnames(gene_score) %in% rownames(clinical_df)]
gene_score <- as.data.frame(gene_score) %>% 
  mutate(futime = clinical_df$overall_survival, fustat = clinical_df$deceased, age = clinical_df$age_at_index, sex = clinical_df$gender)
gene_score = gene_score %>% 
  dplyr::select(futime, fustat, age, sex, everything())

res_cut = surv_cutpoint(data = gene_score, time = 'futime', event = 'fustat', variables = 'gene_score', minprop = .25)
summary(res_cut)
res_cat = surv_categorize(res_cut) # categorize the continuous variable
ggsurvplot(survfit(Surv(futime, fustat) ~ gene_score, data = res_cat), pval = TRUE, palette = rev(c("#377EB8", "#E41A1C")), conf.int = F, risk.table = T, size = .25, censor.size = 0)


# LIRI-JP datasets ------------------------

# create a summerizedexperiment obj

meta_sample = meta_sample %>% 
  column_to_rownames(var = 'icgc_sample_id')

meta_sample = meta_sample[colnames(exp_mat), ]

exp_data = SummarizedExperiment(exp_mat, colData = meta_sample)

# data preparation

## expression data

exp_data <- exp_data[, exp_data$specimen_type == 'Primary tumour - solid tissue']

exp_mat = assay(exp_data)

meta_sample = as.data.frame(colData(exp_data))

## prognosis data

survival_info <- meta_sample %>% 
  dplyr::select(donor_vital_status, donor_survival_time)

survival_info = survival_info %>% 
  mutate(status = case_when(
    donor_vital_status == 'alive' ~ 0, 
    TRUE ~ 1
  ))

survival_info = survival_info %>% 
  filter(!is.na(donor_survival_time)) %>% 
  filter(rownames(survival_info) %in% colnames(exp_data))

### make the order of progonosis data match the exp_mat

patient_order = colnames(exp_mat)

survival_info = survival_info[match(patient_order, rownames(survival_info)), ]

survival_info = survival_info %>% dplyr::select(donor_survival_time, status) %>% 
  dplyr::rename('time' = donor_survival_time) %>% 
  as.matrix()

## scRNA-seq dataset

myeloid_integ = FindNeighbors(myeloid_integ, reduction = 'scanvi', dims = 1:30)

## Run scissors

infos1 <- Scissor(exp_mat, sc_dataset = myeloid_integ, phenotype = survival_info, alpha = 0.05, 
                  family = "cox", Save_file = 'data/Myeloid_cell/Scissor/Scissor_LIRI_survival.RData')

## tune the model
infos2 <- Scissor(exp_mat, sc_dataset = myeloid_integ, phenotype = survival_info, alpha = 0.008, 
                  family = "cox", Load_file = 'data/Myeloid_cell/Scissor/Scissor_LIRI_survival.RData')

Scissor_select <- rep(0, ncol(myeloid_integ))
names(Scissor_select) <- colnames(myeloid_integ)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
myeloid_integ <- AddMetaData(myeloid_integ, metadata = Scissor_select, col.name = "scissor_liri_0.008")
DimPlot(myeloid_integ, reduction = 'umap', group.by = 'scissor_liri_0.008', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
table(myeloid_integ$int_clusters, myeloid_integ$scissor_0.025)

saveRDS(myeloid_integ, file = 'data/Myeloid_cell/Scissor/scissor_seurat.rds')

## visualization

meta = myeloid_integ@meta.data

cell_com = meta %>% 
  dplyr::select(int_clusters, scissor_liri_0.008) %>% 
  filter(scissor_liri_0.008 == 1)

cell_prop = cell_com %>% 
  group_by(int_clusters) %>% 
  summarise(sum = n(), prop = sum/nrow(cell_com)) %>% 
  arrange(desc(prop))

cell_prop$int_clusters = factor(cell_prop$int_clusters, levels = unique(cell_prop$int_clusters))

ggplot(cell_prop, aes(x = int_clusters, y = prop)) + 
  geom_col(aes(fill = int_clusters)) + 
  theme_classic() + 
  scale_fill_manual(values = c('#6DB6FFFF', '#009292FF', '#FF6DB6FF', '#490092FF', '#004949FF', '#FFB6DBFF', '#920000FF', '#B6DBFFFF', '#B66DFFFF', '#924900FF', '#006DDBFF')) + 
  labs(x = '', y = 'Cell proportion', fill = 'Cell type')

markers = FindAllMarkers(seu)
