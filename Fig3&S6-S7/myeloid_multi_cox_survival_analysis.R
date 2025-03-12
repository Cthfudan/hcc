# progonosis analysis on each macrophage cluster

# load libraries
library(Seurat)
library(tidyverse)
library(forestploter)
library(GSVA)
library(SummarizedExperiment)

# perform DEG

Idents(myeloid_integ) = myeloid_integ$int_clusters

DEG_myeloid = FindAllMarkers(myeloid_integ, logfc.threshold = .1, test.use = 'MAST', only.pos = T)

# using first 50 genes as signature (rank by avg_log2FC)

DEG_top50 = DEG_myeloid %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50)

celltype = unique(DEG_top50$cluster)

gene_sig = vector(mode = 'list', length = length(celltype))

for(i in seq_along(celltype)){
  DEG = DEG_top50 %>% 
    filter(cluster == celltype[i])
  gene_sig[[i]] = DEG$gene
}

# score signature

gene_score = gsva(expr = exp_mat, gset.idx.list = gene_sig, method = 'gsva', kcdf = 'Poisson', parallel.sz = 15)

# calculate prognosis
library(survival)
library(survminer)

clinical_df <- meta_sample %>% 
  dplyr::select(patient, vital_status, days_to_death, days_to_last_follow_up, gender, age_at_index, ethnicity, ajcc_pathologic_stage) ## select needed information for survival analysis

clinical_df <- clinical_df %>% 
  mutate(deceased = case_when(
    vital_status == 'Alive' ~ FALSE, 
    TRUE ~ TRUE
  )) # create a var where a patient is alive is FALSE

clinical_df$overall_survival <- ifelse(clinical_df$deceased, clinical_df$days_to_death, clinical_df$days_to_last_follow_up) # define the overall survival time

clinical_df <- clinical_df %>% 
  filter(!is.na(overall_survival) & !is.na(ajcc_pathologic_stage))

gene_score = gene_score[,colnames(gene_score) %in% rownames(clinical_df)]

gene_score <- as.data.frame(t(gene_score)) %>% 
  mutate(futime = clinical_df$overall_survival, fustat = clinical_df$deceased, age = clinical_df$age_at_index, sex = clinical_df$gender, tumor_stage = clinical_df$ajcc_pathologic_stage)

gene_score = gene_score %>% 
  dplyr::select(futime, fustat, age, sex, tumor_stage, everything())

gene_score = gene_score %>% 
  mutate(tumor_stage = case_when(
    tumor_stage %in% c('Stage I') ~ 1, 
    tumor_stage %in% c('Stage II') ~ 2, 
    tumor_stage %in% c('Stage IIIA', 'Stage IIIB', 'Stage IIIC') ~ 3, 
    TRUE ~ 4
  )) %>% 
  mutate(sex = case_when(
    sex == 'male' ~ 1, 
    TRUE ~ 2
  ))

for(i in 1:11){
  score = gene_score[,paste('V', i, sep = '')]
  score_level = rep(2, length = nrow(gene_score))
  score_level[score >= quantile(score, .75)] = 3
  score_level[score <= quantile(score, .25)] = 1
  gene_score[,paste('V', i, sep = '')] = score_level
} # it can also be done by converting it to a factor
  
# multicox analysis

cox = coxph(Surv(futime, fustat) ~ age + sex + tumor_stage + V1 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = gene_score)
p = ggforest(model = cox, data = gene_score, fontsize = .5)

saveRDS(gene_score, file = 'data/Myeloid_cell/prognosis/gene_score_prognosis.rds')

# TREM2+ LAM signature prognosis ------------------------------------------

clinical_df <- meta_sample %>% 
  dplyr::select(patient, vital_status, days_to_death, days_to_last_follow_up, gender, age_at_index, ethnicity) ## select needed information for survival analysis

clinical_df <- clinical_df %>% 
  mutate(deceased = case_when(
    vital_status == 'Alive' ~ FALSE, 
    TRUE ~ TRUE
  ))

clinical_df$overall_survival <- ifelse(clinical_df$deceased, clinical_df$days_to_death, clinical_df$days_to_last_follow_up) # define the overall survival time

clinical_df <- clinical_df %>% 
  filter(!is.na(overall_survival))

gene_score = gsva(expr = exp_mat, gset.idx.list = gene_sig[8], method = 'gsva', kcdf = 'Poisson', parallel.sz = 15)

gene_score = gene_score[,colnames(gene_score) %in% rownames(clinical_df)]

gene_score <- as.data.frame(gene_score) %>% 
  mutate(futime = clinical_df$overall_survival, fustat = clinical_df$deceased, age = clinical_df$age_at_index, sex = clinical_df$gender)

gene_score = gene_score %>% 
  dplyr::select(futime, fustat, age, sex, everything())

res_cut = surv_cutpoint(data = gene_score, time = 'futime', event = 'fustat', variables = 'gene_score', minprop = 0.25)
summary(res_cut)

res_cat = surv_categorize(res_cut) # categorize the continuous variable

res_cat$futime = res_cat$futime/365

ggsurvplot(survfit(Surv(futime, fustat) ~ gene_score, data = res_cat), pval = TRUE, palette = c("#377EB8", "#E41A1C"), conf.int = F, risk.table = T, size = .25)


