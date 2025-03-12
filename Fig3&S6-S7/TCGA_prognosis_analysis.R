library(tidyverse)
library(survival)
library(survminer)

# expression

exp_mat <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/TCGA_validation/data/normalized_exp_mat.rds")


# survival 

clinical_data = readxl::read_excel("clinical/TCGA-CDR.xlsx", sheet = "TCGA-CDR")

clinical_data = clinical_data %>% 
  filter(type == "LIHC")

clinical_data = clinical_data %>% 
  filter(histological_type == "Hepatocellular Carcinoma")

# change the colnames

sample_id = colnames(exp_mat)

sample_id = str_split(sample_id, pattern = "-01A-")

sample_id = map(sample_id, 1) %>% unlist()

sample_id = str_split(sample_id, pattern = "-01B-")

sample_id = map(sample_id, 1) %>% unlist()

colnames(exp_mat) = sample_id

clinical_data = clinical_data %>% 
  filter(!is.na(PFI.time) & is.na(Redaction))

exp_mat = t(exp_mat) %>% as.data.frame()

exp_mat = as.data.frame(exp_mat)

exp_mat$colnames = rownames(exp_mat)

clinical_data = merge(clinical_data, exp_mat, by.x = "bcr_patient_barcode", by.y = "colnames")

# survival analysis

gene = "HMGCR" # test

cut = surv_cutpoint(data = clinical_data, time = 'PFI.time', event = 'PFI', variables = gene, minprop = .25)

res_cat = surv_categorize(cut) # categorize the continuous variable

ggsurvplot(survfit(Surv(PFI.time, PFI) ~ HMGCR, data = res_cat), pval = T, palette = c("#E64B35FF", "#4DBBD5FF"), conf.int = F, risk.table = T, )

# survival analysis on two genes

g1 = "TREM2"
g2 = "SPP1"

cut_g1 = surv_cutpoint(data = clinical_data, time = 'PFI.time', event = 'PFI', variables = g1)
cut_g1 = cut_g1$cutpoint$cutpoint

cut_g2 = surv_cutpoint(data = clinical_data, time = 'PFI.time', event = 'PFI', variables = g2)
cut_g2 = cut_g2$cutpoint$cutpoint

clinical_data = clinical_data %>% 
  dplyr::mutate(group = case_when(
    TREM2 > cut_g1 & SPP1 > cut_g2 ~ "g1hg2h", 
    TREM2 > cut_g1 & SPP1 <= cut_g2 ~ "g1hg2l", 
    TREM2 <= cut_g1 & SPP1 > cut_g2 ~ "g1lg2h", 
    TRUE ~ "g1lg2l"
  ))

ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group, data = clinical_data), pval = TRUE, conf.int = F, risk.table = T, size = .25, censor = F, palette = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))

clinical_data_p = clinical_data %>% 
  filter(group %in% c("g1hg2h", "g1lg2l"))

ggsurvplot(survfit(Surv(PFI.time, PFI) ~ group, data = clinical_data_p), pval = TRUE, conf.int = F, risk.table = T, size = .25, censor = F, palette = c("#E64B35", "#4DBBD5"))
