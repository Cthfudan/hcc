# extract metabolism-tumor specific genes

library(Seurat)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(SummarizedExperiment)
seu = NormalizeData(seu)
seu = subset(seu, archtype != "unassigned")

Idents(seu) = seu$archtype

DotPlot(seu, features = "SQLE")

markers = FindMarkers(seu, ident.1 = "archtype1", only.pos = T, min.pct = 0, logfc.threshold = .05)

markers = markers %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))

# data prep

##################################
#TCGA
gene_int = rownames(markers)

expr_mat = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/TCGA_validation/data/normalized_exp_mat.rds")

clinical_data = readxl::read_excel("~/R/Projects/snRNA_scRNA_hcc/project/TCGA_validation/clinical/TCGA-CDR.xlsx", sheet = "TCGA-CDR")

sample_id = colnames(expr_mat)

sample_id = str_split(sample_id, pattern = "-01A-")

sample_id = map(sample_id, 1) %>% unlist()

sample_id = str_split(sample_id, pattern = "-01B-")

sample_id = map(sample_id, 1) %>% unlist()

colnames(expr_mat) = sample_id

clinical_data = clinical_data %>% 
  filter(!is.na(PFI.time) & is.na(Redaction))

expr_mat = t(expr_mat) %>% as.data.frame()

expr_mat = as.data.frame(expr_mat)

expr_mat$colnames = rownames(expr_mat)

clinical_data = merge(clinical_data, expr_mat, by.x = "bcr_patient_barcode", by.y = "colnames")

clinical_data = clinical_data[, !is.na(colnames(clinical_data))]


##########################################################################

# LIRI-JP

library(SummarizedExperiment)

expr_mat = assay(exp_liri)
clinical_data = metadata(exp_liri)[[1]]
annotation = colData(liri_se) %>% as.data.frame()

annotation = annotation %>% 
  filter(submitted_sample_id %like% "_Cancer$") %>% 
  rownames_to_column(var = "sample_id") %>% 
  select(sample_id, icgc_donor_id)

clinical_data = clinical_data %>% 
  merge(annotation, by = "icgc_donor_id")

expr_mat = expr_mat[, colnames(expr_mat) %in% clinical_data$sample_id] %>% 
  t() %>% as.data.frame() %>% rownames_to_column(var = "sample_id") %>% 
  merge(clinical_data, by = "sample_id")

expr_mat = expr_mat %>% 
  mutate(donor_vital_status = case_when(
    donor_vital_status == "alive" ~ 0, 
    TRUE ~ 1
  ))

expr_mat$donor_survival_time = as.numeric(expr_mat$donor_survival_time)

###########################################################################

# GSE14520
library(SummarizedExperiment)
exp_data <- readRDS("~/R/Projects/recurrence_gset/data/GSE14520/exp_data.rds")
expr_mat = assay(exp_data)
clinical = fread("~/R/public_bulk_datasets/Microarray/GSE14520/HCCDB6.patient.txt", header = F)
clinical = as.data.frame(clinical) %>% 
  column_to_rownames(var = "V1")
clinical = t(clinical) %>% as.data.frame()
annotation = fread("~/R/public_bulk_datasets/Microarray/GSE14520/HCCDB6.sample.txt", header = F)
annotation = as.data.frame(annotation) %>% 
  column_to_rownames(var = "V1")
annotation = t(annotation) %>% as.data.frame()

annotation = annotation %>% 
  filter(TYPE == "HCC") %>% 
  select(SAMPLE_ID, PATIENT_ID)

clinical = clinical %>% 
  mutate(STATUS = case_when(
    STATUS == "Alive" ~ 0, 
    TRUE ~ 1
  ))

clinical = clinical %>% select(-PATIENT) %>% 
  merge(annotation, by = "PATIENT_ID")

clinical$SAMPLE_ID = str_replace(clinical$SAMPLE_ID, pattern = "-", replacement = ".")

expr_mat = t(expr_mat) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SAMPLE_ID") %>% 
  merge(clinical, by = "SAMPLE_ID")

expr_mat$SURVIVAL_TIME = as.numeric(expr_mat$SURVIVAL_TIME)

expr_mat = expr_mat %>% 
  filter(!is.na(SURVIVAL_TIME))

###############################################################
# gse124751

expr_mat = fread("~/R/public_bulk_datasets/GSE124751/HCCDB19_mRNA_level3.txt")

expr_mat = as.data.frame(expr_mat) %>% 
  select(-1) %>% 
  filter(!duplicated(Symbol)) %>% 
  column_to_rownames(var = "Symbol")

clinical = fread("~/R/public_bulk_datasets/GSE124751/HCCDB19.patient.txt") %>% as.data.frame()

annotation = fread("~/R/public_bulk_datasets/GSE124751/HCCDB19.sample.txt", header = F)

annotation = annotation %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "V1") %>% 
  t() %>% as.data.frame()

annotation = annotation %>% 
  select(PATIENT_ID, SAMPLE_ID) %>% 
  merge(clinical, by = "PATIENT_ID")

expr_mat = expr_mat %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SAMPLE_ID") %>% 
  merge(annotation, by = "SAMPLE_ID")

expr_mat = mutate(expr_mat, OS_STATUS = case_when(
  OS_STATUS == "Alive" ~ 0, 
  TRUE ~ 1
)) %>% 
  mutate(OS = as.numeric(OS))

expr_mat = mutate(expr_mat, DFS_STATUS = case_when(
  DFS_STATUS == "Yes" ~ 0, 
  TRUE ~ 1
)) %>% 
  mutate(DFS = as.numeric(DFS))

##########################################################
# OEP000321
expr_mat = fread("~/R/public_bulk_datasets/OEP000321/HCCDB25_mRNA_level3.txt")

expr_mat = as.data.frame(expr_mat) %>% 
  select(-1) %>% 
  column_to_rownames(var = "Symbol")

clinical = fread("~/R/public_bulk_datasets/OEP000321/HCCDB25.patient.txt") %>% as.data.frame()

annotation = fread("~/R/public_bulk_datasets/OEP000321/HCCDB25.sample.txt", header = F)

annotation = annotation %>% 
  as.data.frame() %>% 
  column_to_rownames(var = "V1") %>% 
  t() %>% as.data.frame()

annotation = annotation %>% 
  filter(TYPE == "HCC") %>% 
  select(SAMPLE_ID, PATIENT_ID) %>% 
  merge(clinical, by = "PATIENT_ID")

expr_mat = as.data.frame(expr_mat) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "SAMPLE_ID")

expr_mat = expr_mat %>% 
  merge(annotation, by = "SAMPLE_ID")

expr_mat = expr_mat %>% 
  mutate(OS_STATUS = case_when(
    OS_STATUS == "Alive" ~ 0, 
    TRUE ~ 1
  )) %>% 
  mutate(OS = as.numeric(OS))

### pipline

calculate_gene_hrp = function(expr_clinical, gene_int, time_column, status_column, minprop = .25){
  gene_surv = list()
  
  for(g in gene_int){
    
    if(! g %in% colnames(expr_clinical)){
      next
    }
    
    cut_g = tryCatch(
      {surv_cutpoint(data = expr_clinical, 
                     time = as.character(rlang::enexpr(time_column)), 
                     event = as.character(rlang::enexpr(status_column)), 
                     variables = g, 
                     minprop = minprop)}, 
      error = {function(err){
        return(NA)
      }})
    
    if(is.na(cut_g[[1]][[1]])){
      next
    }
    
    cut_g = cut_g$cutpoint$cutpoint
    
    expr_g = expr_clinical %>% 
      mutate(group = case_when(
        !!rlang::sym(g) > cut_g ~ 2, 
        TRUE ~ 1
      ))
    
    params = list(fun.time = substitute(time_column),
                  fun.event = substitute(status_column), 
                  fun.group = substitute(group),
                  fun.data = substitute(expr_g))
    expr = substitute(coxph(Surv(time = fun.time, event = fun.event) ~ fun.group, 
                               data = fun.data), params)
    
    cox_model = eval.parent(expr)
    
    HR = round(exp(coef(cox_model)), 2)
    
    p = summary(cox_model)$coefficients[, 5]
    
    gene_surv[[g]] = c(HR, p)
  }
  
  return(gene_surv)
}

gene_surv = calculate_gene_hrp(expr_clinical = expr_mat, gene_int = gene_int, time_column = donor_survival_time, status_column = donor_vital_status)

gene_surv = as.data.frame(gene_surv)
gene_surv = t.data.frame(gene_surv)
gene_surv = as.data.frame(gene_surv)

colnames(gene_surv) = c("HR", "P")

gene_surv = arrange(gene_surv, desc(HR), P)

gene_surv = gene_surv %>% 
  filter(HR > 1)

gene_surv["SQLE", ]

saveRDS(gene_surv, file = "surv_TCGA.rds")
saveRDS(gene_surv, file = "surv_LIRI.rds")
saveRDS(gene_surv, file = "surv_gse14520.rds")
saveRDS(gene_surv, file = "surv_gse124751.rds")
saveRDS(gene_surv, file = "surv_oep000321.rds")

gene = "SQLE"

cut = surv_cutpoint(data = expr_mat, time = 'OS_YEAR', event = 'OS_STATUS', variables = gene)

res_cat = surv_categorize(cut) # categorize the continuous variable

ggsurvplot(survfit(Surv(OS_YEAR, OS_STATUS) ~ SQLE, data = res_cat), pval = T, palette = c("#E64B35FF", "#4DBBD5FF"), conf.int = F, risk.table = T, )
