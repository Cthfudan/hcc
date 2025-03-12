library(BayesPrism)
library(data.table)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(AnnoProbe)
library(Seurat)

'%notin%' = Negate('%in%')
# bulk RNA-seq data
exp_mat = fread("GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt") ## TCGA
exp_mat = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/LIRI-JP/exp_mat_LIRI-JP.rds") ## LIRI-JP

###################################################################################
# TCGA-LIHC 
exp_mat = as.data.frame(exp_mat) %>% 
  column_to_rownames(var = "V1")

exp_data <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/TCGA_validation/data/TCGA_RNAseq_data.rds")

sample_anno = as.data.frame(colData(exp_data))

sample_anno = sample_anno %>% 
  filter(primary_diagnosis == "Hepatocellular carcinoma, NOS" & sample_type == "Primary Tumor")

exp_mat = exp_mat[, colnames(exp_mat) %in% rownames(sample_anno)]

saveRDS(exp_mat, file = "TCGA_HCC_exp_counts.rds")
####################################################################################
# LIRI-JP
meta = readRDS("~/R/Projects/snRNA_scRNA_hcc/project/LIRI-JP/LIRI-JP_metadata.rds")

meta = meta %>% 
  filter(specimen_type == "Primary tumour - solid tissue")

exp_mat = exp_mat[, colnames(exp_mat) %in% meta$icgc_sample_id]

saveRDS(exp_mat, file = "LIRIJP_HCC_exp_counts.rds")
##################################################################################

bk.dat = t(exp_mat)

# single-cell RNA-seq data

seu_ref <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Spatial/data/scRNA_ref_high_res_new.rds")

seu_ref = subset(seu_ref, condition %in% c("PT", "RT", "tumor"))

seu_ref = subset(seu_ref, celltype %notin% c("apCAF", "Proliferating EC", "Proliferating T")) # clusters with too few cells are discarded

sc.dat = as.matrix(seu_ref@assays$RNA@counts) %>% t()

meta = seu_ref@meta.data

meta$cell_type_label = meta$celltype

meta = meta %>% 
  mutate(cell_type_label = case_when(
    new_celltype == "B cells" ~ "B cells", 
    new_celltype == "CAF" ~ "CAF", 
    new_celltype == "Endothelial cells" ~ "Endothelial cells", 
    new_celltype == "Myeloid cells" ~ "Myeloid cells", 
    new_celltype == "NK cells" ~ "NK cells", 
    new_celltype == "T cells" ~ "T cells", 
    new_celltype %in% c("archtype1", "archtype2", "inflammatory", "Proliferation") ~ "Malignant", 
    new_celltype == "Plasma cells" ~ "Plasma cells", 
  ))

table(meta$cell_type_label)

cell.type.labels = meta$celltype

cell.state.labels = meta$celltype

# QC on cell type and cell state labels
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell state correlation",
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))

plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)

# filter genes with low specifity on scRNA-seq and bulk dataset

sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

head(bk.stat)

## do the filtering
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

# check the concordance of gene expression in bulk and scRNA-seq dataset

plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

# construct a prism object

myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key="Hepatocytes",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# run bayes prism

bp.res <- run.prism(prism = myPrism, n.cores=16)

saveRDS(bp.res, file = "bp_result_sc_high_res_liri.rds")

slotNames(bp.res)

theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")

theta_state = get.fraction (bp=bp.res,
                            which.theta="first",
                            state.or.type="state")