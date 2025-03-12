library(sceasy)
library(tidyverse)
library(Seurat)
'%notin%' = Negate('%in%')

mouse_fname = "~/R/Projects/snRNA_scRNA_hcc/project/mouse_landscape/sc_mice_annotation.h5ad"

sc_mouse = convertFormat(mouse_fname, from = 'anndata', to = 'seurat', main_layer = 'counts')

gene_mouse = rownames(sc_mouse)

sc_human <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/landscape_QC/seurat_integration.rds")

# convert mouse gene to human conterparts ---------------------------------

mouse_human_genes = read.csv("HOM_MouseHumanSequence.rpt.txt",sep="\t")

convert_mouse_to_human <- function(gene_list){
  require('tidyverse')
  output = c()
  mouse_output = c()
  
  mouse_genes = mouse_human_genes %>% 
    filter(Common.Organism.Name == 'mouse, laboratory')
  
  human_genes = mouse_human_genes %>% 
    filter(Common.Organism.Name == 'human')
  
  for(gene in gene_list){
    
    class_key = (mouse_genes %>% filter(Symbol == gene))[['DB.Class.Key']]
    
    if(length(class_key) > 1){
      print(paste0(gene, " is a gene with multiple "))
      next
    }
    
    if(!identical(class_key, integer(0))){
      
      human_gene = (human_genes %>% filter(DB.Class.Key == class_key))[,"Symbol"]
      
      # genes without mapping or with multiple mappings were discarded 
      
      if((length(human_gene) > 1) | (length(human_gene) == 0)){ 
        next
      }
      
      output = append(output, human_gene)
      mouse_output = append(mouse_output, gene)
      names(output) = mouse_output # returns a one-one matching between mouse and human genes
      
      }
  }
  return(output)
}

human_mouse_con = convert_mouse_to_human(gene_list = gene_mouse)

human_mouse_con = human_mouse_con[!duplicated(human_mouse_con)]

# remove cellcycle genes to avoid biases

cell_cycle = cc.genes.updated.2019

human_mouse_con = human_mouse_con[human_mouse_con %notin% c(cell_cycle$s.genes, cell_cycle$g2m.genes)]

human_mouse_con = human_mouse_con[str_detect(human_mouse_con, pattern = "^MT-", negate = T)]

human_mouse_con = human_mouse_con[str_detect(human_mouse_con, pattern = "^RP[LS]", negate = T)]

# human training data preparation -----------------------------------------

# subset human seurat obj

sc_human = sc_human[rownames(sc_human) %in% human_mouse_con, ]

# saveRDS(sc_human, file = 'project/convert_mouse_to_human/human_mye_filtered.rds')

# subset the seurat object to desired celltype

sc_human = subset(sc_human, sampletype == 'scRNA-seq')

sc_human = subset(sc_human, clusters_2 %notin% c('Cycling cells_H', 'Cycling cells_M', 'Cycling cells_T'))

sc_human$clusters_2 = as.character(sc_human$clusters_2)
sc_human$clusters_2[sc_human$clusters_2 == "Neutrophils"] = "Myeloid cells"
sc_human$clusters_2[sc_human$clusters_2 %in% c("T cells", "NK cells")] = "T/NK"

sc_human$clusters_2 = factor(sc_human$clusters_2, levels = unique(sc_human$clusters_2))

Idents(sc_human) = sc_human$clusters_2

# select the 300 top hvg for svm training

sc_human = sc_human %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 300)

var_features = sc_human@assays$RNA@var.features

sc_human_sub = subset(sc_human, downsample = 607)

sc_human_sub = sc_human_sub[rownames(sc_human_sub) %in% var_features, ]

sc_human_train = subset(sc_human_sub, downsample = 607*0.80)

sc_human_test = sc_human_sub[,colnames(sc_human_sub) %notin% colnames(sc_human_train)]

saveRDS(sc_human_train, file = 'data/sc_human_train.rds')
saveRDS(sc_human_test, file = 'data/sc_human_test.rds')
convertFormat(sc_human_train, from = 'seurat', to = 'anndata', outFile = 'data/sc_human_train.h5ad', main_layer = 'data')
convertFormat(sc_human_test, from = 'seurat', to = 'anndata', outFile = 'data/sc_human_test.h5ad', main_layer = 'data')

# harmonize the mouse genes -----------------------------------------------

human_hvg = rownames(sc_human_train)

mouse_hvg = names(human_mouse_con[match(human_hvg, human_mouse_con)])

sc_mouse = NormalizeData(sc_mouse)

sc_mouse_train = sc_mouse[rownames(sc_mouse) %in% mouse_hvg, ]

# ensure the order is the same

mat = sc_mouse_train@assays$RNA@data

mat = mat[match(mouse_hvg, rownames(mat)), ]
head(rownames(mat))

sc_mouse_train = CreateSeuratObject(counts = mat, meta.data = sc_mouse_train@meta.data)

convertFormat(sc_mouse_train, from = 'seurat', to = 'anndata', outFile = 'data/sc_mouse_all.h5ad', main_layer = 'data')
