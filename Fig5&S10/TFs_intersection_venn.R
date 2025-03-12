# TFs intersection

library(ggvenn)
library(tidyverse)

# read TFs
sig_tfs_mouse <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/mouse_pyscenic/sig_tfs_mouse_trem2.rds")
sig_tfs_human <- readRDS("~/R/Projects/snRNA_scRNA_hcc/project/Myeloid_cell/pyscenic/sig_tfs_trem2_human.rds")
# convert mouse to human

mouse_human_genes = read.csv("~/R/Projects/snRNA_scRNA_hcc/project/svm/HOM_MouseHumanSequence.rpt.txt",sep="\t")

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

sig_tfs_mouse_human = convert_mouse_to_human(gene_list = sig_tfs_mouse)

# intersection
library(ggvenn)
library(RColorBrewer)

gene_list = list(`Human TREM2+ TAM` = sig_tfs_human, `Mouse Trem2+ TAM` = sig_tfs_mouse_human)
myCol <- brewer.pal(3, "Pastel1")

ggvenn(gene_list, fill_color = myCol, stroke_color = NA)

a = base::intersect(sig_tfs_human, sig_tfs_mouse_human)
a
