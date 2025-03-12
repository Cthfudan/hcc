library(ggvenn)

gene_list = list(signature_gene = rownames(surv_TCGA), 
                 metabolism_gene = metabolic_gene, 
                 up_gene = gene_up)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

ggvenn(gene_list, fill_color = myCol, stroke_color = NA)
