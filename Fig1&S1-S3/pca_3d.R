# install.packages("plotly")
library(plotly)
library(Seurat)
library(reticulate)

seu_tumor = subset(seu_tumor, archtype != "unassigned")
Idents(seu_tumor) = seu_tumor$archtype
seu_tumor_sub = subset(seu_tumor, downsample = 2000)

seu_tumor_sub = seu_tumor_sub %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

pca_data = seu_tumor_sub@reductions$pca@cell.embeddings

pca_data = pca_data[, 1:3]

pca_data = as.data.frame(pca_data)

pca_data$group = seu_tumor_sub$archtype

pca_data$patient = seu_tumor_sub$patient_id

fig = plot_ly(pca_data, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~group, 
              marker = list(size = 5))

fig = fig %>% add_markers()

fig = fig %>% layout(scene = list(xaxis = list(showticklabels = T, title = "PC1"), 
                                  yaxis = list(showticklabels = T, title = "PC2"), 
                                  zaxis = list(showticklabels = T, title = "PC3")))

fig = fig %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = 800,
      height = 600
    )) # gain a pdf file rather than png

fig

# 3D plot with color scaling on gene expression

gene = "SAA1" # genes to plot

exp = seu_tumor_sub@assays$RNA@data

exp_gene = exp[gene, ]

pca_data[["exp"]] = exp_gene

fig = plot_ly(pca_data, x = ~PC_1, y = ~PC_2, z = ~PC_3, 
              marker = list(color = ~exp, colorscale = "Reds", showscale = TRUE, size = 5))

fig = fig %>% add_markers()

fig = fig %>% layout(scene = list(xaxis = list(showticklabels = F, title = ""), 
                                  yaxis = list(showticklabels = F, title = ""), 
                                  zaxis = list(showticklabels = F, title = "")))

fig = fig %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = 600,
      height = 700
    )) # gain a svg file rather than png

fig

