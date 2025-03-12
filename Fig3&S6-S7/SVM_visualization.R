# svm_visualization

library(tidyverse)
library(ggsci)

svm_proba = read_csv('data/cell_prob_landscape.csv')
source('plotting_functions.R')

color = c("#378C4F", "#7BBC5E" , "#E2A7CC", "#D9579B", "#A59ACB", "#7464AA", "#006ddb", "#6db6ff", "#b6dbff", "#f5cdcd")
names(color) = c("B cell", "CAF", "Cycling cells", "Endothelial cells", "Hepatocytes" , "NK cells", "Mono/marco", "Granulocytes", "Plasma cells", "T cells")


# data manipulation -------------------------------------------------------

svm_proba = svm_proba %>% 
  select(-...1)

svm_proba_list = split.data.frame(svm_proba, f = ~ annotation)

svm_list = map2(svm_proba_list, names(svm_proba_list), function(x, y){
  df = x[,y]
  colnames(df) = 'score'
  df$celltype = y
  return(df)
})

svm_data = data.frame()

for(df in svm_list){
  svm_data = rbind(svm_data, df)
}

svm_data = svm_data %>% 
  arrange(desc(score))

svm_data = svm_data %>% 
  mutate(celltype = factor(celltype, levels = c("Plasma cells", "CAF", "Hepatocytes", "Endothelial cells", "Myeloid cells", "T/NK", "B cells")))
plot_boxplot(data = svm_data, x = celltype, y = score, group = celltype, theme = 'classic', compare = F, point = F, jitter_width = .2, point_sz = 1,) + geom_hline(yintercept = 0.4)
