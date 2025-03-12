# snRNA tumor --- NMF with refined process: refined normalization method and NMF algorithum

# Load libraries
library(NMF)
library(data.table)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(AnnoProbe)
library(Seurat)
library(parallel)
library(ggdendro)
library(corrplot)
library(ComplexHeatmap)
library(RColorBrewer)
library(pagoda2)
library(paletteer)
nmf.options(maxIter = 12000, verbose = 1)

# read data and metadata

snRNA_list <- SplitObject(snRNA_tumor, split.by = "patient_id")

patient <- names(snRNA_list)

# Here, we need to normalize, scale each sample SAPERATELY, in case to involve bias 
# run nmf
topn <- 10000
rank <- 10
for(i in patient){
  # Create directory
  if (!dir.exists("data/cnmf4")){
    dir.create("data/cnmf4")
  }
  if (!dir.exists(paste("data/cnmf4", i, sep = "/"))){
    dir.create(paste("data/cnmf4", i, sep = "/"))
  }
  
  # nmf matrix processing
  
  ## First, subset the seu obj to only contain protein genes and rm MT, Rbio genes
  snRNA <- snRNA_list[[i]]

  ## MT and Ribo genes
  mito.genes <- grep(pattern = "^MT-", x = rownames(snRNA), value = TRUE)
  rbl.genes <- grep(pattern = "^RB-", x = rownames(snRNA), value = TRUE)
  rsl.genes <- grep(pattern = "^RS-", x = rownames(snRNA), value = TRUE)
  rpl.genes <- grep(pattern = "^RPL-", x = rownames(snRNA), value = TRUE)
  rbl.genes <- grep(pattern = "^RBL-", x = rownames(snRNA), value = TRUE)
  rps.genes <- grep(pattern = "^RPS-", x = rownames(snRNA), value = TRUE)
  rbs.genes <- grep(pattern = "^RBS-", x = rownames(snRNA), value = TRUE)
  rbl1.genes <- grep(pattern = "^RB", x = rownames(snRNA), value = TRUE)
  rsl1.genes <- grep(pattern = "^RS", x = rownames(snRNA), value = TRUE)
  rpl1.genes <- grep(pattern = "^RPL", x = rownames(snRNA), value = TRUE)
  rbl1.genes <- grep(pattern = "^RBL", x = rownames(snRNA), value = TRUE)
  rps1.genes <- grep(pattern = "^RPS", x = rownames(snRNA), value = TRUE)
  rbs1.genes <- grep(pattern = "^RBS", x = rownames(snRNA), value = TRUE)
  
  counts <- snRNA@assays$RNA@counts
  counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
  counts <- counts[rowSums(counts) > 0,] # keep only detected genes
  snRNA <- subset(snRNA, features = rownames(counts))
  
  ## Second, doing SCTransform normalization
  ## We won't need to contain nCount_RNA or mt_ratio, because it has been removed from our matrix
  snRNA <- SCTransform(snRNA, vst.flavor = "v2", method = "glmGamPoi", verbose = F, return.only.var.genes = F)
  
  ## Third, we will process the nmf matrix
  nmf_mat <- GetAssayData(snRNA, slot = "scale.data")
  nmf_mat[nmf_mat < 0] <- 0 # non-negative
  nmf_mat <- nmf_mat[rowSums(nmf_mat) > 0,] # remove zero-sum rows
  gene_sd <- apply(nmf_mat, 1, sd)
  top_gene <- gene_sd[order(gene_sd, decreasing = T)][1:topn]
  nmf_mat <- nmf_mat[rownames(nmf_mat) %in% names(top_gene),] # select top genes
  ina <- which(colSums(is.na(nmf_mat)) == 0)
  nmf_mat <- nmf_mat[,ina] # remove column with NA values
  
  # run nmf
  ## We will use "nndsvd" method to perform concensus NMF
  res <- nmf(nmf_mat, rank = rank, method = "brunet", seed = "nndsvd", .options = "p20v1")
  
  ## better, we will next try to compare multiple algorithms in nmf, and choose the best one
  
  # extract signature
  signature <- NMF::basis(res)
  colnames(signature) <- paste(i, 1:rank, sep = "_")
  signature <- as.data.frame(signature)
  
  # save data
  write.table(signature, paste0("data/cnmf4/", i, "/signature_", topn, ".txt"), sep = "\t")
  saveRDS(res, file = paste0("data/cnmf4/", i, "/result_", rank, ".rds"))
  print(paste0("NMF for ", i, " is done!"))
}

#################################################################################
# Extract NMF program

topn <- 10000
ranks <- 10
topRank <- 100
programG <- list()
patient <- c('PT1', 'PT2', 'PT3', 'PT4', 'PT5', 'PT6', 'PT7', 'PT8', 'PT9', 'PT10', 'PT11', 'PT12')
## extract program
for (i in seq_along(patient)){
  filedir <- paste0("./nmf_res/", patient[i], "/signature_", topn, ".txt")
  geneloading <- read.table(filedir, header = T, sep = "\t")
  geneloading$maxC <- apply(geneloading, 1, which.max) %>% 
    paste0(patient[i], "_", .)
  
  topgenelist <- rownames_to_column(geneloading, var = "gene") %>%
    pivot_longer(., cols = starts_with(c("P", "R")), 
                 names_to = "program", values_to = "loading")
  
  topgenelist <- dplyr::filter(topgenelist, maxC == program) %>% 
    group_by(maxC) %>% top_n(n = topRank, wt = loading)
  topgenelist <- split(topgenelist$gene, topgenelist$maxC)
  programG <- c(programG, topgenelist)
}

snRNA <- snRNA_tumor

## MT and Ribo genes
mito.genes <- grep(pattern = "^MT-", x = rownames(snRNA), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(snRNA), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(snRNA), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(snRNA), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(snRNA), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(snRNA), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(snRNA), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(snRNA), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(snRNA), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(snRNA), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(snRNA), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(snRNA), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(snRNA), value = TRUE)

counts <- snRNA@assays$RNA@counts
counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
counts <- counts[rowSums(counts) > 0,] # keep only detected genes
snRNA <- subset(snRNA, features = rownames(counts))

# Calculate the program score, using pogoda2
snRNA <- SCTransform(snRNA, vst.flavor = "v2", method = "glmGamPoi", verbose = F, return.only.var.genes = F)

exp_mat <- GetAssayData(snRNA, slot = "data")
exp_mat <- as.matrix(exp_mat)

score_list <- list()
score_list <- lapply(programG, function(x){
  score <- score.cells.puram(data = t(exp_mat), signature = x)
  return(score)
})

score_mat <- sapply(score_list, FUN = function(x) x, simplify = T)

saveRDS(score_mat, file = 'score_mat.rds')

###############################################################################

# filter signature with low sds -- optional

## calculate sample-specific sds
sds <- c()
sds <- apply(score_mat, 2, FUN = "sd")
sds <- as.data.frame(sds)
ggplot(sds) + 
  geom_histogram(aes(x = sds))

## select the signature that have sd > 0.1
sds$group <- str_split(rownames(sds), pattern = "_") %>% 
  sapply("[", 1)
sds$program <- rownames(sds)
sds$group <- strsplit(sds$program, split = "_", fixed = T) %>% 
  sapply(FUN = "[", 1)

sds_f <- sds %>% 
  group_by(group) %>% 
  slice_max(order_by = sds, n = 7)

## subset the matrix 
score_f <- score_mat[,colnames(score_mat) %in% sds_f$program]

#######################################################################################

# calculate the correlation and visualization

library(corrplot)
pear_cor <- cor(x = score_mat, method = "pearson")

corrplot(corr = pear_cor, 
         method = "color", 
         order = "hclust", 
         hclust.method = "ward.D2", 
         addrect = 9, 
         tl.pos = "n", 
         col = colorRampPalette(rev(brewer.pal(name = 'RdBu', 11)))(50)) ## remove one signature with no correlation

cororder <- corrMatOrder(pear_cor, order = "hclust", hclust.method = "ward.D2")
pear_cor_hc <- pear_cor[cororder, cororder]
pear_dist <- as.dist(1 - pear_cor_hc)

# Extract program

tree <- hclust(as.dist(1 - pear_cor_hc), method = "ward.D2")
clus <- cutree(tree, 9)
table(clus)

# 提取signature
ProSig <- split(names(clus), clus) 
names(ProSig) <- paste0("tumorsig", names(ProSig))
ProSig <- lapply(ProSig, function(z){
  programG[which(names(programG) %in% z)] %>% unlist() %>% as.character() %>% 
    unique()
})
sapply(ProSig, length)

metalist <- split(names(clus), clus) 
patientLoading <- lapply(patient, function(z){
  filedir <- paste0("nmf_res/", z, "/signature_", topn, ".txt")
  geneloading <- read.table(filedir, header = T, sep = "\t")
  data.frame(Gene = rownames(geneloading), geneloading)
})
AllLoading <- Reduce(function(x, y)merge(x = x, y = y, by = "Gene", all = T), patientLoading)
head(AllLoading)


metaGene <- list()
for (i in 1:length(metalist)){
  programs <- metalist[[i]]
  Selgene <- ProSig[[i]]
  metaGene[[i]] <- AllLoading[AllLoading$Gene %in% Selgene, 
                              colnames(AllLoading) %in% c("Gene",  programs)] %>%
    pivot_longer(cols = starts_with(c("P", "R")), 
                 names_to = "program", values_to = "loading")  %>% na.omit() %>%
    dplyr::group_by(Gene) %>% summarise(Avgloading = mean(loading)) %>% 
    top_n(n = 30, wt = Avgloading) %>% pull(Gene)
}
metaGene_top50[[3]]

saveRDS(metaGene, file = 'metaGene_top50.rds')
saveRDS(metaGene, file = "metaGene_top30.rds")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

enrich_list <- list()

for(i in seq_along(metaGene)){
  enrichres <- enrichGO(gene = metaGene[[i]], OrgDb = org.Hs.eg.db, 
                        keyType = "SYMBOL", ont = 'ALL', pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  enrichres <- as.data.frame(enrichres)
  write_csv(enrichres, file = paste0('nmf_res/enrich_sig', i, '.csv'))
}

#####################################################################################
# visualization
pear_cor_hc[pear_cor_hc < 0] <- 0
pear_cor_hc[pear_cor_hc > 0.9] <- 0.9

patients = rownames(pear_cor_hc) %>% 
  strsplit(split = "_", fixed = T) %>% 
  sapply(FUN = "[", 1)

col_anno <- data.frame(row.names = rownames(pear_cor_hc), patients = patients)
cluster = clus[rownames(col_anno)]
col_anno$cluster = as.character(cluster)
anno_color = as.character(paletteer_d("colorBlindness::paletteMartin"))[2:13]
names(anno_color) = patient
pheatmap(mat = pear_cor_hc, 
         color = colorRampPalette(c('white', '#F9F6B7', '#F4E228', '#E61A13', '#330204'))(100),
         clustering_distance_rows = pear_dist, 
         clustering_distance_cols = pear_dist, 
         clustering_method = "ward.D2",
         show_colnames = F, 
         show_rownames = F, 
         annotation_col = col_anno, 
         annotation_colors = list(patients = anno_color))
)

pheatmap(mat = pear_cor_hc, 
         color = colorRampPalette(c('white', '#FCFFC9', '#F4E228', '#E61A13', '#5A1833', '#1D0B14'))(50),
         clustering_distance_rows = pear_dist, 
         clustering_distance_cols = pear_dist, 
         clustering_method = "ward.D2",
         show_colnames = F, 
         show_rownames = F, 
         annotation_col = col_anno
)

"%notin%" = Negate("%in%")
col_anno = col_anno %>% 
  dplyr::count(cluster, patients) %>% 
  filter(cluster %notin% c("9"))
ggplot(col_anno) + 
  geom_bar(aes(x = cluster, y = n, fill = patients), stat = "identity", position = "fill") + 
  scale_fill_manual(values = as.character(paletteer_d("colorBlindness::paletteMartin"))[2:13]) + 
  theme_bw()
########################################################################################

