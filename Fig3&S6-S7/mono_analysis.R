library(tidyverse)
library(Seurat)
library(data.table)
library(harmony)

seu_mono = subset(seu_anno, celltype %in% c("Myeloid cells"))

table(seu_mono$orig.ident)

seu_mono = NormalizeData(seu_mono) %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = "mt_ratio") %>% 
  RunPCA() %>% 
  RunUMAP(reduction = "pca", dims = 1:20)

seu_mono = RunHarmony(seu_mono, group.by.vars = "orig.ident", max_iter = 100, plot_convergence = T)

seu_mono = RunUMAP(seu_mono, reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = .8)

DimPlot(seu_mono, group.by = "orig.ident")
DimPlot(seu_mono)

FeaturePlot(seu_mono, features = c("TREM2", "SPP1"))

VlnPlot(seu_mono, features = c("SPP1", "TREM2"))

DimPlot(seu_mono, group.by = "response", split.by = "response")

VlnPlot(seu_mono, features = "TREM2", group.by = "response")

DimPlot(seu_mono)

markers = FindAllMarkers(seu_mono, only.pos = T)

seu_mono_anno = RenameIdents(seu_mono, "0" = "TREM2+ SPP1+ macrophage", 
                             "1" = "FOLR2+ macrophage", 
                             "2" = "CD1c+ DC", 
                             "3" = "Monocytes", 
                             "4" = "RTM-like TAM", 
                             "5" = "Hepatocytes", 
                             "6" = "NK cells", 
                             "7" = "CELC9A+ DC", 
                             "8" = "Proliferating macrophage", 
                             "9" = "LAMP3+ DC", 
                             "10" = "EC", 
                             "11" = "pDC")

DimPlot(seu_mono_anno)

seu_mono_anno$celltype = seu_mono_anno@active.ident

seu_mono_anno = subset(seu_mono_anno, celltype != "Hepatocytes" & celltype != "NK cells" & celltype != "EC")

seu_mono_anno = RunUMAP(seu_mono_anno, reduction = "harmony", dims = 1:30)

DimPlot(seu_mono, split.by = "condition")

### calculate OR value

seu_mono = subset(seu_mono, celltype != "pDC")
cellInfo.tb = seu_pub_mye@meta.data
cellInfo.tb = cellInfo.tb %>% 
  dplyr::rename("meta.cluster" = DefineTypes, "loc" = condition) %>% 
  as.data.table()

test.dist.table <- function(count.dist,min.rowSum=0)
{
  library(plyr)
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  #count.dist.melt.ext.tb[adj.p.value < 0.05,]
  #count.dist.melt.ext.tb[p.value < 0.05 & OR > 0,]
  
  return(count.dist.melt.ext.tb)
  
}

do.tissueDist <- function(cellInfo.tb,out.prefix,pdf.width=3,pdf.height=5,verbose=0)
{
  library("Startrac")
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb[,meta.cluster:=as.character(meta.cluster)]
  loc.avai.vec <- unique(cellInfo.tb[["loc"]])
  # loc.avai.vec <- intersect(c("P","N","T"),loc.avai.vec)
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  startrac.dist <- unclass(calTissueDist(cellInfo.tb,colname.cluster="meta.cluster"))
  startrac.dist <- startrac.dist[,loc.avai.vec]
  
  cuts <- c(0, 0.8, 1.2,Inf)
  startrac.dist.bin.values <- factor(c("-", "+/-", "+"),levels=c("-", "+/-", "+"))
  startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                              ncol=ncol(startrac.dist))
  colnames(startrac.dist.bin) <- colnames(startrac.dist)
  rownames(startrac.dist.bin) <- rownames(startrac.dist)
  
  #	sscClust:::plot.matrix.simple(freq.dist.bin,
  #								  col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),
  #												   names=0:9 )),
  #								  par.legend=list(labels=rev(sprintf("%s%%~%s%%",10*(0:9),c(10*(1:9),100) )),
  #												  at=0:9),
  #								  out.prefix=sprintf("%s.freq.dist",out.prefix),
  #								  show.number=F,clust.row=T,exp.name=expression(italic(Freq)),
  #								  #palatte=(brewer.pal(n = 7,name = "Blues")),
  #								  #palatte=viridis::viridis(7),
  #								  pdf.width = 4.5, pdf.height = pdf.height)
  
  #	sscClust:::plot.matrix.simple(startrac.dist.bin,
  #								  col.ht=rev(structure(viridis::viridis(3),
  #													   names=levels(startrac.dist.bin.values))),
  #								  out.prefix=sprintf("%s.startrac.dist.bin",out.prefix),
  #								  show.number=F,clust.row=T,exp.name=expression(italic(R)[o/e]),
  #								  pdf.width = pdf.width, pdf.height = pdf.height)
  
  sscVis::plotMatrix.simple(startrac.dist,
                            out.prefix=sprintf("%s.startrac.dist",out.prefix),
                            show.number=F,
                            clust.row=T,
                            #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(R)[o/e]),
                            z.hi=2,
                            #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            #clust.row=T,
                            #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
  
}

OR.dist.mtx = do.tissueDist(cellInfo.tb = cellInfo.tb, out.prefix = "./Fig3")

OR.dist.mtx[OR.dist.mtx > 3] = 3

OR.dist.mtx = as.data.frame(OR.dist.mtx) %>% 
  rownames_to_column(var = "celltype") %>% 
  arrange(NR, desc(NR)) %>% 
  column_to_rownames(var = "celltype")

OR.dist.mtx = OR.dist.mtx[, c(2, 1)]


pheatmap::pheatmap(OR.dist.mtx, cluster_rows = F, cluster_cols = F, color = viridis::viridis(n = 50), border_color = NA)

# use a category visualization

OR.dist.mtx = mutate(OR.dist.mtx, R_c = case_when(
  R < 1 ~ 0, # +-
  R > 1 & R <= 1.5 ~ 1, # +
  R > 1.5 & R <= 2 ~ 2, # ++
  R > 2 ~ 3 # +++
), NR_c = case_when(
  NR < 1 ~ 0, 
  NR > 1 & NR <= 1.5 ~ 1, 
  NR > 1.5 & NR <= 2 ~ 2, 
  NR > 2 ~ 3
))

OR.dist.mtx_c = OR.dist.mtx[, c(3, 4)]

pheatmap::pheatmap(OR.dist.mtx_c, cluster_rows = F, cluster_cols = F, color = colorRampPalette(brewer.pal(n = 8, name = "Reds"))(60), border_color = NA)
FeaturePlot(seu_mono, features = "C1QA")
DimPlot(seu_mono, split.by = "condition") + scale_color_manual(values = c('#006DDBFF', '#F5CDCD', '#E2A7CC', '#378C4F', '#7464AA',  '#B6DBFFFF', '#6DB6FFFF', '#924900FF'))


seu_macro = subset(seu_mono, celltype %in% c("TREM2+ SPP1+ macrophage", "FOLR2+ macrophage", "RTM-like TAM"))

VlnPlot(seu_macro, features = c("TREM2", "SPP1", "GPNMB", "CTSC", "CD9", "MARCO", "FOLR2", "CD163", "CD5L"), group.by = "condition", cols = c("#E41A1C", "#377EB8"))
