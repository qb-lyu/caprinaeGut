library(ggplot2) 
library(reshape2)
library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)

dt = as.matrix(read.table("profile.txt",sep = "\t", header = T,row.names = 1,check.names = F))
dt <- t(dt)
map = read.table("P2Fmap.txt",sep = "\t", header = T,check.names = F)
md = read.table("module.txt",sep = "\t", header = T,check.names = F)

kk = heatmap(dt)
row_order = rownames(dt)[kk$rowInd]
col_order = colnames(dt)[kk$colInd]
dt = dt[row_order, col_order]

rownames(map)=map$family
map = map[colnames(dt),] 
csd = data.frame(enriched=map$phylum)
csd$enriched = factor(csd$enriched, unique(csd$enriched))

rownames(md)=md$modules
md  = md[rownames(dt),]
rsd = data.frame(enriched=md$B)
rsd$enriched = factor(rsd$enriched, unique(rsd$enriched))

mycols = colorRamp2(breaks=c(0, 1, 2),
                    colors=c("#ece2f0","#a6bddb", "#1c9099"))


Heatmap(as.matrix(dt),bottom_annotation = ha2,
        cluster_rows=FALSE, cluster_columns=FALSE, # 是否启用聚类
        row_split=rsd, column_split=csd, # 行/列分割
        row_title_side='right',column_title_side = 'bottom', # 文字显示方位
        row_title_rot = 0, column_title_rot=90, # 文字旋转角度
        border=T,
        show_row_names=F, show_column_names = T, 
        row_gap = unit(0, 'mm'), column_gap = unit(0, 'mm'),
        col=mycols
)