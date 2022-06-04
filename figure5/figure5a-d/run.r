library(ggplot2) 
library(ggpubr)
library(reshape2) 
library(vegan) 
library(ggrepel) 
library(aplot) 
library(gridExtra)

########导入数据#####
raw <- read.table("profile.reads.norm",sep="\t", header=T,row.names = 1)
group <- read.table("map.txt",sep="\t",header=T, row.names = 1)
raw_t <- t(raw)
raw_t <- raw_t[rownames(group),]

########PCoA#####

dist <- raw_t %>% vegdist(method="bray")
pcoa <- dist %>% cmdscale(eig=TRUE, k = 10)
temp <- pcoa$points %>% as.data.frame()
pc_importance <- round(pcoa$eig/sum(pcoa$eig)*100, digits = 2)

plot_data <- merge(temp, group, by.x = "row.names",by.y = "row.names")

Pcoa = function(plot_data, group_name, pc1, pc2){
  First_pc <- paste("V",pc1,sep="")
  Second_pc <- paste("V",pc2,sep="")
  ggscatter(plot_data, x = First_pc, y = Second_pc, color = group_name, ellipse = T, 
            ellipse.type = 'confidence', mean.point = T, star.plot = T)+
    labs(x=paste("PCoA ",pc1," (", pc_importance[pc1],digits=4,"%)", sep=""), 
         y=paste("PCoA ",pc2," (", pc_importance[pc2],digits=4, "%)", sep=""))+
    theme_bw()#+theme(legend.position="none") ## 取消图例
}

p1 = 1
p2 = 2
Group_title <- c("group")

mycols<-c("#8dd3c7","#fedf61","#bebada","#fb8072","#80b1d3","#fdb462")

pcoa_plot <- Pcoa(plot_data, Group_title, p1, p2)+scale_color_manual(values = mycols)
pcoa_plot

adonis_sample <-adonis(raw_t~group[,1], na.rm = T) # R2 = 0.33203    
adonis_sample

########解释量#####
sig_list <- t(combn(unique(plot_data$group),2)) #显著性两两比较
compire<-list()
for(i in 1:nrow(sig_list)){
  compire[[i]]<-as.character(sig_list[i,])
}

compire <- list(c('AH','XZ'),c('AH','GX'),c('GX','JL'),c('GX','SD'),c('GX','SX')
                ,c('JL','SD'),c('JL','XZ'),c('SD','XZ'),c('SX','XZ'))

compire2 <- list(c('AH','XZ'),c('GX','XZ'))

plot_data$Row.names <- factor(plot_data$Row.names,levels = plot_data$Row.names) 

lv_box = function(dt,pc,compire){
  ggplot(dt, aes(x=group,y=pc))+     
    geom_boxplot(aes(fill=group))+
    scale_fill_manual(values=c("#8dd3c7","#fedf61","#bebada","#fb8072","#80b1d3","#fdb462"))+
    geom_signif(comparisons = compire,step_increase=0.1,map_signif_level = T,test = wilcox.test)+
    theme_bw()
}

p1 <- lv_box(plot_data,plot_data$V1,compire)+theme(legend.position="none")+coord_flip()
p2 <- lv_box(plot_data,plot_data$V2,compire2)+theme(legend.position="none")

pca.var = cbind(pc_importance[c(1:10)]) %>% as.data.frame()
colnames(pca.var)[1] <- 'var'
pca.var$pc <- c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
pca.var$pc <- factor(pca.var$pc, levels = pca.var$pc)

pcoa.10 <- ggplot(pca.var, aes(pc,var)) +
  geom_bar(stat = 'identity')+
  ylim(0, 16)+
  geom_text(aes(label = paste(var,' %', sep = "")), nudge_y = 0.5)+
  labs(x = '主坐标',y = '解释量（%）')+
  theme_bw()+
  theme(axis.text.x = element_text(angle=70,vjust=.5))
pcoa.10
plot1 <- grid.arrange(arrangeGrob(p2,pcoa_plot, widths=c(1,2)),
                      arrangeGrob(pcoa.10,p1, widths=c(1,2)),heights=c(2,1))

########多样性#####
sig_list <- t(combn(unique(group$group),2)) #显著性两两比较
compire<-list()
for(i in 1:nrow(sig_list)){
  compire[[i]]<-as.character(sig_list[i,])
}

Diversity <- data.frame(matrix(ncol= 7, nrow=30))
colnames(Diversity) <- c('ID','Observed','Chao1','ACE','Shannon','Simpson','Invsimpson')

Diversity$ID <- rownames(raw_t)
Diversity$Observed <- estimateR(ceiling(raw_t))[1, ]
Diversity$Chao1 <- estimateR(ceiling(raw_t))[2, ]
Diversity$ACE <- estimateR(ceiling(raw_t))[4, ]
Diversity$Shannon <- diversity(raw_t, index = "shannon")
Diversity$Simpson <- diversity(raw_t, index = "simpson")
Diversity$Invsimpson <- diversity(raw_t, index = "invsimpson")

Diversity.plot <- merge(Diversity, group, by.x = "ID", by.y = "row.names")

lv_box2 = function(dt,pc,compire){
  ggplot(dt, aes(x=group,y=pc))+     
    geom_boxplot(aes(fill=group))+
    scale_fill_manual(values=c("#8dd3c7","#fedf61","#bebada","#fb8072","#80b1d3","#fdb462"))+
    geom_signif(comparisons = compire,step_increase=0.1,map_signif_level = T,test = wilcox.test)+
    theme_bw()
}

Observed <- lv_box2(Diversity.plot, Diversity.plot$Observed,compire)+labs(title="Observed")
Chao1 <- lv_box2(Diversity.plot, Diversity.plot$Chao1,compire)+labs(title="Chao1")
ACE <- lv_box2(Diversity.plot, Diversity.plot$ACE,compire)+labs(title="ACE")
Shannon <- lv_box2(Diversity.plot, Diversity.plot$Shannon,compire)+labs(title="Shannon")
Simpson <- lv_box2(Diversity.plot, Diversity.plot$Simpson,compire)+labs(title="Simpson")
Invsimpson <- lv_box2(Diversity.plot, Diversity.plot$Invsimpson,compire)+labs(title="Invsimpson")

# write.table(Diversity.plot,"Diversity.plot",sep = "\t",quote=F)

grid.arrange(Shannon, Simpson, Invsimpson, 
             Observed, Chao1,ACE,
             nrow = 2, ncol = 3)

result = rbind()
for(i in runif(1000,1,10000)){
  set.seed(i)
  xz.div <- subset(Diversity.plot, group=="XZ")[,c(5)]
  plain.div <- subset(Diversity.plot, group!="XZ")[,c(5)]
  random.div <- sample(plain.div, 5)
  p <- wilcox.test(xz.div,random.div)$p.value
  temp_result = data.frame(pvalue=p)
  result = rbind(result, temp_result)
}


########物种组成#####
library(ggsci)
group2 <- read.table("map.txt",sep="\t",header=T)
top10 <- read.table("profile.family10.txt", sep="\t", header=T,row.names = 1)
top10.d <- vegdist(t(top10))
top10.hc <- hclust(top10.d)
lable_order <- top10.hc$labels[top10.hc$order]

top10$family <- rownames(top10)
top10.m <- melt(top10,ID = "family")
top10.m <- merge(top10.m,group2,by.x = "variable",by.y = "Sample")

top10.m$family <- factor(top10.m$family, levels = unique(top10.m$family))
top10.m$variable <- factor(top10.m$variable, levels = unique(lable_order))

ggplot(top10.m, aes(variable, y = value, fill = family)) +
  geom_bar(position = 'stack', width = 1,stat="identity") +
  labs(x = '', y = 'Relative Abundance(%)') + 
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_wrap(group~.,scales = 'free_x')


########UpSet######
library(UpSetR)
dt <- read.table("upset.txt",sep = "\t",header = T,check.names = F) 
metadata <- data.frame(Sets = c("AH","GX","JL","SD","SX","XZ"),
                       Cities = c("AH","GX","JL","SD","SX","XZ"))

upset(dt, point.size = 3,keep.order = TRUE,order.by = "degree",nsets = 6,nintersects = 70,
      line.size = 1,
      set.metadata = list(
        data = metadata, 
        plots = list(
          list(
            type = "matrix_rows", 
            column = "Cities", 
            colors = c(
              AH = "#8dd3c7", 
              GX = "#fedf61", 
              SD = "#fb8072",
              SX = "#80b1d3",
              JL = "#bebada",
              XZ = "#fdb462"),
            alpha = 0.5)
        )
      ),
      queries = list(list(query = elements,
                          params = list("phylum","Firmicutes"),
                          color = "#fb9a99",active = F),
                     list(query = elements,
                          params = list("phylum","Bacteroidetes"),
                          color = "#1f78b4",active = F),
                     list(query = elements,
                          params = list("phylum","Proteobacteria"),
                          color = "#ff7f00",active = F),
                     list(query = elements,
                          params = list("phylum","Actinobacteria"),
                          color = "#a6cee3",active = F),
                     list(query = elements,
                          params = list("phylum","Verrucomicrobia"),
                          color = "#6a3d9a",active = F)
      )
)




