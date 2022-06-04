library(ggplot2)
library(reshape2)
library(vegan)
library(ggpubr)

### annotation ratio
dt <- read.table("ratio.txt",sep = "\t",header = T,check.names = F) 

dt$Phylum <- factor(dt$Phylum,levels = unique(dt$Phylum))

ggplot(dt, aes(x=Phylum,y=r))+     
  geom_boxplot(aes(fill=Phylum))+
  theme_bw()+theme(axis.text.x = element_text(angle=30,vjust=.7))

### pcoa

raw <- read.table("profile.kegg.txt", header=T, sep="\t",row.names = 1)
raw_t <- t(raw) 
group <- read.table("map.txt",sep="\t",header=T) 

raw.dist <- vegdist(raw_t,method="bray") 
raw_pcoa <- cmdscale(raw.dist,eig=TRUE)
temp1 <- raw_pcoa$points[,1:2]
temp1 <- as.data.frame(temp1) 
temp1$ID <- row.names(temp1)
pc_importance <- round(raw_pcoa$eig/sum(raw_pcoa$eig)*100,digits = 2)
data <- merge(temp1,group,by="ID")

raw <- raw[,group$ID]
adonis(t(raw)~group[,2], na.rm=T) # R2 = 0.29406

ggscatter(data, x="V1", y="V2", color = "phylum", 
          ellipse = T, ellipse.type = 'confidence', 
          # shape = "cyl", palette = c("#00AFBB", "#E7B800", "#FC4E07"),                  
          mean.point = T, star.plot = T)

## KEGG lvB

dt <- read.table("profile.lvB.txt",sep = "\t",header = T,check.names = F) 
dt.m <- melt(dt)

ggplot(dt.m, aes(B, value, fill = variable)) +
  geom_bar(position = 'dodge', width = 0.8,stat="identity") +
  labs(x = '', y = 'Average gene copy number') + 
  scale_fill_brewer(palette = "Set3")+
  theme_bw()+coord_flip()

## 
