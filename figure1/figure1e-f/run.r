library(ggplot2)
library(forcats)
library(ggpubr)

# Family 
dt <- read.table("profile.ANI.txt",header = T,sep = "\t")
dt.f <- dt[-c(grep("other",dt$Family)),]
dt.f$Family <- fct_inorder(dt.f$Family)

plot.f <- ggplot(dt.f,aes(x=Family,y=size))+
            geom_boxplot(aes(fill=Phylum))+
            coord_flip()+
            labs(title="",x="", y = "",fill = "")+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.f

# size GC geom 

dt.s <- read.table('profile.size.txt',header = T,sep = '\t',check.names = T)

ggscatter(dt.s, x="Estimated.genomic.size", y="GC.content", color = "Phylum", shape = 19,
          ellipse = T, ellipse.type = 'confidence',size = 2,
          ellipse.level = 0.95,ellipse.alpha = 0.4,
          mean.point = T, star.plot = F,alpha=0.5,xlim = c(0,7))+
      theme_bw()

