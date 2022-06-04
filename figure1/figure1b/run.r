library(ggplot2)
library(viridis)
library(ggridges)
library(forcats)
library(ggpubr)
library(aplot)
library(ggExtra)
library(vegan)
library(ggsci)
library(ggrepel)

# N50 and Depth

dt.d <- read.table('profile.dep.txt',header = T, sep = '\t')

plot_dep <- ggplot(dt.d, aes(x=log10(N50.length), y=log10(Depth), colour=group2))+
  geom_point(size=3,shape = 19,alpha = 0.5)+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_bw()

ggMarginal(plot_dep, type="density",groupColour = T,groupFill = T)

# size GC geom 

dt.s <- read.table('profile.size.txt',header = T,sep = '\t',check.names = T)

ggscatter(dt.s, x="Estimated.genomic.size", y="GC.content", color = "Phylum", shape = 19,
          ellipse = T, ellipse.type = 'confidence',size = 2,
          ellipse.level = 0.95,ellipse.alpha = 0.4,
          mean.point = T, star.plot = F,alpha=0.5,xlim = c(0,7))+
      theme_bw()

