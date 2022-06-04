rm(list=ls())
getwd()

library(maptools)
library(ggplot2)
library(plyr)
china_map <- readShapePoly("bou2_4p.shp")
china_map1 <- fortify(china_map)
china_data <- read.table("clipboard",header = T, sep = "\t")

ggplot()+
  geom_polygon(data=china_map1, aes(x=long, y=lat, group=group), fill="grey95", colour="grey60")+ 
  geom_point(data=china_data, aes(x = jd,y = wd, size=zb, fill=zb2), shape=21, colour="black")+ 
  scale_size_area(max_size=10)+         
  scale_fill_gradient2(low="DarkCyan", high="red", 
                       midpoint=median(na.omit(china_data$zb2)))+   
  coord_map("polyconic") +ggtitle("Heat&Bubble plot")+
  #geom_text(aes(x=jd+5,y=wd,label=city),size =3,data=china_data)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    #legend.position = "none"
  )  


