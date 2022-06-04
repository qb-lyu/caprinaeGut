library(UpSetR)
library(ggplot2)
library(ggpubr)
# pie

dt.pie <- read.table("profile.pie.txt",header = T, sep = "\t")

pie_hg = function(dt){
  dt$X <- factor(dt$X,levels = dt$X)
  ggplot(dt, aes(x = 2,y = HG, fill = X))+
    geom_bar(stat = "identity",color = "white")+
    
    coord_polar(theta = "y",start = 0)+
    labs(title="HG")+
    theme_void()
}
pie_nhp = function(dt){
  dt$X <- factor(dt$X,levels = dt$X)
  ggplot(dt, aes(x = 2,y = NHP, fill = X))+
    geom_bar(stat = "identity",color = "white")+
    coord_polar(theta = "y",start = 0)+
    labs(title="NHP")+
    theme_void()
}
pie_rug = function(dt){
  dt$X <- factor(dt$X,levels = dt$X)
  ggplot(dt, aes(x = 2,y = RUG, fill = X))+
    geom_bar(stat = "identity",color = "white")+
    coord_polar(theta = "y",start = 0)+
    labs(title="RUG")+
    theme_void()
}
pie_sheep = function(dt){
  dt$X <- factor(dt$X,levels = dt$X)
  ggplot(dt, aes(x = 2,y = Sheep, fill = X))+
    geom_bar(stat = "identity",color = "white")+
    coord_polar(theta = "y",start = 0)+
    labs(title="Sheep")+
    theme_void()
}


hg <- pie_hg(dt.pie)
nhp <- pie_nhp(dt.pie)
rug <- pie_rug(dt.pie)
sheep <- pie_sheep(dt.pie)

ggarrange(hg,nhp,rug,sheep, ncol = 2, nrow = 2)

### UpSet
dt <- read.table("profile.upset.txt",sep = "\t",header = T,check.names = F) 

upset(dt,sets=c("Sheep_MAG","NHP","RUG","Uhgg","RGIG"),keep.order = TRUE,
      order.by = "degree",number.angles = 30,
      point.size = 5,line.size = 1.5, text.scale = c(2, 2, 2, 2, 2, 2))
