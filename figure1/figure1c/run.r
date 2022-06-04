library(ggplot2) 
library(vegan) #vegan 2.5-7

raw <- read.table("profile.reads.norm",sep="\t", header=T,row.names = 1)
raw_t <- t(raw)

get_plot_data=function(dd, nm=nm){
  dd.curve=specaccum(dd, method = "random")
  dd.curve.data=data.frame(Sites=dd.curve$sites, Richness=dd.curve$richness, SD=dd.curve$sd)
  dd.curve.data$label=rep(nm, nrow(dd.curve.data))
  dd.curve.data
}

total.data = get_plot_data(raw_t, nm="Sample")

rara.plot <- ggplot(total.data, aes(x=Sites, y=Richness, color=label))+
  geom_line() + 
  geom_errorbar(aes(ymax = Richness + SD, ymin = Richness - SD), width = .5)+
  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(name="Number of samples")+
  scale_y_continuous(name="") #
rara.plot
