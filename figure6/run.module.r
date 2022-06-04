
pacman::p_load(vegan,dplyr,tibble,tidyr,ggplot2,ggpubr,ggrepel,randomForest,pROC,stringr)
options(stringsAsFactors = FALSE)

K <- read.delim("ko_abund.txt", sep = "\t", header = TRUE, row.names = NULL)
group <- read.delim("group.txt", sep = "\t", header = TRUE, row.names = NULL)
load("kegg_module_database-2021.Rdata")
K <- column_to_rownames(K, var = "name")


# 对每一个KO计算出现在三个样本以上的wilcox秩和检验
sample_xz <- group$sample[group$group == "XZ"]
sample_fxz <- group$sample[group$group == "FXZ"]
count_xz <- apply(K[,sample_xz], 1, function(x) sum(x > 0))
count_fxz <- apply(K[,sample_fxz], 1, function(x) sum(x > 0))


K_t <- data.frame(t(K))
p_value <- c()
for (i in seq_len(ncol(K_t))) {
  data_i <- data.frame(sample = rownames(K_t), abund = K_t[,i])
  data_i$group <- group$group[match(data_i$sample, group$sample)]
  vec1 <- data_i$abund[data_i$group%in%"XZ"]
  vec2 <- data_i$abund[data_i$group%in%"FXZ"]
  wilcox <- wilcox.test(vec1, vec2, paired = FALSE)
  p_value <- c(p_value, wilcox$p.value)
}


# 根据样本>3过滤
p_value_filtered <- data.frame(KO = colnames(K_t), p_value = p_value, count_xz = count_xz, count_fxz = count_fxz) %>%
                           subset(count_xz >= 3 & count_fxz >= 3)
p_value_filtered$p_adj <- p.adjust(p_value_filtered$p_value, method = "BH")
p_value_filtered <- p_value_filtered %>% select(KO, p_value, p_adj)


# 差异倍数foldchange
fc <- p_value_filtered
fc$xz_mean <- rowMeans(K[match(fc$KO,rownames(K)), sample_xz])
fc$fxz_mean <- rowMeans(K[match(fc$KO,rownames(K)), sample_fxz])

foldchange <- c()
direction <- c()
for (i in 1:nrow(fc)) {
  if (fc$xz_mean[i] > fc$fxz_mean[i]) {
  foldchange[i] <- fc$xz_mean[i] / fc$fxz_mean[i]
  direction[i] <- 1
  } else {
  foldchange[i] <- fc$fxz_mean[i] / fc$xz_mean[i]
  direction[i] <- -1
  }
}
fc <- cbind(fc, foldchange, direction)


# 筛选
sum(fc$p_adj < 0.05 & fc$foldchange >= 4)
fc_filtered <- fc %>% subset(fc$p_adj < 0.05 & fc$foldchange >= 4)
fc_filtered$zs <- -qnorm(fc_filtered$p_adj/2)
fc_filtered$rzs <- ifelse(fc_filtered$direction == 1, fc_filtered$zs, -1*fc_filtered$zs)
write.csv(fc_filtered, "wilcox_KO_rzs.csv", quote = F, row.names = F)


# 获得注释到的K参与的所有模块
module <- na.omit(unique(kegg_module_database$Module[match(rownames(K), kegg_module_database$lvE)]))

module_rzs <- rbind()
for (i in 1:length(module)) {
  module_i_KO <- kegg_module_database$lvE[kegg_module_database$Module%in%module[i]]
  module_i_KO <- intersect(module_i_KO, fc_filtered$KO)
  if (!length(module_i_KO) == 0) {
    module_i_KO_rzs <- fc_filtered$rzs[match(module_i_KO, fc_filtered$KO)]
    module_i_rzs <- (1/sqrt(length(module_i_KO)))*sum(module_i_KO_rzs)
    module_rzs <- rbind(module_rzs, data.frame(module = module[i], rzs = module_i_rzs))
    } else {
      next
    }
  }


# 筛选rz大于1.7的item
module_rzs <- module_rzs %>% subset(abs(rzs) >= 1.7)
module_rzs$module_des <- kegg_module_database$Module_des[match(module_rzs$module,kegg_module_database$Module)]
module_rzs$module_des_rectify <- gsub(".{1}$","",unlist(lapply((strsplit(module_rzs$module_des, split = "\\[")), "[", 1)))
module_rzs$lvC <- kegg_module_database$lvC[match(module_rzs$module,kegg_module_database$Module)]
module_rzs$lvB <- kegg_module_database$lvB[match(module_rzs$module,kegg_module_database$Module)]
module_rzs$label <- paste0(module_rzs$module," ", module_rzs$module_des_rectify, "[", round(module_rzs$rzs,digits = 6), "]")
module_rzs$type <- ifelse(module_rzs$rzs > 0, "Tibet", "Plain")
module_rzs <- module_rzs[module_rzs$lvB%in%unique(module_rzs$lvB)[grep("[M|m]etabolism",unique(module_rzs$lvB))],]
module_rzs <- module_rzs %>% select(lvB,lvC,module,module_des_rectify,rzs,type,label)
write.csv(module_rzs, "module_rzs_info.csv", row.names = FALSE)


# 画图
plot_data <- select(module_rzs, module, rzs, label, type)

ggplot(data = plot_data, aes(x = reorder(label, rzs), y = (rzs), fill = type)) +
  geom_bar(stat = "identity", width = 0.9, ) + 
  coord_flip() +
  scale_y_continuous(breaks = seq(-5,8,2), label = as.character(seq(-5,8,2)), limits = c(-5,8)) + 
  scale_x_discrete(position = "top") + 
  theme(panel.grid.major.x = element_line(color = 'gray', size = 0.2, linetype = 8), 
        panel.background = element_rect(fill = "transparent", color = "black"),
        panel.border = element_rect(fill = NA, color = "black", size = 1), 
        legend.key = element_rect(fill = 'transparent'),
        plot.title = element_text(hjust = 0), 
        legend.position = "right",
        legend.title = element_blank())
