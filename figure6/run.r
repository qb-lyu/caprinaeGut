
pacman::p_load(vegan,dplyr,tibble,tidyr,ggplot2,ggpubr,ggrepel,randomForest,pROC,stringr)

diff <- read.delim("deseq2_diff_sgb.txt", sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = FALSE)
rc <- read.delim("profile.reads.norm", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
rc <- rc[diff$name,]
rc_relative <- apply(rc, 2, function(x) x/sum(x))
rc_relative <- data.frame(rc_relative)


row_annotation <- data.frame(Enrich = factor(diff$Enrich), Phylum = diff$Phylum, Order = diff$Order, row.names = diff$name)
row_annotation$Enrich <- gsub("XZ","Tibet",row_annotation$Enrich)
row_annotation$Enrich <- gsub("plain","Plain",row_annotation$Enrich)


color_phylum <- c(paletteer::paletteer_d("ggsci::default_locuszoom")[-7], "#E7C9C6")
names(color_phylum) <- unique(diff$Phylum)
color_order <- c(paletteer::paletteer_dynamic("cartography::pastel.pal", 20),"#bebebe")
names(color_order) <- unique(diff$Order)


colors <- list(Enrich = c(Tibet = "#10CF9B", Plain = "#A5C249"),
               Phylum = color_phylum,
               Order = color_order)

pheatmap::pheatmap(rc_relative,
                   scale = "row",
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(c("#089392","#EADD97","#CF597E"))(100),
                   border_color = NA,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   cellwidth = 15, 
                   cellheight = 3,
                   treeheight_col = FALSE,
                   treeheight_row = 30,
                   fontsize_col = 5, 
                   fontsize_row = 1,
                   gaps_col = 25,
                   cutree_rows = 2,
                   annotation_row = row_annotation,
                   annotation_colors = colors,
                   # annotation_names_col = T,
                   main = "Heatmap")


