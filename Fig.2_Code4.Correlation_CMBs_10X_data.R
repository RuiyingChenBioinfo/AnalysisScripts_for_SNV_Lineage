library(Seurat)
library(ggplot2)
library(png)
library(dplyr)
library(reshape2)
library(tibble)
library(ggpubr)

rds <- readRDS("Processed_tomato_10X_CMBinfo.rds")

AvgExp <- as.data.frame(AverageExpression(rds, group.by = "CMB_info", assays = "SCT"))

corr1 <- cor(AvgExp, method = "pearson")

library(corrplot)
pdf("Corr_CMB_Heatmap.pdf", height = 9, width = 10)
corrplot(as.matrix(corr1), title = "CMB Pearson correlation",
         method = 'square', type = "lower",
         order = 'original', 
         addCoef.col = 'black',
         number.cex = 0.8, diag = T, 
         col.lim = c(min(corr1),1),
         #col = rev(COL2('RdBu', 100))[10:90], 
         col = COL2('RdYlBu', 100)[0:88], 
         cl.pos = 'b', 
         tl.cex = 0.8, 
         tl.col = "black")
dev.off()
