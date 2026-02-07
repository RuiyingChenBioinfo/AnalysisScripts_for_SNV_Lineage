library(sp)
library(SeuratObject)
library(Seurat)
library(glmGamPoi)
library(ggplot2)
library(ggalluvial)
library(tidyverse)
library(stringr)
library(dplyr)
library(forcats)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(ComplexHeatmap)
library(circlize)

ref   <- readRDS("GSE152766_Ground_Tissue_Atlas.rds")
query <- readRDS("AT_Multiome_RootTip.RNA.rds")
cat("ref:", nrow(ref), "genes x", ncol(ref), "cells\n")
cat("query:", nrow(query), "genes x", ncol(query), "cells\n")

label_col <- "time.anno.fine"
cat("\nref time.anno.fine :\n")
print(table(ref[[label_col]]))

target_celltypes <- c("QC", "Endodermis", "Cortex")
cat("\n query :\n")
print(table(query$Cell.type.RNA))

query_subset <- subset(query, subset = Cell.type.RNA %in% target_celltypes)
cat("\nquery :", ncol(query_subset), "\n")
print(table(query_subset$Cell.type.RNA))

DefaultAssay(ref)   <- "RNA"
DefaultAssay(query_subset) <- "RNA"
genes_common <- intersect(rownames(ref[["RNA"]]), rownames(query_subset[["RNA"]]))

ref <- subset(ref, features = genes_common)
query_subset <- subset(query_subset, features = genes_common)
cat("subset  ref:", nrow(ref), "genes\n")
cat("subset  query_subset:", nrow(query_subset), "genes\n")

ref   <- NormalizeData(ref, normalization.method = "LogNormalize", verbose = FALSE)
query_subset <- NormalizeData(query_subset, normalization.method = "LogNormalize", verbose = FALSE)

ref   <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
query_subset <- FindVariableFeatures(query_subset, nfeatures = 2000, verbose = FALSE)
features_use <- intersect(VariableFeatures(ref), VariableFeatures(query_subset))
cat("HVG :", length(features_use), "\n")

ref   <- ScaleData(ref, features = features_use, verbose = FALSE)
query_subset <- ScaleData(query_subset, features = features_use, verbose = FALSE)

npc_ref   <- min(50, length(features_use) - 1, ncol(ref) - 1)
npc_query <- min(50, length(features_use) - 1, ncol(query_subset) - 1)

ref[["pca"]] <- NULL
query_subset[["pca"]] <- NULL

ref   <- RunPCA(ref, features = features_use, npcs = npc_ref, verbose = FALSE)
query_subset <- RunPCA(query_subset, features = features_use, npcs = npc_query, verbose = FALSE)

cat("ref PCA:", dim(Embeddings(ref, "pca")), "\n")
cat("query_subset PCA:", dim(Embeddings(query_subset, "pca")), "\n")

max_dims <- min(30, npc_ref, npc_query)
use_dims <- 1:max_dims

cat("\n FindTransferAnchors.. .\n")
anchors <- FindTransferAnchors(
  reference = ref,
  query = query_subset,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  features = features_use,
  reduction = "pcaproject",
  verbose = TRUE
)

n_anchors <- nrow(anchors@anchors)
cat("\n anchors :", n_anchors, "\n")

k_wt <- max(1, min(10, n_anchors - 1))
ref_labels <- as.character(ref[[label_col]][, 1])

cat("\n TransferData.. .\n")
preds <- TransferData(
  anchorset = anchors,
  refdata = ref_labels,
  prediction.assay = FALSE,
  dims = use_dims,
  k.weight = k_wt,
  verbose = TRUE
)

cat("Done\n")

query_full <- readRDS("AT_Multiome_RootTip.RNA.rds")

meta_to_add <- data.frame(time_anno_fine_pred = preds$predicted.id,
                          row.names = rownames(preds))

meta_to_add <- meta_to_add[colnames(query_full), , drop = FALSE]
query_full <- AddMetaData(query_full, metadata = meta_to_add)

print(table(query_full$time_anno_fine_pred))

saveRDS(query_full, "root_stress_time_anno_new.rds")

crosstab <- table(query_full$Cell.type.RNA, query_full$time_anno_fine_pred)
crosstab_matrix <- as.matrix(crosstab)
crosstab_prop <- prop.table(crosstab_matrix, margin = 1) * 100
crosstab_prop[is.na(crosstab_prop)] <- 0
crosstab_prop[is.infinite(crosstab_prop)] <- 0

col_fun_prop <- colorRamp2(c(0, 50, 100), 
                           c("#F7FBFF", "#6BAED6", "#08306B"))

ht_prop <- Heatmap(crosstab_prop,
                   name = "Percentage (%)",
                   col = col_fun_prop,
                   
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.1f%%", crosstab_prop[i, j]), 
                              x, y, gp = gpar(fontsize = 9, col = ifelse(crosstab_prop[i, j] > 50, "white", "black")))
                   },
                   
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_method_rows = "complete",
                   
                   row_names_gp = gpar(fontsize = 11, fontface = "bold"),
                   column_names_gp = gpar(fontsize = 11, fontface = "bold"),
                   column_names_rot = 45,
                   border = TRUE,
                   
                   column_title = "Cell Type vs Time Stage (Percentage with Clustering)",
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                   
                   heatmap_legend_param = list(
                     title = "Percentage (%)",
                     title_gp = gpar(fontsize = 11, fontface = "bold"),
                     labels_gp = gpar(fontsize = 10),
                     legend_height = unit(5, "cm"),
                     at = c(0, 25, 50, 75, 100),
                     border = "black"
                   )
)

pdf("Cell_type_vs_time_heatmap_prop_cluster.pdf", width = 10, height = 8)
draw(ht_prop, heatmap_legend_side = "right")
dev.off()