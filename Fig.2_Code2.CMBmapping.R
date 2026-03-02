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

tenxRDS <- readRDS("Processed_tomato_10X.rds")
bmkRDS <- readRDS("Processed_tomato_BMK.rds")

meta_10x_csv <- read.csv("10x_cmb_27_08.csv", row.names = 1)
meta_bmk_csv <- read.csv("bmk_cmb_63_05.csv", row.names = 1)

tenxRDS <- AddMetaData(tenxRDS, metadata = meta_10x_csv)
head(tenxRDS@meta.data)
bmkRDS <- AddMetaData(bmkRDS, metadata = meta_bmk_csv)
head(bmkRDS@meta.data)

bmkRDS <- UpdateSeuratObject(bmkRDS)
if (!"misc" %in% slotNames(bmkRDS)) {
  bmkRDS@misc <- list()
}

bmkRDS <- SCTransform(
  bmkRDS,
  assay = "Spatial",
  new.assay.name = "SCT",
  verbose = FALSE
)
DefaultAssay(bmkRDS) <- "SCT"

counts_10x <- LayerData(tenxRDS, layer = "counts")
rownames(counts_10x) <- gsub("^gene:", "", rownames(counts_10x))
tenxRDS_clean <- CreateSeuratObject(
  counts = counts_10x,
  meta.data = tenxRDS@meta.data,
  assay = "RNA"
)
tenxRDS_clean <- SCTransform(tenxRDS_clean, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
DefaultAssay(tenxRDS_clean) <- "SCT"

ref_obj   <- tenxRDS_clean
query_obj <- bmkRDS

ref_obj$Broad_Type   <- ref_obj$annotation   
query_obj$Broad_Type <- query_obj$annotation 

common_types <- intersect(
  na.omit(unique(ref_obj$Broad_Type)),
  na.omit(unique(query_obj$Broad_Type))
)

common_features <- intersect(
  VariableFeatures(ref_obj),
  VariableFeatures(query_obj)
)

npcs_ref   <- min(50, ncol(ref_obj) - 1, length(common_features) - 1)
npcs_query <- min(50, ncol(query_obj) - 1, length(common_features) - 1)

ref_obj   <- RunPCA(ref_obj, features = common_features, npcs = npcs_ref, approx = FALSE, verbose = FALSE)
query_obj <- RunPCA(query_obj, features = common_features, npcs = npcs_query, approx = FALSE, verbose = FALSE)

all_predictions <- list()

for (ctype in common_types) {
  
  ref_subset   <- subset(ref_obj, Broad_Type == ctype)
  query_subset <- subset(query_obj, Broad_Type == ctype)
  
  max_dims <- min(30,
                  ncol(Embeddings(ref_subset, "pca")),
                  ncol(Embeddings(query_subset, "pca")))
  use_dims <- 1:max_dims
  
  anchors <- FindTransferAnchors(
    reference = ref_subset,
    query = query_subset,
    normalization.method = "SCT",
    reference.assay = DefaultAssay(ref_subset),
    query.assay = DefaultAssay(query_subset),
    dims = use_dims,
    reduction = "rpca",
    verbose = FALSE
  )
  
  n_anchors <- nrow(anchors@anchors)
  if (is.null(n_anchors) || n_anchors < 3) {
    message(paste("  -Done :", ctype, "（anchors=", n_anchors, "）"))
    next
  }
  
  k_wt <- max(1, min(10, n_anchors - 1)) 
  
  preds <- TransferData(
    anchorset = anchors,
    refdata = ref_subset$CMB,
    prediction.assay = FALSE,
    dims = use_dims,
    k.weight = k_wt,
    verbose = FALSE
  )
  
  all_predictions[[ctype]] <- preds
  message(paste("  -Done :", ctype, "（anchors=", n_anchors, ", k.weight=", k_wt, ")"))
}

pred_list <- lapply(all_predictions, function(df) {
  df[, "predicted.id", drop = FALSE]
})

final_predictions <- do.call(rbind, pred_list)

rownames(final_predictions) <- sub("^[^.]+\\.", "", rownames(final_predictions))
aligned_preds <- data.frame(
  row.names    = colnames(query_obj),
  predicted.id = rep(NA_character_, length(colnames(query_obj)))
)
ov <- intersect(colnames(query_obj), rownames(final_predictions))
aligned_preds[ov, "predicted.id"] <- final_predictions[ov, "predicted.id"]
bmkRDS <- AddMetaData(bmkRDS, metadata = aligned_preds)
message("  -Done   ", length(ov))
