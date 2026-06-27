library(Seurat) # v5.4.0
library(dplyr)
library(patchwork)
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3) 
set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung <- readRDS("lung_xenium_merge_sketch_ProjectData_annotation.rds")

# extraction of epitehlial cell sub-clusters
lung <- subset(lung, subset = merge_annotation %in% c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar"))

DefaultAssay(lung) <- "Xenium"
lung

# Save the object
saveRDS(lung, file = "lung_sub.rds")

rm(list = ls())
gc();gc()
