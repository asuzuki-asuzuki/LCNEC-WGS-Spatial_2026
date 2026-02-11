library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

options(future.globals.maxSize = 100 * 1024^3)

set.seed(1234)
packageVersion("Seurat")

# Load a Seurat object
lung.merge <- readRDS("lung_visium_merge_harmony.rds")

# Spatial assay
DefaultAssay(lung.merge) <- "Spatial"
lung.merge

# Join layers
lung.merge[["Spatial"]] <- JoinLayers(lung.merge[["Spatial"]])

# DEG analysis
lung.markers <- FindAllMarkers(lung.merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10_merge_harmony.csv",append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)
write.table(lung.markers, "all_merge_harmony.csv",append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

rm(list = ls())
gc();gc()
