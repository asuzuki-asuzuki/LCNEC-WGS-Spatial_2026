library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)
packageVersion("Seurat")

library(future)
plan("multisession", workers = 8)

options(future.globals.maxSize = 100 * 1024^3)

# Load a Seurat object
lung <- readRDS("lung_sub_sketch.rds")
lung

# Sketch assay
DefaultAssay(lung) <- "sketch"

# Join layers
lung[["sketch"]] <- JoinLayers(lung[["sketch"]])

# DEG analysis
lung.markers <- FindAllMarkers(lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- lung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "sub_top10_xenium.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)
write.table(lung.markers, "sub_all_xenium.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

rm(list = ls())
gc();gc()
