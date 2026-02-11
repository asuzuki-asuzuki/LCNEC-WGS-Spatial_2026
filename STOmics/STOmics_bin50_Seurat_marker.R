library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
lung <- readRDS("lung_stomics_bin50.rds")

# Marker expression patterns in spatial plots
pdf("marker.pdf")
SpatialFeaturePlot(lung, features = c("ASCL1", "ASCL2", "YAP1", "SFTPC", "POU2F3", "NEUROD1"), ncol = 3, image.alpha = 0, pt.size.factor = 70, stroke=NA)
dev.off()

rm(list = ls())
gc();gc()
