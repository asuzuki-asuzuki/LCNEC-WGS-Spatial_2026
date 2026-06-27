library(Seurat) # v5.4.0
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1234)

# Load the object
lung <- readRDS("lung_xenium.rds")

#Set raw counts to data layer for visualization of raw counts
lung[["Xenium"]]$data <- lung[["Xenium"]]$counts

DefaultAssay(lung) <- "Xenium"

SIZE <- 0.1

# Spatial plots
pdf("marker_basic.pdf", width = 10, height = 10)
p = ImageFeaturePlot(lung, features = c("ASCL1", "NEUROD1", "POU2F3", "YAP1"), size = SIZE, cols = c("white", "red"), border.size = NA, axes = FALSE, max.cutoff = "q90") & NoGrid()
p + plot_layout(ncol = 2, widths = c(1, 1), heights = c(1, 1))
dev.off()

rm(list = ls())
gc();gc()
