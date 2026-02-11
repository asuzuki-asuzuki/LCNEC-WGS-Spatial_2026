library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung <- readRDS("lung_sub_sketch_ProjectData_subannotation.rds")

# Xenium assay
DefaultAssay(lung) <- "Xenium"
lung

# Set colors
COL <- c("#0000FF", "#33CCFF", "#000099", "#0099FF", "#003399", "#3366FF", "#3399FF", "#6699FF", "#0066FF", "#FF0000", "#660000", "#990000", "#CC3300", "#FF0033", "#FF6600", "#FF9900", "#FFCC00", "#663399", "#9900FF", "#660099", "#6600CC", "#FF00CC", "#FF0099", "#330033", "#00FF00", "#99FF33")

# Violin plots for ASCL1, ASCL2 and YAP1
pdf("violinplot_subannotation_full_raw_limit_ASCL1.pdf", width=9, height=5)
VlnPlot(lung, features = c("ASCL1"), y.max = 50, group.by = "subannotation", layer="counts", ncol=1, pt.size=0, cols = COL) & NoLegend()
dev.off()

pdf("violinplot_subannotation_full_raw_limit_ASCL2.pdf", width=9, height=5)
VlnPlot(lung, features = c("ASCL2"), y.max = 10, group.by = "subannotation", layer="counts", ncol=1, pt.size=0, cols = COL) & NoLegend()
dev.off()

pdf("violinplot_subannotation_full_raw_limit_YAP1.pdf", width=9, height=5)
VlnPlot(lung, features = c("YAP1"), y.max = 15, group.by = "subannotation", layer="counts", ncol=1, pt.size=0, cols = COL) & NoLegend()
dev.off()

rm(list = ls())
gc();gc()
