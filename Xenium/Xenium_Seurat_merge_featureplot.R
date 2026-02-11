library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)

set.seed(1234)

packageVersion("Seurat")

# Load a merged Seurat object
lung <- readRDS("lung_xenium_merge_sketch_ProjectData.rds")

# Sketch assay
DefaultAssay(lung) <- "sketch"
lung

pdf("featureplot_sketch_raw.pdf", width=32, height=28)
FeaturePlot(lung, features = c("ASCL1", "CHGA", "CHGB", "YAP1", "FOXI1", "KIT", "KRT8", "KRT5", "KRT14", "SCGB1A1", "COL1A1", "PECAM1", "CD3E", "CD19", "CD68", "MZB1"), max.cutoff = "q90", slot="counts", reduction = "umap", ncol=4)
dev.off()

pdf("featureplot_sketch_raw_2.pdf", width=32, height=28)
FeaturePlot(lung, features = c("PTPRC", "CD163", "ACTA2", "CDH1", "SCGB3A2", "SFTPC", "SFTPB", "SFTPD", "CALCA", "TP63", "CCNO", "HEPACAM2", "MUC5AC", "MUC5B", "FOXJ1"), max.cutoff = "q90", slot="counts", reduction = "umap", ncol=4)
dev.off()


# Xenium assay (full dataset)
DefaultAssay(lung) <- "Xenium"
lung

pdf("featureplot_full_raw.pdf", width=32, height=28)
FeaturePlot(lung, features = c("ASCL1", "CHGA", "CHGB", "YAP1", "FOXI1", "KIT", "KRT8", "KRT5", "KRT14", "SCGB1A1", "COL1A1", "PECAM1", "CD3E", "CD19", "CD68", "MZB1"), slot = "counts", max.cutoff = "q90", reduction = "full.umap", ncol=4, alpha=0.1)
dev.off()

pdf("featureplot_full_raw_2.pdf", width=32, height=28)
FeaturePlot(lung, features = c("PTPRC", "CD163", "ACTA2", "CDH1", "SCGB3A2", "SFTPC", "SFTPB", "SFTPD", "CALCA", "TP63", "CCNO", "HEPACAM2", "MUC5AC", "MUC5B", "FOXJ1"), slot = "counts", max.cutoff = "q90", reduction = "full.umap", ncol=4, alpha=0.1)
dev.off()

rm(list = ls())
gc();gc()
