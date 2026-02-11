library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)

set.seed(1234)

packageVersion("Seurat")

# Load output of Cell Ranger
cell.data <- Read10X(data.dir = "<DIR>/filtered_feature_bc_matrix")

# Creat a Seurat object
cell <- CreateSeuratObject(counts = cell.data, project = "LCNEC", min.cells = 0, min.features = 0)

# QC
cell[["percent.mt"]] <- PercentageFeatureSet(cell, pattern = "^MT-")
length(cell@meta.data$nFeature_RNA)

pdf("VlnPlot_QC.pdf", width = 10, height = 5)
VlnPlot(cell, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pdf("FeatureScatter.pdf", width = 10, height = 5)
plot1 <- FeatureScatter(cell, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

rm(list = ls())
gc();gc()
