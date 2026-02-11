library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)

library(RColorBrewer)
library(scales)
library(ggplot2)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
cell <- readRDS("lung_stomics_bin50.rds")

# Number of bins
length(cell@meta.data$nFeature_Spatial)

# Load the module member list
DATA <- read.table("module_genes_kME.txt", header=T) # the list of module member genes created by hdWGCNA (from Visium)

# Extraction of member genes in the DP-associated modules
Module17 <- DATA$gene_name[which(DATA$module == "Module17" & DATA$kME_Module17 > 0.4)]
Module17
Module13 <- DATA$gene_name[which(DATA$module == "Module13" & DATA$kME_Module13 > 0.4)]
Module13
Module11 <- DATA$gene_name[which(DATA$module == "Module11" & DATA$kME_Module11 > 0.4)]
Module11

# Extraction of member genes in the oxidative stress response/HNF4A module
Module9 <- DATA$gene_name[which(DATA$module == "Module9" & DATA$kME_Module9 > 0.4)]
Module9
# Extraction of member genes in the ASCL1 module
Module14 <- DATA$gene_name[which(DATA$module == "Module14" & DATA$kME_Module14 > 0.4)]
Module14

# Add module scores
cell <- AddModuleScore(cell, features = list(Module17), name = "Module17")
colnames(cell@meta.data) <- gsub(x = colnames(cell@meta.data), pattern = "Module171", replacement = "Module17")
cell <- AddModuleScore(cell, features = list(Module13), name = "Module13")
colnames(cell@meta.data) <- gsub(x = colnames(cell@meta.data), pattern = "Module131", replacement = "Module13")
cell <- AddModuleScore(cell, features = list(Module11), name = "Module11")
colnames(cell@meta.data) <- gsub(x = colnames(cell@meta.data), pattern = "Module111", replacement = "Module11")
cell <- AddModuleScore(cell, features = list(Module9), name = "Module9")
colnames(cell@meta.data) <- gsub(x = colnames(cell@meta.data), pattern = "Module91", replacement = "Module9")
cell <- AddModuleScore(cell, features = list(Module14), name = "Module14")
colnames(cell@meta.data) <- gsub(x = colnames(cell@meta.data), pattern = "Module141", replacement = "Module14")

# Module scores in UMAP and spatial plots
pdf("Module_projection_FeaturePlot.pdf", width = 15, height = 8)
FeaturePlot(cell, features = c("Module17", "Module13", "Module11", "Module9", "Module14"), ncol = 3) & scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(-0.1, 0.1), oob=squish)
dev.off()

pdf("Module_projection_SpatialFeaturePlot.pdf", width = 21, height = 14)
SpatialFeaturePlot(cell, features = c("Module17", "Module13", "Module11", "Module9", "Module14"), ncol = 3, image.alpha = 0, pt.size.factor = 70, stroke=NA) & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(-0.1, 0.1), oob=squish)
dev.off()

rm(list = ls())
gc();gc()
