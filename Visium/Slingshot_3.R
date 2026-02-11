library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(monocle3)
library(magrittr)
library(slingshot) # v2.10.0
library(viridis)
library(fields)

set.seed(1234)

packageVersion("Seurat")
packageVersion("slingshot")

# Load a Seurat object
lung <- readRDS("lung_sub_merge_harmony_tumor.rds")

# Spatial assay
DefaultAssay(lung) <- "Spatial"

# Join layers
lung[["Spatial"]] <- JoinLayers(lung[["Spatial"]])
lung

# Load a object
sce <- readRDS("lung_slingshot_1_harmony_tumor.rds")

# Set colors
DATA <- read.table("id_type_color.txt", header = T, comment.char = "") # please prepare a list of sample names, subtypes and colors
# Set colors for 19 samples
COL <- DiscretePalette(19, palette = "polychrome", shuffle = FALSE)
names(COL) <- sort(DATA$ID)

# Plot for a sample
COL_ok <- rep("grey", length(COL))
names(COL_ok) <- sort(DATA$ID)
COL_ok["<SAMPLE>"] <- DATA$TypeCol[which(DATA$ID == "<SAMPLE>")] # please specify a sample

foreground <- c("<SAMPLE>") # please specify a sample
is_foreground <- lung$original %in% foreground
idx <- order(is_foreground)
sub <- lung$original

pdf("slingshot_tumor_<SAMPLE>.pdf", width = 7, height = 7) # please specify a sample
plot(reducedDims(sce)$UMAP.HARMONY[idx,], col = COL_ok[sub[idx]], cex = 0.5, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

rm(list = ls())
gc();gc()
