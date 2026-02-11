library(Seurat) # v5.0.3
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA) # v0.3.3

library(RColorBrewer)
library(scales)

theme_set(theme_cowplot())
set.seed(1234)

packageVersion("Seurat")
packageVersion("hdWGCNA")

# Load a object
lung <- readRDS("lung_hdwgcna.rds")
lung

# Extract module info
modules <- GetModules(lung)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# Spatial plots of MEs
for (i in 1:length(mods)) {
    pdf(paste0("hdwgcna_spatialplot/", "hdwgcna_spatialplot_", mods[i], ".pdf"), width = 16, height = 25)
    p <- SpatialFeaturePlot(lung, features = mods[i], alpha = c(1), ncol = 4, pt.size.factor = 2.5, image.alpha = 0, stroke = NA)
    p <- p & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(min(lung[[mods[i]]]), max(lung[[mods[i]]])), oob=squish)
    print(p)
    dev.off()
}

for (i in 1:length(mods)) {
    pdf(paste0("hdwgcna_spatialplot/", "hdwgcna_spatialplot_image_", mods[i], ".pdf"), width = 16, height = 25)
    p <- SpatialFeaturePlot(lung, features = mods[i], alpha = c(1), ncol = 4, pt.size.factor = 2.5, stroke = NA)
    p <- p & scale_fill_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(min(lung[[mods[i]]]), max(lung[[mods[i]]])), oob=squish)
    print(p)
    dev.off()
}

rm(list = ls())
gc();gc()
