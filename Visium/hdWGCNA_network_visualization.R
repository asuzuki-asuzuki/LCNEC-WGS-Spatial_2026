library(Seurat) # v5.0.3
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA) # v0.3.3

theme_set(theme_cowplot())
set.seed(1234)

packageVersion("Seurat")
packageVersion("hdWGCNA")

# Load a object
lung <- readRDS("lung_hdwgcna.rds")
lung

# UMAP for network modules
lung <- RunModuleUMAP(lung, n_hubs = 3, n_neighbors=15, min_dist=0.3, spread=1)

# Module UMAP plot
pdf("networkplot.pdf")
ModuleUMAPPlot(lung, edge.alpha=0.5, sample_edges=TRUE, keep_grey_edges=FALSE, edge_prop=0.075, label_hubs=3)
dev.off()

pdf("networkplot_2.pdf")
umap_df <- GetModuleUMAP(lung)
Name <- umap_df$gene
Name[!Name %in% umap_df$gene[umap_df$hub == "hub"]] <- ""
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label=Name)) + geom_point(color=umap_df$color, size=umap_df$kME*2) + umap_theme() + geom_text_repel(max.overlaps = Inf)
dev.off()

# Network visualization of representative module members in each module
ModuleNetworkPlot(lung, outdir = "ModuleNetworks")

rm(list = ls())
gc();gc()
