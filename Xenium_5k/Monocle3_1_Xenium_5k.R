library(Seurat) # v5.4.0
library(ggplot2)
library(patchwork)
library(dplyr)

library(SeuratWrappers)
library(monocle3)
library(magrittr)

set.seed(1234)

packageVersion("Seurat")
packageVersion("monocle3")

# Load a Seurat object (epithelial cells)
lung <- readRDS("lung_sub.rds")

# Sketch assay
DefaultAssay(lung) <- "sketch"
lung

# Set cds
genedata <- as.data.frame(x = rownames(lung), row.names = rownames(lung))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(lung[["sketch"]]$counts, cell_metadata = lung[[]][colnames(lung[["sketch"]]),], gene_metadata = genedata)

rm(lung)

# Construction of trajectories
cds <- preprocess_cds(cds, num_dim = 10)
cds <- align_cds(cds, alignment_group = "original")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Set colors for merge_annotation
COL <- DiscretePalette(14, palette = "glasbey", shuffle = FALSE)
COL2 <- c(COL[1:13], "brown", COL[14])
names(COL2) <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC/Pericyte", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Mastocyte", "Others")

# Set colors for 13 cases
big_palette <- c(
  DiscretePalette(78, palette = "parade"),
  DiscretePalette(26, palette = "alphabet"),
  DiscretePalette(26, palette = "alphabet2"),
  DiscretePalette(36, palette = "polychrome"),
  DiscretePalette(32, palette = "glasbey")
)

COL3 <- big_palette[1:13]

# Set colors for annotation in each section
COL4 <- c(DiscretePalette(14, palette = "glasbey", shuffle = FALSE), "cyan3", "darkorange", "hotpink", "cyan4")
names(COL4) <- c("Tumor cell, NE", "Tumor cell, non-NE", "Basal epithelium", "Bronchiolar epithelium", "Alveolar epithelium", "Fibroblast", "SMC", "EC", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others", "DC", "Mastocyte", "Pericyte", "NK cell")

# Plots
pdf("monocle3_1.pdf", width = 13, height = 10)
p1 <- plot_cells(cds, color_cells_by = "partition")
p2 <- plot_cells(cds, color_cells_by = "annotation", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_principal_points = TRUE) + scale_colour_manual(values = COL4)
p3 <- plot_cells(cds, color_cells_by = "merge_annotation", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + scale_colour_manual(values = COL2)
p4 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + scale_colour_manual(values = COL3)
wrap_plots(p1, p2, p3, p4, ncol = 2)
dev.off()

# Save the object
saveRDS(cds, file = "lung_monocle3_1.rds")

rm(list = ls())
gc();gc()
