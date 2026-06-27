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

# Load a CDS object
cds <- readRDS("lung_monocle3_1.rds")
cds

# Plots of expression levels on the trajectories
pdf("monocle3_marker-1.pdf")
plot_cells(cds, genes=c("SFTPC", "SFTPB", "AGER", "HOPX"), label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_segment_size = 0.5)
dev.off()

pdf("monocle3_marker-2.pdf")
plot_cells(cds, genes=c("ASCL1", "NEUROD1", "POU2F3", "YAP1"), label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_segment_size = 0.5)
dev.off()


# Calculate pseudotime
cds <- order_cells(cds, root_pr_nodes='Y_72')


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

# Set colors for subtypes (ANPY)
COL4 <- c("#FFC000", "#FF0000", "#00B0F0", "#00B050", "#FF0000", "#FF0000", "#00B0F0", "#00B0F0", "#FF0000", "#FF0000", "#00B050", "#FF0000", "#00B050")
names(COL4) <- c(<Sample IDs>)

# Set colors for subtypes
COL5 <- c("#FF0000", "#002060", "#9933FF", "#C00000", "#0000FF", "#6600CC", "#002060", "#0000FF", "#0000FF", "#002060", "#C00000", "#0000FF", "#002060")
names(COL5) <- c(<Sample IDs>)


# Plots
pdf("monocle3_2.pdf", width = 13, height = 15)
p1 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL4)
p2 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL5)
p3 <- plot_cells(cds, color_cells_by = "merge_annotation", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL2)
p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
p5 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL3)
wrap_plots(p1, p2, p3, p4, p5, ncol = 2)
dev.off()

# Save the object
saveRDS(cds, file = "lung_monocle3_2.rds")

rm(list = ls())
gc();gc()
