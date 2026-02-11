library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(monocle3) # v1.3.7
library(magrittr)

set.seed(1234)

packageVersion("Seurat")
packageVersion("monocle3")

# Load a object
cds <- readRDS("lung_monocle3_1_all.rds")
cds

# Plots of marker expression patterns
pdf("monocle3_marker-1_all.pdf")
plot_cells(cds, genes=c("SFTPC", "SFTPB", "AGER", "HOPX"), label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_segment_size = 0.5)
dev.off()

pdf("monocle3_marker-2_all.pdf")
plot_cells(cds, genes=c("ASCL1", "NEUROD1", "POU2F3", "YAP1"), label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_segment_size = 0.5)
dev.off()

# Calculate pseudotime
cds <- order_cells(cds, root_pr_nodes='Y_146')

# Set colors for 6 samples
COL <- DiscretePalette(6, palette = "parade", shuffle = FALSE)
names(COL) <- c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")

# Set colors for 6 samples in each subtype (NE, non-NE (ASCL2), non-NE (YAP1), Mixed, DP)
COL2 <- c("#7030a0", "#7030a0", "#0000ff", "#7030a0", "#0000ff", "#c00000")
names(COL2) <- c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")

# Set colors for annotation
COL3 <- c("pink", "#FF4FB5", "orange", "yellow", "green", "blue", "brown", "purple", "grey")
names(COL3) <- c("Tumor cell", "Alveolar epithelial cell", "Fibroblast", "Endothelial cell", "T cell", "B cell", "Plasma cell", "Myeloid cell", "Others")

# Plots
pdf("monocle3_2_all.pdf", width = 13, height = 10)
p1 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL)
p2 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL2)
p3 <- plot_cells(cds, color_cells_by = "annotation", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL3)
p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
wrap_plots(p1, p2, p3, p4, ncol = 2)
dev.off()

# Save the object
saveRDS(cds, file = "lung_monocle3_2_all.rds")

rm(list = ls())
gc();gc()
