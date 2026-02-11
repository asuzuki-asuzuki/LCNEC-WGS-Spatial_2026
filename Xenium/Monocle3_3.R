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

# Load a CDS object
cds <- readRDS("lung_monocle3_2.rds")
cds

# Marker gene expression plots
pdf("monocle3_marker.pdf")
plot_cells(cds, genes=c("SFTPC", "YAP1", "ASCL1", "ASCL2"), label_cell_groups = FALSE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, trajectory_graph_segment_size = 0.5)
dev.off()


#######################
# Plots of Mixed cases#
#######################
# Set sample names of the Mixed type
Sample <- c("<SAMPLE1>", "<SAMPLE2>", "<SAMPLE3>", "<SAMPLE4>", "<SAMPLE5>", "<SAMPLE6>")

# Set colors
ANNO <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others", "No")
COL <- c(DiscretePalette(14, palette = "glasbey", shuffle = FALSE), "grey")
names(COL) <- ANNO

cds$OK <- "No"
cds$OK[which(cds$original %in% Sample)] <- as.vector(cds$annotation[which(cds$original %in% Sample)])

SUB <- c("NE-1", "NE-2", "NE-3", "NE-4", "NE-5", "NE-6", "NE-7", "NE-8", "NE-9", "Non-NE-p-1", "Non-NE-p-2", "Non-NE-p-3", "Non-NE-p-4", "Non-NE-p-5", "Non-NE-y-1", "Non-NE-y-2", "Non-NE-y-3", "NE/Non-NE-p-1", "NE/Non-NE-y-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3", "Alveolar-1", "Alveolar-2", "Bronchiolar-1", "Basal-1", "Fibrotic-1", "T cell-1", "T cell-2", "T cell-3", "T cell-4", "T cell-5", "T cell-6", "T cell-7", "B cell-1", "Plasma cell-1", "Plasma cell-2", "Plasma cell-3", "Plasma cell-4", "Macrophage-1", "Macrophage-2", "Macrophage-3", "Macrophage-4", "Macrophage-5", "Granulocyte-1", "Fibroblast-1", "Fibroblast-2",
"Fibroblast-3", "Fibroblast-4", "Fibroblast-5", "Fibroblast-6", "Fibroblast-7", "Endothelial cell-1", "Endothelial cell-2", "Endothelial cell-3", "Endothelial cell-4", "LEC-1", "SMC-1", "Others", "No")
COL2 <- c("#0000FF", "#33CCFF", "#000099", "#0099FF", "#003399", "#3366FF", "#3399FF", "#6699FF", "#0066FF", "#FF0000", "#660000", "#990000", "#CC3300", "#FF0033", "#FF6600", "#FF9900", "#FFCC00", "#663399", "#9900FF", "#660099", "#6600CC", "#FF00CC", "#FF0099", "#330033", "#00FF00", "#99FF33", "#33FFFF", "#99FFFF", "#66FFFF", "#00FFFF", "#00FCCF", "#99FCCF", "#33FF99", "#7742c0", "#109093", "#006666", "#336666", "#669999", "#FFCCFF", "#FF99FF", "#FF99CC", "#FF999F", "#FFCCCC", "#663300", "#006600", "#003300", "#336633", "#336600", "#006633", "#009933", "#339900",  "#00CCFF", "#0099FF", "#66CCFF", "#3399FF", "#a16057", "#FFD300", "#B1CC71", "grey")
names(COL2) <- SUB

cds$OK2 <- "No"
cds$OK2[which(cds$original %in% Sample)] <- as.vector(cds$subannotation[which(cds$original %in% Sample)])

# Plots
pdf("monocle3_mixed.pdf", width = 13, height = 10)
p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "OK", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL) + theme(legend.position = "none")
p3 <- plot_cells(cds, color_cells_by = "OK2", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL2) + theme(legend.position = "none")
p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
wrap_plots(p1, p2, p3, p4, ncol = 2)
dev.off()


######################
# Plots of SCLC cases#
######################
# Set sample names of SCLC cases
Sample <- c("<SAMPLE1>", "<SAMPLE2>", "<SAMPLE3>", "<SAMPLE4>", "<SAMPLE5>", "<SAMPLE6>")

# Set colors
ANNO <- c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others", "No")
COL <- c(DiscretePalette(14, palette = "glasbey", shuffle = FALSE), "grey")
names(COL) <- ANNO

cds$OK <- "No"
cds$OK[which(cds$original %in% c("tkbhgnec26", "tkbhgnec21", "tkbhgnec15", "tkbhgnec17", "tkbhgnec18", "tkbhgnec24", "tkbhgnec19-a", "tkbhgnec19-b"))] <- as.vector(cds$annotation[which(cds$original %in% c("tkbhgnec26", "tkbhgnec21", "tkbhgnec15", "tkbhgnec17", "tkbhgnec18", "tkbhgnec24", "tkbhgnec19-a", "tkbhgnec19-b"))])

SUB <- c("NE-1", "NE-2", "NE-3", "NE-4", "NE-5", "NE-6", "NE-7", "NE-8", "NE-9", "Non-NE-p-1", "Non-NE-p-2", "Non-NE-p-3", "Non-NE-p-4", "Non-NE-p-5", "Non-NE-y-1", "Non-NE-y-2", "Non-NE-y-3", "NE/Non-NE-p-1", "NE/Non-NE-y-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3", "Alveolar-1", "Alveolar-2", "Bronchiolar-1", "Basal-1", "Fibrotic-1", "T cell-1", "T cell-2", "T cell-3", "T cell-4", "T cell-5", "T cell-6", "T cell-7", "B cell-1", "Plasma cell-1", "Plasma cell-2", "Plasma cell-3", "Plasma cell-4", "Macrophage-1", "Macrophage-2", "Macrophage-3", "Macrophage-4", "Macrophage-5", "Granulocyte-1", "Fibroblast-1", "Fibroblast-2",
"Fibroblast-3", "Fibroblast-4", "Fibroblast-5", "Fibroblast-6", "Fibroblast-7", "Endothelial cell-1", "Endothelial cell-2", "Endothelial cell-3", "Endothelial cell-4", "LEC-1", "SMC-1", "Others", "No")
COL2 <- c("#0000FF", "#33CCFF", "#000099", "#0099FF", "#003399", "#3366FF", "#3399FF", "#6699FF", "#0066FF", "#FF0000", "#660000", "#990000", "#CC3300", "#FF0033", "#FF6600", "#FF9900", "#FFCC00", "#663399", "#9900FF", "#660099", "#6600CC", "#FF00CC", "#FF0099", "#330033", "#00FF00", "#99FF33", "#33FFFF", "#99FFFF", "#66FFFF", "#00FFFF", "#00FCCF", "#99FCCF", "#33FF99", "#7742c0", "#109093", "#006666", "#336666", "#669999", "#FFCCFF", "#FF99FF", "#FF99CC", "#FF999F", "#FFCCCC", "#663300", "#006600", "#003300", "#336633", "#336600", "#006633", "#009933", "#339900",  "#00CCFF", "#0099FF", "#66CCFF", "#3399FF", "#a16057", "#FFD300", "#B1CC71", "grey")
names(COL2) <- SUB

cds$OK2 <- "No"
cds$OK2[which(cds$original %in% c("tkbhgnec26", "tkbhgnec21", "tkbhgnec15", "tkbhgnec17", "tkbhgnec18", "tkbhgnec24", "tkbhgnec19-a", "tkbhgnec19-b"))] <- as.vector(cds$subannotation[which(cds$original %in% c("tkbhgnec26", "tkbhgnec21", "tkbhgnec15", "tkbhgnec17", "tkbhgnec18", "tkbhgnec24", "tkbhgnec19-a", "tkbhgnec19-b"))])

# Plots
pdf("monocle3_SCLC.pdf", width = 13, height = 10)
p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "OK", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL) + theme(legend.position = "none")
p3 <- plot_cells(cds, color_cells_by = "OK2", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE) + scale_colour_manual(values = COL2) + theme(legend.position = "none")
p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
wrap_plots(p1, p2, p3, p4, ncol = 2)
dev.off()

rm(list = ls())
gc();gc()
