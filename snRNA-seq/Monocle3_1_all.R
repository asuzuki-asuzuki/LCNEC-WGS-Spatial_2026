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

# Load a merged Seurat object
lung <- readRDS("lung_merge_annotation.rds")

# RNA assay
DefaultAssay(lung) <- "RNA"
lung

# Join layers
lung[["RNA"]] <- JoinLayers(lung[["RNA"]])

# Preprocessing
genedata <- as.data.frame(x = rownames(lung), row.names = rownames(lung))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(lung[["RNA"]]$counts, cell_metadata = lung[[]][colnames(lung[["RNA"]]),], gene_metadata = genedata)
rm(lung)

# Calculation of %mitochondorial genes
mito_genes <- grep("^MT-", rownames(cds), value = TRUE, ignore.case = TRUE)
colData(cds)$pct_mito <- Matrix::colSums(exprs(cds)[mito_genes, ]) / Matrix::colSums(exprs(cds)) * 100

# Construction of trajectories
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "original", residual_model_formula_str = "~ pct_mito")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Set colors for annotation
COL <- c("pink", "#FF4FB5", "orange", "yellow", "green", "blue", "brown", "purple", "grey")
names(COL) <- c("Tumor cell", "Alveolar epithelial cell", "Fibroblast", "Endothelial cell", "T cell", "B cell", "Plasma cell", "Myeloid cell", "Others")

# Set colors for 6 samples
COL2 <- DiscretePalette(6, palette = "parade", shuffle = FALSE)
names(COL2) <- c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")

# Plots
pdf("monocle3_1_all.pdf", width = 12, height = 10)
p1 <- plot_cells(cds, color_cells_by = "partition")
p2 <- plot_cells(cds, color_cells_by = "annotation", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = TRUE, label_principal_points = TRUE) + scale_colour_manual(values = COL)
p3 <- plot_cells(cds, color_cells_by = "original", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + scale_colour_manual(values = COL2)
wrap_plots(p1, p2, p3, ncol = 2)
dev.off()

#Save the object
saveRDS(cds, file = "lung_monocle3_1_all.rds")

rm(list = ls())
gc();gc()
