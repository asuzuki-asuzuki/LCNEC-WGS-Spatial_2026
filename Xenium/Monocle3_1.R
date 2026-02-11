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

# Load a Seurat object
lung <- readRDS("lung_sub_sketch_ProjectData_subannotation.rds")

# Sketch assay
DefaultAssay(lung) <- "sketch"
lung

# Join layers
lung[["sketch"]] <- JoinLayers(lung[["sketch"]])
lung

# Set cds
genedata <- as.data.frame(x = rownames(lung), row.names = rownames(lung))
colnames(genedata) <- "gene_short_name"
cds <- new_cell_data_set(lung[["sketch"]]$counts, cell_metadata = lung[[]][colnames(lung[["sketch"]]),], gene_metadata = genedata)

rm (lung)

# Construction of trajectories
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Save the object
saveRDS(cds, file = "lung_monocle3_1.rds")

rm(list = ls())
gc();gc()
