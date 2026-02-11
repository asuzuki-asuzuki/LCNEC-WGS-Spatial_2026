library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a merged Seurat object
lung.merge <- readRDS("lung_xenium_merge_sketch_ProjectData.rds")

# Xenium assay
DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# 34 clusters
lung.merge$cluster_full <- factor(x = lung.merge$cluster_full, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33"))

# Assign annotations into clusters
Idents(object = lung.merge) <- "cluster_full"
new.clutser.id <- c("Others", "T cell", "Macrophage", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Plasma cell", "Fibroblast", "Endothelial cell", "Epithelial cell, NE tumor", "Epithelial cell, alveolar", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, bronchiolar", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "B cell", "Epithelial cell, non-NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "T cell", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "LEC", "SMC", "Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, NE tumor", "Epithelial cell, NE tumor", "Macrophage", "Epithelial cell, non-NE tumor")
names(new.clutser.id) <- levels(lung.merge)
lung.merge <- RenameIdents(lung.merge, new.clutser.id)

# Set metadata "annotation"
lung.merge[["annotation"]] <- Idents(lung.merge)

# Set colors for 37 sections
COL <- c(DiscretePalette(36, palette = "polychrome", shuffle = FALSE), "#0000FF")

# Set colors for 14 cell types (annotation)
lung.merge$annotation <- factor(x = lung.merge$annotation, levels = c("Epithelial cell, NE tumor", "Epithelial cell, non-NE tumor", "Epithelial cell, basal", "Epithelial cell, bronchiolar", "Epithelial cell, alveolar", "Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others"))
COL2 <- DiscretePalette(14, palette = "glasbey", shuffle = FALSE)

# UMAP plots
pdf("full_UMAP_merge_cluster_annotation.pdf", width = 18, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("annotation"), alpha = 0.1, cols = COL2)
p2 <- DimPlot(lung.merge, reduction = "full.umap", group.by = c("original"), alpha = 0.1, cols = COL)
p1 + p2
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge_sketch_ProjectData_annotation.rds")

# Count cells for each celltype
id <- unique(lung.merge$original)
Cluster <- names(table(lung.merge$annotation))
OK <- data.frame(Cluster)

for(i in 1:37) {
      tt <- table(lung.merge$annotation[which(lung.merge$original == id[i])])      
      OK[,id[i]] <- as.vector(tt)
}

write.table(OK, file = "numCell_annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

rm(list = ls())
gc();gc()
