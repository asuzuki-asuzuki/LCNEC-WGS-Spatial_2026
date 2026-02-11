library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(1234)

packageVersion("Seurat")

options(future.globals.maxSize = 100 * 1024^3)

# Load a merged Seurat object
lung.merge <- readRDS("lung_sub_sketch_ProjectData.rds")

# Xenium assay
DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# 26 subclusters
lung.merge$cluster_full <- factor(x = lung.merge$cluster_full, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25"))

# Projection of full data to the sketch assay
Idents(object = lung.merge) <- "cluster_full"
new.clutser.id <- c("Non-NE-p-1", "NE-1", "Non-NE-p-2", "NE-2", "NE-3", "NE-4", "Non-NE-p-3", "Alveolar-1", "Alveolar-2", "NE-5", "NE-6", "NE/Non-NE-y-1", "Non-NE-y-1", "Non-NE-p-4", "NE-7", "Non-NE-p-5", "NE-8", "Non-NE-y-2", "NE/Non-NE-p-1", "NE-9", "Fibrotic-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3", "Bronchiolar-1", "Non-NE-y-3", "Basal-1")
names(new.clutser.id) <- levels(lung.merge)
lung.merge <- RenameIdents(lung.merge, new.clutser.id)

# Set metadata "subannotation"
lung.merge[["subannotation"]] <- Idents(lung.merge)

## Set colors for 26 cell subtypes (subannotation)
lung.merge$subannotation <- factor(x = lung.merge$subannotation, levels = c("NE-1", "NE-2", "NE-3", "NE-4", "NE-5", "NE-6", "NE-7", "NE-8", "NE-9", "Non-NE-p-1", "Non-NE-p-2", "Non-NE-p-3", "Non-NE-p-4", "Non-NE-p-5", "Non-NE-y-1", "Non-NE-y-2", "Non-NE-y-3", "NE/Non-NE-p-1", "NE/Non-NE-y-1", "NE/Non-NE-y-2", "NE/Non-NE-y-3", "Alveolar-1", "Alveolar-2", "Bronchiolar-1", "Basal-1", "Fibrotic-1"))
COL <- c("#0000FF", "#33CCFF", "#000099", "#0099FF", "#003399", "#3366FF", "#3399FF", "#6699FF", "#0066FF", "#FF0000", "#660000", "#990000", "#CC3300", "#FF0033", "#FF6600", "#FF9900", "#FFCC00", "#663399", "#9900FF", "#660099", "#6600CC", "#FF00CC", "#FF0099", "#330033", "#00FF00", "#99FF33")

# UMAP plot
pdf("full_UMAP_sub_cluster_subannotation.pdf", width = 9, height = 8)
p1 <- DimPlot(lung.merge, reduction = "full.sub_umap", group.by = c("subannotation"), label = TRUE, alpha = 0.1, repel=TRUE, cols = COL)
p1
dev.off()

# Save the object
saveRDS(lung.merge, file = "lung_sub_sketch_ProjectData_subannotation.rds")

# Count cells for each subtype
id <- unique(lung.merge$original)
Cluster <- names(table(lung.merge$subannotation))
OK <- data.frame(Cluster)

for(i in 1:37) {
      tt <- table(lung.merge$subannotation[which(lung.merge$original == id[i])])
      OK[,id[i]] <- as.vector(tt)
}

write.table(OK, file = "numCell_subannotation.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

rm(list = ls())
gc();gc()
