library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(reshape2)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object for the sample
lung <- readRDS("lung_annotation.rds")
lung

# Subset microenvironmental cells
lung <- subset(lung, subset = annotation %in% c("Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others"))

# Annotation and subannotation
ANNO <- c("Fibroblast", "SMC", "Endothelial cell", "LEC", "T cell", "B cell", "Plasma cell", "Macrophage", "Others")
lung$annotation <- factor(x = lung$annotation, levels = ANNO)
SUB <- c("T cell-1", "T cell-2", "T cell-3", "T cell-4", "T cell-5", "T cell-6", "T cell-7", "B cell-1", "Plasma cell-1", "Plasma cell-2", "Plasma cell-3", "Plasma cell-4", "Macrophage-1", "Macrophage-2", "Macrophage-3", "Macrophage-4", "Macrophage-5", "Granulocyte-1", "Fibroblast-1", "Fibroblast-2", "Fibroblast-3", "Fibroblast-4", "Fibroblast-5", "Fibroblast-6", "Fibroblast-7", "Endothelial cell-1", "Endothelial cell-2", "Endothelial cell-3", "Endothelial cell-4", "LEC-1", "SMC-1", "Others")
lung$subannotation <- factor(x = lung$subannotation, levels = SUB)

# constructing a cellular niche matrix and conducting K-means clustering
lung <- BuildNicheAssay(object = lung, fov = "fov", group.by = "annotation", niches.k = 12, neighbors.k = 25)
lung$niches <- factor(x = lung$niches, levels = sort(as.numeric(unique(lung$niches))))

# Spatial plot of cellular niche patterns
pdf("spatialplot_niche_TME.pdf", width = 9, height = 8)
p1 <- ImageDimPlot(lung, group.by = "niches", size = 0.3, dark.background = F, cols = "parade", border.size = NA, axes = TRUE) + NoGrid()
p1
dev.off()

# Output of tables about cellular niche and annotation/subannotation
table(lung$subannotation, lung$niches)
data <- table(lung$subannotation, lung$niches)
write.table(data, file = "niche_subannotation_TME.txt", sep = "\t", col.names = NA, row.names = T, append = F, quote = F)

table(lung$annotation, lung$niches)
data <- table(lung$annotation, lung$niches)
write.table(data, file = "niche_annotation_TME.txt", sep = "\t", col.names = NA, row.names = T, append = F, quote = F)
data2 <- melt(data/sum(data), varnames=c("CellType", "Niche"), value.name = "Rate")

# Dotplot
N <- length(sort(as.numeric(unique(lung$niches))))
COL <- DiscretePalette(N, palette = "parade")

pdf("dot_plot_niche_TME.pdf", width = 7, height = 6)
ggplot(data2, aes(x = Niche, y =  fct_rev(CellType))) +  labs(x = "Niche", y = "Cell-type group") +
	      geom_point(aes(size = Rate, colour = as.factor(Niche))) +
	      scale_size_area(breaks = c(0, 0.01, 0.02, 0.05, 0.1), max_size = 10, oob = scales::squish, limits = c(0, 0.1), labels = c("0","0.01","0.02","0.05","0.1+")) +
	      scale_colour_manual(values = COL) +
	      theme_bw() +
	      theme(panel.grid=element_blank(), axis.ticks = element_blank()) +
	      scale_x_continuous(breaks = c(1:12))
dev.off()

# Save the object
saveRDS(lung, file = "lung_niche_TME.rds")

rm(list = ls())
gc();gc()
