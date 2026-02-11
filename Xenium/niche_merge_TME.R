library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(reshape2)

set.seed(1234)

packageVersion("Seurat")

# Set number of clusters
niche.k <- 20

# Load a merged Seurat object with the merged niche matrix
lung <- readRDS("lung_niche_TME_merge.rds")

# niche assay
DefaultAssay(lung) <- "niche"

# K-means clustering for the merged niche matrix
lung <- ScaleData(lung)
results <- kmeans(x = t(lung[["niche"]]@scale.data), centers = niche.k, nstart = 30)
lung[["mergedniches"]] <- results[["cluster"]]

# 20 cellular niches for the merged dataset
lung$mergedniches <- factor(x = lung$mergedniches, levels = sort(as.numeric(unique(lung$mergedniches))))


# Output of tables about merged niches in each annotation/subannotation
table(lung$subannotation, lung$mergedniches)
data <- table(lung$subannotation, lung$mergedniches)
write.table(data, file = "merged_niche_subannotation_TME.txt", sep = "\t", col.names = NA, row.names = T, append = F, quote = F)

table(lung$annotation, lung$mergedniches)
data <- table(lung$annotation, lung$mergedniches)
write.table(data, file = "merged_niche_annotation_TME.txt", sep = "\t", col.names = NA, row.names = T, append = F, quote = F)


# Dot plot
N <- length(sort(as.numeric(unique(lung$mergedniches))))
COL <- DiscretePalette(N, palette = "parade")
names(COL) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")

data <- table(lung$annotation, lung$mergedniches)
data2 <- melt(data/sum(data), varnames=c("CellType", "Niche"), value.name = "Rate")

pdf("dot_plot_merged_niche_TME_re.pdf", width = 7, height = 8)
ggplot(data2, aes(x = Niche, y =  fct_rev(CellType))) +  labs(x = "Niche", y = "Cell-type group") +
              geom_point(aes(size = Rate, colour = as.factor(Niche))) +
              scale_size_area(breaks = c(0, 0.01, 0.02, 0.05, 0.1), max_size = 10, oob = scales::squish, limits = c(0, 0.1), labels = c("0","0.01","0.02","0.05","0.1+")) +
              scale_colour_manual(values = COL) +
              theme_bw() +
              theme(panel.grid=element_blank(), axis.ticks = element_blank()) +
              scale_x_continuous(breaks = c(1:20))
dev.off()

saveRDS(lung, file = "lung_niche_TME_merge_kmeans.rds")

rm(list = ls())
gc();gc()
