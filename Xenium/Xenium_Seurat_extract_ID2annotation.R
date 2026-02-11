library(Seurat) # v5.0.3
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1234)

packageVersion("Seurat")
options(future.globals.maxSize = 100 * 1024^3)

# Load a Seurat object
lung.merge <- readRDS("lung_xenium_merge_sketch_ProjectData_annotation.rds")

data <- data.frame(case_id = as.vector(lung.merge$original), cell_id = names(lung.merge$annotation), group=as.vector(lung.merge$annotation))
write.table(data, file = "extract_ID2Cluster.txt", sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = ls())
gc();gc()
