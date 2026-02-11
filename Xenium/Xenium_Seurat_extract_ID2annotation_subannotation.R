library(Seurat) # v5.0.3
library(patchwork)
library(dplyr)
library(ggplot2)

set.seed(1234)

packageVersion("Seurat")
options(future.globals.maxSize = 100 * 1024^3)

# Load a Seurat object with cell type annotation
lung.merge <- readRDS("lung_xenium_merge_sketch_ProjectData_annotation.rds")

data <- data.frame(case_id = as.vector(lung.merge$original), cell_id = names(lung.merge$annotation), annotation = as.vector(lung.merge$annotation))
write.table(data, file = "extract_ID2annotation.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Load a Seurat object with cell subtype annotation (Epithelial clusters)
lung.merge <- readRDS("<Epithelial DIR>/lung_sub_sketch_ProjectData_subannotation.rds")

data <- data.frame(case_id = as.vector(lung.merge$original), cell_id = names(lung.merge$subannotation), subannotation = as.vector(lung.merge$subannotation))
write.table(data, file = "extract_ID2subannotation_Epithelial.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# Load a Seurat object with cell subtype annotation (Stromal clusters)
lung.merge <- readRDS("<Stromal DIR>/lung_sub_sketch_ProjectData_subannotation.rds")

data <- data.frame(case_id = as.vector(lung.merge$original), cell_id = names(lung.merge$subannotation), subannotation = as.vector(lung.merge$subannotation))
write.table(data, file = "extract_ID2subannotation_Stromal.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Load a Seurat object with cell subtype annotation (Immune clusters)
lung.merge <- readRDS("<Immune DIR>/lung_sub_sketch_ProjectData_subannotation.rds")

data <- data.frame(case_id = as.vector(lung.merge$original), cell_id = names(lung.merge$subannotation), subannotation = as.vector(lung.merge$subannotation))
write.table(data, file = "extract_ID2subannotation_Immune.txt", sep = "\t", row.names = F, col.names = T, quote = F)

rm(list = ls())
gc();gc()
