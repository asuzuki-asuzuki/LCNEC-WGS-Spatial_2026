library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(GenomicRanges)

set.seed(1234)
packageVersion("Seurat")

# Load Seurat objects
lung_t4 <- readRDS("<DIR4>/scRNA_annotation.rds")
lung_t6 <- readRDS("<DIR6>/scRNA_annotation.rds")
lung_t9 <- readRDS("<DIR9>/scRNA_annotation.rds")
lung_t11 <- readRDS("<DIR11>/scRNA_annotation.rds")
lung_j23 <- readRDS("<DIR23>/scRNA_annotation.rds")
lung_j28 <- readRDS("<DIR28>/scRNA_annotation.rds")
# Total of 6 objects

lung_t4$original <- "<SAMPLE4>"
lung_t6$original <- "<SAMPLE6>"
lung_t9$original <- "<SAMPLE9>"
lung_t11$original <- "<SAMPLE11>"
lung_j23$original <- "<SAMPLE23>"
lung_j28$original <- "<SAMPLE28>"
# Total of 6 samples

# Merge
lung.merge <- merge(lung_t4, y = c(lung_t6, lung_t9, lung_t11, lung_j23, lung_j28), add.cell.ids = c("lung_t4", "lung_t6", "lung_t9", "lung_t11", "lung_j23", "lung_j28"), project = "LCNEC")

# Number of cells
length(lung.merge@meta.data$original)
table(lung.merge@meta.data$original)

# Save the object
saveRDS(lung.merge, file = "lung_merge_annotation.rds")

rm(list = ls())
gc();gc()
