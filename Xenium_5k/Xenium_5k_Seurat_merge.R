library(Seurat) # v5.4.0
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)
library(SeuratObject)

set.seed(1234)
packageVersion("Seurat")

# Please read the sample list (Sample ID, Subtype)
DATA <- read.table("<LIST>", header=T)

# Load all Seurat objects
data <- list()
dir.create("bpcells_counts")

for (i in DATA$Case) {
    data[[i]] <- readRDS(paste0("../", i, "/lung_xenium_annotation.rds"))
    write_matrix_dir(mat = data[[i]][["Xenium"]]$counts, dir = paste0("bpcells_counts/", i))
    data[[i]][["Xenium"]]$counts <- open_matrix_dir(dir = paste0("bpcells_counts/", i))
    data[[i]]$original <- i
    data[[i]]$subtype <- DATA$Subtype[which(DATA$Case == i)]
    DefaultAssay(data[[i]]) <- "Xenium"
    message(i)
}

# Total of 13 objects


lung.merge <- merge(data[[1]], y = data[-1], add.cell.ids = names(data), project = "LCNEC")

# number of cells
length(lung.merge@meta.data$original)
table(lung.merge@meta.data$original)

DefaultAssay(lung.merge) <- "Xenium"
lung.merge

# Save the object
saveRDS(lung.merge, file = "lung_xenium_merge.rds")

rm(list = ls())
gc();gc()
