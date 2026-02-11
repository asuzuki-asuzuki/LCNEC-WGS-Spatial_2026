library(Seurat) # v5.0.3
library(dplyr)
library(patchwork)

set.seed(1234)

packageVersion("Seurat")

# Load output of Cell Ranger
cell.data <- Read10X(data.dir = "<DIR>/filtered_feature_bc_matrix")

# Creat a Seurat object
cell <- CreateSeuratObject(counts = cell.data, project = "LCNEC", min.cells = 0, min.features = 0)

# Filtering
cell[["percent.mt"]] <- PercentageFeatureSet(cell, pattern = "^MT-")
cell <- subset(cell, subset = nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 10)

length(cell@meta.data$nFeature_RNA)

pdf("VlnPlot_QC_subset.pdf", width = 10, height = 5)
VlnPlot(cell, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Normalization, dimensional reduction, clustering and UMAP
cell <- NormalizeData(cell, normalization.method = "LogNormalize", scale.factor = 10000)
cell <- FindVariableFeatures(cell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cell)
cell <- ScaleData(cell, features = all.genes)
cell <- RunPCA(cell, verbose = FALSE)
cell <- FindNeighbors(cell, dims = 1:30, verbose = FALSE)
cell <- FindClusters(cell, verbose = FALSE)
cell <- RunUMAP(cell, dims = 1:30, verbose = FALSE)

# UMAP plot
pdf("UMAP_clsuter.pdf", width = 7, height = 6)
DimPlot(cell, reduction = "umap", label = TRUE)
dev.off()

# Save the object
saveRDS(cell, file = "scRNA.rds")

# DEG analysis
cell.markers <- FindAllMarkers(cell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "top10.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)
write.table(cell.markers, "all.csv", append = FALSE, quote = TRUE, sep = ",", row.names = TRUE, col.names = NA)

# Heatmap with DEGs
pdf("DoHeatmap_top10.pdf", width = 30, height = 20)
DoHeatmap(cell, features = top10$gene) + NoLegend()
dev.off()

rm(list = ls())
gc();gc()
