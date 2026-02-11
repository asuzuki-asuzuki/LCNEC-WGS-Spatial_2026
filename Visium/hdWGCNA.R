library(Seurat) # v5.0.3
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA) # v0.3.3

theme_set(theme_cowplot())
set.seed(1234)
enableWGCNAThreads(nThreads = 8)

packageVersion("Seurat")
packageVersion("hdWGCNA")

# Load a Seurat object
lung <- readRDS("lung_sub_merge_harmony.rds")

# Spatial assay
DefaultAssay(lung) <- "Spatial"

# Normlaization, dimentional reduction, clustering and UMAP
lung <- NormalizeData(lung)
lung <- FindVariableFeatures(lung)
lung <- ScaleData(lung)
lung <- RunPCA(lung)
lung <- FindNeighbors(lung, dims = 1:30)
lung <- FindClusters(lung, verbose = TRUE)
lung <- RunUMAP(lung, dims = 1:30)

# Processiong image data
image_df <- do.call(rbind, lapply(names(lung@images), function(x){
  lung@images[[x]]@coordinates
}))

new_meta <- merge(lung@meta.data, image_df, by='row.names')

rownames(new_meta) <- new_meta$Row.names
ix <- match(as.character(colnames(lung)), as.character(rownames(new_meta)))
new_meta <- new_meta[ix,]

lung@meta.data <- new_meta
head(image_df)

# Join layers
lung[["Spatial"]] <- JoinLayers(lung[["Spatial"]])

# Setup for hdWGCNA
lung <- SetupForWGCNA(lung, gene_select = "fraction", fraction = 0.05, wgcna_name = "vis")
lung <- MetaspotsByGroups(lung, group.by = c("original"), ident.group = "original", assay = 'Spatial', min_spots = 25)
lung <- NormalizeMetacells(lung)

m_obj <- GetMetacellObject(lung)
m_obj

lung <- SetDatExpr(lung, group.by=NULL, group_name = NULL)

# Soft-thresholding powers
lung <- TestSoftPowers(lung)
plot_list <- PlotSoftPowers(lung)

pdf("hdwgcna_softpower.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

# Network construction
lung <- ConstructNetwork(lung, tom_name='hdwgcna', overwrite_tom=TRUE, minModuleSize = 30)

pdf("hdwgcna_dendrogram.pdf")
PlotDendrogram(lung, main='Spatial hdWGCNA dendrogram')
dev.off()

# Calculation of module eigengenes (MEs)
lung <- ModuleEigengenes(lung)
lung <- ModuleConnectivity(lung)
lung <- ResetModuleNames(lung, new_name = "Module")

# Set MEs to metadata
MEs <- GetMEs(lung)
modules <- GetModules(lung)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
lung@meta.data <- cbind(lung@meta.data, MEs)

pdf("hdwgcna_dotplot.pdf")
p <- DotPlot(lung, features=mods, group.by = 'original', dot.min=0.1)
p <- p + coord_flip() + RotatedAxis() + scale_color_gradient2(high='red', mid='grey95', low='blue') + xlab('') + ylab('')
p
dev.off()

# Extract gene members and kMEs in each modules
genes <- lung@misc$vis$wgcna_modules
write.table(genes, "module_genes_kME.txt",append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

# Save the object
saveRDS(lung, file = "lung_hdwgcna.rds")

rm(list = ls())
gc();gc()
