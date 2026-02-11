library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(monocle3)
library(magrittr)
library(slingshot) # v2.10.0
library(viridis)
library(fields)

set.seed(1234)

packageVersion("Seurat")
packageVersion("slingshot")

# Load a object
sce <- readRDS("lung_slingshot_1_harmony_tumor.rds")

# Calculating average pseudotime if a cell was assigned to multiple trajectory paths.
combined_pseudo <- rowMeans(as.matrix(slingPseudotime(sce)), na.rm = TRUE)

# Plot with pseudotime
col <- viridis(100, option = "D")
pdf("slingshot_tumor_pseudotime.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 7), xpd = TRUE)
plot(reducedDims(sce)$UMAP.HARMONY, col = col[cut(combined_pseudo, breaks=100)], cex = 0.4, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
image.plot(legend.only = TRUE, zlim = range(combined_pseudo, na.rm = TRUE), col = col, legend.lab = "Pseudotime", legend.line = 3)
dev.off()

# Plot with NEUROD1 expression
gene_expr <- assays(sce)$logcounts["NEUROD1", ]
upper_limit <- quantile(gene_expr, 0.99, na.rm = TRUE)
gene_expr[gene_expr > upper_limit] <- upper_limit
idx <- order(gene_expr)
col_palette <- viridis(100, option = "A")

pdf("slingshot_tumor_NEUROD1.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 7), xpd = TRUE)
plot(reducedDims(sce)$UMAP.HARMONY[idx, ], col = col_palette[cut(gene_expr[idx], breaks=100)], cex = 0.4, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
image.plot(legend.only = TRUE, zlim = range(gene_expr, na.rm = TRUE), col = col_palette, legend.lab = "Expression", legend.line = 3)
dev.off()

# Plot with POU2F3 expression
gene_expr <- assays(sce)$logcounts["POU2F3", ]
upper_limit <- quantile(gene_expr, 0.99, na.rm = TRUE)
gene_expr[gene_expr > upper_limit] <- upper_limit
idx <- order(gene_expr)
col_palette <- viridis(100, option = "A")

pdf("slingshot_tumor_POU2F3.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 7), xpd = TRUE)
plot(reducedDims(sce)$UMAP.HARMONY[idx, ], col = col_palette[cut(gene_expr[idx], breaks=100)], cex = 0.4, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
image.plot(legend.only = TRUE, zlim = range(gene_expr, na.rm = TRUE), col = col_palette, legend.lab = "Expression", legend.line = 3)
dev.off()

# Plot with ASCL1 expression
gene_expr <- assays(sce)$logcounts["ASCL1", ]
upper_limit <- quantile(gene_expr, 0.99, na.rm = TRUE)
gene_expr[gene_expr > upper_limit] <- upper_limit
idx <- order(gene_expr)
col_palette <- viridis(100, option = "A")

pdf("slingshot_tumor_ASCL1.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 7), xpd = TRUE)
plot(reducedDims(sce)$UMAP.HARMONY[idx, ], col = col_palette[cut(gene_expr[idx], breaks=100)], cex = 0.4, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
image.plot(legend.only = TRUE, zlim = range(gene_expr, na.rm = TRUE), col = col_palette, legend.lab = "Expression", legend.line = 3)
dev.off()

# Plot with YAP1 expression
gene_expr <- assays(sce)$logcounts["YAP1", ]
upper_limit <- quantile(gene_expr, 0.99, na.rm = TRUE)
gene_expr[gene_expr > upper_limit] <- upper_limit
idx <- order(gene_expr)
col_palette <- viridis(100, option = "A")

pdf("slingshot_tumor_YAP1.pdf", width = 7, height = 7)
par(mar = c(5, 4, 4, 7), xpd = TRUE)
plot(reducedDims(sce)$UMAP.HARMONY[idx, ], col = col_palette[cut(gene_expr[idx], breaks=100)], cex = 0.4, pch = 16)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
image.plot(legend.only = TRUE, zlim = range(gene_expr, na.rm = TRUE), col = col_palette, legend.lab = "Expression", legend.line = 3)
dev.off()

rm(list = ls())
gc();gc()
