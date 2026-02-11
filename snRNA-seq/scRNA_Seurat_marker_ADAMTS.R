library(Seurat) # v5.0.3
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)

set.seed(1234)

packageVersion("Seurat")

# Load a Seurat object
lung <- readRDS("lung_merge_annotation_harmony.rds")
lung

# UMAP plot with ADAMTS gene family expressions
pdf("marker_ADAMTS.pdf", width = 20, height = 5)
p1 <- FeaturePlot(lung, features = c("ADAMTS20", "ADAMTS9", "ADAMTS12", "ADAMTS2"), ncol=4, reduction="umap.harmony")
p1
dev.off()

# Violin plot with ADAMTS gene family expressions
mods <- c("ADAMTS20", "ADAMTS9", "ADAMTS12", "ADAMTS2")

# Set colors for 6 samples
lung$original <- factor(x = lung$original, levels = c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")) # please enter the sample names into the vector
COL <- c("#ff0000", "#ff0000", "#00b0f0", "#00b0f0", "#00b050", "#ffc000")
names(COL) <- c("<SAMPLE4>", "<SAMPLE6>", "<SAMPLE9>", "<SAMPLE11>", "<SAMPLE23>", "<SAMPLE28>")

# Join layers
lung[["RNA"]] <- JoinLayers(lung[["RNA"]])

# Subset tumor cells only
lung <- subset(lung, subset = annotation %in% c("Tumor cell"))

# Violin plots
pdf("ADAMTS_Vlnplot_ANPY.pdf", width = 20, height = 3)
p <- VlnPlot(lung, features = mods, ncol = 4, cols = COL, group.by = "original", pt.size = 0)
p & geom_jitter(aes(color = lung$original), size = 0.1, alpha = 0.5) & scale_color_manual(values = COL) & xlab(NULL)
dev.off()

# Violin plots for subtypes
A <- c("<SAMPLE9>", "<SAMPLE11>") # ASCL1 type cases
N <- c("<SAMPLE23>", "<SAMPLE4>") # NEUROD1 type cases
P <- c("<SAMPLE28>") # POU2F3 type cases
Y <- c("<SAMPLE6>") # YAP1 type cases

lung$subtype <- case_when(
 lung$original %in% A ~ "ASCL1",
 lung$original %in% N ~ "NEUROD1",
 lung$original %in% P ~ "POU2F3",
 lung$original %in% Y ~ "YAP1",
 TRUE ~ NA_character_
)

target_groups <- c("ASCL1", "NEUROD1", "POU2F3", "YAP1")
comparisons <- combn(target_groups, 2, simplify = FALSE)

COL2 <- c("#ff0000", "#00b0f0", "#00b050", "#ffc000")
names(COL2) <- target_groups

pdf("ADAMTS_Vlnplot_ANPY_p.pdf", width = 20, height = 5)
plots <- lapply(mods, function(gene) {
      actual_max <- max(FetchData(lung, vars = gene)[,1]) #max
      p <- VlnPlot(lung, features = gene, cols = COL2, group.by = "subtype", pt.size = 0, y.max = actual_max * 2.2)
      p + geom_jitter(aes(color = ident), size = 0.1, alpha = 0.5) + scale_color_manual(values = COL2) + xlab(NULL) + NoLegend() +
      stat_compare_means(comparisons = comparisons, label = "p.format", step.increase = 0.2, method = "wilcox.test")
})
CombinePlots(plots, ncol = 4)
dev.off()

rm(list = ls())
gc();gc()
