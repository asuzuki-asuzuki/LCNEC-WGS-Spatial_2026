library(Seurat) # v5.0.3
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA) # v0.3.3
library(enrichR)
library(GeneOverlap)

theme_set(theme_cowplot())
set.seed(1234)

packageVersion("Seurat")
packageVersion("hdWGCNA")

# Load a object
lung <- readRDS("lung_hdwgcna.rds")
lung

# GO analysis
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
lung <- RunEnrichr(lung, dbs=dbs, max_genes = 50)
enrich_df <- GetEnrichrTable(lung)
write.table(enrich_df, "enrichment_GO_2023.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)

# Barplots
EnrichrBarPlot(lung, outdir = "enrichr_plots", n_terms = 10, plot_size = c(10,7), logscale=TRUE)

# Dotplots
pdf("dotplot_enrichment_GO_P.pdf", width = 12, height = 12)
EnrichrDotPlot(lung, mods = "all", database = "GO_Biological_Process_2023", n_terms=2, term_size=6, p_adj = FALSE) + scale_color_stepsn(colors=rev(viridis::magma(256)))
dev.off()

rm(list = ls())
gc();gc()
