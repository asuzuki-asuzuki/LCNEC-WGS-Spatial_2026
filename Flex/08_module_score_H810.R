# =============================================================================
# 08_module_score_H810.R
# Score H810-derived pseudotime gene modules on the tissue UMAP
# (run after a tissue Seurat object `merged` has been loaded; see config below)
# =============================================================================

library(Seurat)
library(ggplot2)
library(RColorBrewer)

# =====================================================================
# Config (set before running)
# =====================================================================
# `merged`     : a pre-loaded Seurat object (tissue data, e.g. Xenium 5k)
#                containing the reduction "umap_harmony_1.5".
# `output_dir` : directory to write the PDF into.
# MODULE_TABLE : gene_cluster_table.csv produced by 07_pseudotime_modules.R
#                for the H810 cell line.

# output_dir   <- "08_module_score_result"
# merged       <- readRDS("path/to/tissue_merged.rds")
MODULE_TABLE <- "07_pseudotime_modules_result/siRNA_1/group_1/lineage1/gene_cluster_table.csv"

pdf_out <- file.path(output_dir, "sketch_harmony1.5_H810_ModuleScore.pdf")

# =====================================================================
# Load cluster gene list
# =====================================================================
DATA <- read.csv(MODULE_TABLE, header = TRUE)

cluster_ids <- sort(unique(DATA$cluster))
message("Clusters found: ", paste(cluster_ids, collapse = ", "))

gene_lists <- lapply(cluster_ids, function(cl) {
  DATA$gene[DATA$cluster == cl]
})
names(gene_lists) <- paste0("Cluster", cluster_ids)

for (nm in names(gene_lists)) {
  message(nm, ": ", length(gene_lists[[nm]]), " genes")
}

# =====================================================================
# Keep only genes present in the Seurat object
# =====================================================================
features_in_object <- rownames(merged)

gene_lists_used <- lapply(names(gene_lists), function(nm) {
  used <- gene_lists[[nm]][gene_lists[[nm]] %in% features_in_object]
  message(nm, ": ", length(gene_lists[[nm]]), " -> ", length(used), " genes used")
  used
})
names(gene_lists_used) <- names(gene_lists)

# Drop clusters with no remaining genes
gene_lists_used <- gene_lists_used[lengths(gene_lists_used) > 0]

if (length(gene_lists_used) == 0) {
  stop("No genes from the H810 gene lists were found in the Seurat object.")
}

# =====================================================================
# AddModuleScore
# =====================================================================
prefix <- "H810_"

for (nm in names(gene_lists_used)) {

  score_name <- paste0(prefix, nm)

  # Remove the score column if it already exists, then recompute
  existing_cols <- grep(
    paste0("^", score_name, "$"),
    colnames(merged@meta.data),
    value = TRUE
  )

  if (length(existing_cols) > 0) {
    merged@meta.data[, existing_cols] <- NULL
    message("Removed existing columns: ", paste(existing_cols, collapse = ", "))
  }

  merged <- AddModuleScore(
    object   = merged,
    features = list(gene_lists_used[[nm]]),
    name     = score_name
  )

  # AddModuleScore appends "1" (e.g. H810_Cluster11); rename to H810_Cluster1
  colnames(merged@meta.data) <- gsub(
    paste0(score_name, "1$"),
    score_name,
    colnames(merged@meta.data)
  )

  message("AddModuleScore done: ", score_name)
}

# =====================================================================
# FeaturePlot (on harmony_1.5 UMAP)
# =====================================================================
umap_cells <- rownames(Embeddings(merged, "umap_harmony_1.5"))

spectral_11 <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

score_names <- paste0(prefix, names(gene_lists_used))

# Existence check
missing_scores <- setdiff(score_names, colnames(merged@meta.data))
if (length(missing_scores) > 0) {
  stop("Missing score columns: ", paste(missing_scores, collapse = ", "))
}

pdf(pdf_out, width = 5 * length(score_names), height = 4)

print(
  FeaturePlot(
    merged,
    features  = score_names,
    reduction = "umap_harmony_1.5",
    cells     = umap_cells,
    order     = TRUE,
    ncol      = length(score_names)
  ) &
    scale_color_gradientn(
      colours = spectral_11
    )
)

dev.off()

message("PDF saved: ", pdf_out)

# =====================================================================
# Check
# =====================================================================
message("H810 module score columns:")
print(grep("^H810_Cluster", colnames(merged@meta.data), value = TRUE))

message("Duplicated metadata column names:")
print(colnames(merged@meta.data)[duplicated(colnames(merged@meta.data))])
