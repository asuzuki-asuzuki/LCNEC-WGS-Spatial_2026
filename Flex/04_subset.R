# =============================================================================
# 04_subset.R
# Subset by selected clusters -> recompute UMAP -> save RDS
# Multiple (possibly overlapping) subset patterns can be specified
# Output: 04_subset_result/siRNA_2/subset_name/
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
set.seed(1234)

# =====================================================================
# User settings
# =====================================================================

# --- Target samples ---
SAMPLES <- c("siRNA_1")

# --- Input source ---
CC_CORRECTION <- TRUE   # choose CC-corrected / non-corrected from 02_filter_umap_result
INPUT_BASE    <- "02_filter_umap_result"

# --- Subset patterns (overlap allowed, multiple entries) ---
# Give each a name and the cluster IDs (as strings) to keep
SUBSET_PATTERNS <- list(
#  all_clusters   = NULL,          # NULL = all clusters (no subsetting)
  group_1       = c("0", "1", "2", "3", "4")
)

# --- Recompute UMAP? ---
RECOMPUTE <- TRUE   # TRUE: recompute / FALSE: reuse original UMAP coordinates

# --- Preprocessing parameters (used only when RECOMPUTE = TRUE) ---
N_FEATURES <- 2000
N_PCS      <- 20
RESOLUTION <- 0.2

# --- Time-point order ---
TIME_ORDER <- c(
  "0h", "Cont_24h", "Cont_48h", "Cont_72h", "Cont_96h", "Cont_168h",
  "siRNA_1_24h", "siRNA_1_48h", "siRNA_1_72h", "siRNA_1_96h", "siRNA_1_168h",
  "siRNA_2_24h", "siRNA_2_48h", "siRNA_2_72h", "siRNA_2_96h", "siRNA_2_168h",
  "siRNA_3_24h", "siRNA_3_48h", "siRNA_3_72h", "siRNA_3_96h", "siRNA_3_168h"
)

# --- Sample colors ---
red_pal   <- colorRampPalette(c("#FFAAAA", "#8B0000"))(5)
blue_pal  <- colorRampPalette(c("#AAAAFF", "#00008B"))(5)
green_pal <- colorRampPalette(c("#AAFFAA", "#006400"))(5)
cont_pal  <- colorRampPalette(c("#FFE066", "#8B6914"))(5)

SAMPLE_COLORS <- c(
  "0h"           = "#AAAAAA",
  "Cont_24h"     = cont_pal[1], "Cont_48h"     = cont_pal[2],
  "Cont_72h"     = cont_pal[3], "Cont_96h"     = cont_pal[4],
  "Cont_168h"    = cont_pal[5],
  "siRNA_1_24h"  = red_pal[1],  "siRNA_1_48h"  = red_pal[2],
  "siRNA_1_72h"  = red_pal[3],  "siRNA_1_96h"  = red_pal[4],
  "siRNA_1_168h" = red_pal[5],
  "siRNA_2_24h"  = blue_pal[1], "siRNA_2_48h"  = blue_pal[2],
  "siRNA_2_72h"  = blue_pal[3], "siRNA_2_96h"  = blue_pal[4],
  "siRNA_2_168h" = blue_pal[5],
  "siRNA_3_24h"  = green_pal[1],"siRNA_3_48h"  = green_pal[2],
  "siRNA_3_72h"  = green_pal[3],"siRNA_3_96h"  = green_pal[4],
  "siRNA_3_168h" = green_pal[5]
)

OUTPUT_BASE <- "04_subset_result"

# =====================================================================
# Main
# =====================================================================

for (sample_name in SAMPLES) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  # Input RDS
  dir_suffix <- if (CC_CORRECTION) "_cc_corrected" else ""
  input_rds  <- file.path(INPUT_BASE,
                           paste0(sample_name, dir_suffix),
                           "seurat_processed.rds")

  if (!file.exists(input_rds)) {
    cat(sprintf("  Skip: %s not found\n", input_rds))
    next
  }

  # Load full object
  cat("  Loading full object...\n")
  cell_full <- readRDS(input_rds)
  cat(sprintf("  %d cells, clusters: %s\n",
              ncol(cell_full),
              paste(levels(Idents(cell_full)), collapse = ", ")))

  # Order orig.ident by time point
  present_levels <- intersect(TIME_ORDER, unique(cell_full$orig.ident))
  cell_full$orig.ident <- factor(cell_full$orig.ident, levels = present_levels)
  present_colors <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% present_levels]

  # -------------------------------------------------------------------
  # Process each subset pattern
  # -------------------------------------------------------------------
  for (subset_name in names(SUBSET_PATTERNS)) {

    target_clusters <- SUBSET_PATTERNS[[subset_name]]

    cat(sprintf("\n  -- Pattern: %s --\n", subset_name))

    # Run subsetting
    if (is.null(target_clusters)) {
      cell <- cell_full
      cat("    Using all clusters\n")
    } else {
      # Existence check
      all_clusters <- levels(Idents(cell_full))
      missing_cl   <- setdiff(target_clusters, all_clusters)
      if (length(missing_cl) > 0) {
        cat(sprintf("    [!] Missing clusters: %s -> skip\n",
                    paste(missing_cl, collapse = ", ")))
        next
      }

      cell <- subset(cell_full, idents = target_clusters)
      cat(sprintf("    Extracted clusters %s: %d cells\n",
                  paste(target_clusters, collapse = ", "), ncol(cell)))
    }

    # Output directory
    out_dir <- file.path(OUTPUT_BASE, sample_name, subset_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # -----------------------------------------------------------------
    # Recompute UMAP or reuse original coordinates
    # -----------------------------------------------------------------
    if (RECOMPUTE) {
      cat("    Preprocessing & recomputing UMAP...\n")

      cell <- NormalizeData(cell, verbose = FALSE)
      cell <- FindVariableFeatures(cell, nfeatures = N_FEATURES, verbose = FALSE)

      # Match cell-cycle regression to CC_CORRECTION
      if (CC_CORRECTION) {
        cell <- CellCycleScoring(cell,
                                 s.features   = cc.genes.updated.2019$s.genes,
                                 g2m.features = cc.genes.updated.2019$g2m.genes,
                                 set.ident    = FALSE)
        cell <- ScaleData(cell,
                          vars.to.regress = c("S.Score", "G2M.Score"),
                          verbose         = FALSE)
      } else {
        cell <- ScaleData(cell, verbose = FALSE)
      }

      cell <- RunPCA(cell, npcs = N_PCS, verbose = FALSE)
      cell <- RunUMAP(cell, reduction = "pca", dims = 1:N_PCS,
                      reduction.name = "umap.pca", verbose = FALSE)
      cell <- FindNeighbors(cell, reduction = "pca", dims = 1:N_PCS,
                            verbose = FALSE)
      cell <- FindClusters(cell, resolution = RESOLUTION, verbose = FALSE)
      cat(sprintf("    Clusters after recompute: %d\n", nlevels(Idents(cell))))

    } else {
      cat("    Reusing original UMAP coordinates\n")
    }

    # -----------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------
    umap_name <- "umap.pca"

    # UMAP by cluster
    p_cluster <- DimPlot(cell,
                         reduction = umap_name,
                         label     = TRUE,
                         pt.size   = 0.4) +
      ggtitle(sprintf("%s - %s\nClusters", sample_name, subset_name)) +
      NoLegend()

    # UMAP by sample
    p_sample <- DimPlot(cell,
                        reduction = umap_name,
                        group.by  = "orig.ident",
                        cols      = present_colors,
                        pt.size   = 0.4) +
      ggtitle("Samples")

    pdf(file.path(out_dir, "umap_overview.pdf"), width = 12, height = 5)
    print(p_cluster | p_sample)
    dev.off()

    pdf(file.path(out_dir, "umap_clusters.pdf"), width = 6, height = 5)
    print(p_cluster)
    dev.off()

    pdf(file.path(out_dir, "umap_samples.pdf"), width = 7, height = 5)
    print(p_sample)
    dev.off()

    # -----------------------------------------------------------------
    # Save RDS
    # -----------------------------------------------------------------
    out_rds <- file.path(out_dir, "seurat_subset.rds")
    saveRDS(cell, out_rds)
    cat(sprintf("    Saved: %s\n", out_rds))
  }
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
