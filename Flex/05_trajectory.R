# =============================================================================
# 05_trajectory.R
# Slingshot trajectory analysis
# - Start cluster set per sample
# - Works with or without subset (runs all patterns when use_subset = TRUE)
# - Run space (UMAP / PCA) is switchable
# Output: 05_trajectory_result/siRNA_2/group_NE/ etc.
# =============================================================================

library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(viridis)
library(fields)
set.seed(1234)

# =====================================================================
# User settings
# =====================================================================

# --- Per-sample config ---
# start_cluster : Slingshot start cluster (string)
# use_subset    : TRUE -> use 04_subset_result, FALSE -> use 02_filter_umap_result
SAMPLE_CONFIG <- list(
  siRNA_1   = list(start_cluster = "0",  use_subset = TRUE)
)

# --- Samples to run ---
RUN_SAMPLES <- c("siRNA_1")

# --- Slingshot run space ---
# "umap" -> run in UMAP space (also draws trajectory lines)
# "pca"  -> run in PCA space (pseudotime only, visualized on UMAP)
REDUCTION <- "umap"

# --- Cell-cycle correction (selects the input RDS) ---
CC_CORRECTION <- TRUE

# --- Genes to plot ---
GENES_TO_PLOT <- c("ASCL1", "NEUROD1", "CALCA", "CHGA", "INSM1", "NKX2-1",
                   "YAP1", "COL1A1", "CA9", "HEY1", "NRARP", "FGL1",
                   "NOTCH1", "DLL3")

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

INPUT_BASE_FULL   <- "02_filter_umap_result"
INPUT_BASE_SUBSET <- "04_subset_result"
OUTPUT_BASE       <- "05_trajectory_result"

# =====================================================================
# Helper functions
# =====================================================================

# safe_cut: keep cut() working even when the range is zero
safe_cut <- function(x, n = 100) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) rng <- rng + c(-0.01, 0.01)
  cut(x, breaks = seq(rng[1], rng[2], length.out = n + 1),
      include.lowest = TRUE)
}

# run_slingshot_and_plot: run Slingshot + plots + save RDS
run_slingshot_and_plot <- function(cell, start_cluster, reduction,
                                   out_dir, label) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -------------------------------------------------------------------
  # Convert to SCE & add reductions
  # -------------------------------------------------------------------
  sce <- as.SingleCellExperiment(cell, assay = "RNA")
  reducedDim(sce, "UMAP.PCA") <- Embeddings(cell, "umap.pca")
  reducedDim(sce, "PCA")      <- Embeddings(cell, "pca")

  # Add cluster labels to SCE
  sce$clusters <- as.character(Idents(cell))

  # Check the start cluster exists
  if (!start_cluster %in% unique(sce$clusters)) {
    cat(sprintf("    [!] Start cluster '%s' not found. Available: %s\n",
                start_cluster,
                paste(sort(unique(sce$clusters)), collapse = ", ")))
    return(invisible(NULL))
  }

  # -------------------------------------------------------------------
  # Run Slingshot
  # -------------------------------------------------------------------
  redDim_name <- if (reduction == "umap") "UMAP.PCA" else "PCA"
  cat(sprintf("    Running Slingshot (space: %s, start: %s)...\n",
              redDim_name, start_cluster))

  sce <- slingshot(sce,
                   clusterLabels = "clusters",
                   reducedDim    = redDim_name,
                   start.clus    = start_cluster)

  n_lineages <- length(slingLineages(sce))
  cat(sprintf("    Number of lineages: %d\n", n_lineages))
  print(slingLineages(sce))

  # Sanity-check the start cluster pseudotime
  pseudo_mat  <- as.matrix(slingPseudotime(sce))
  cl_idx      <- which(sce$clusters == start_cluster)
  mean_pt_start <- mean(pseudo_mat[cl_idx, ], na.rm = TRUE)
  max_pt        <- max(pseudo_mat, na.rm = TRUE)
  cat(sprintf("    Start cluster '%s' mean pseudotime: %.3f / max: %.3f %s\n",
              start_cluster, mean_pt_start, max_pt,
              if (mean_pt_start < max_pt * 0.2) "OK" else "[!] check start cluster"))

  # -------------------------------------------------------------------
  # Coordinates & colors for plotting
  # -------------------------------------------------------------------
  umap_coords <- reducedDims(sce)$UMAP.PCA
  clus        <- factor(sce$clusters,
                        levels = sort(unique(sce$clusters)))
  n_clus      <- nlevels(clus)
  clus_colors <- setNames(hcl.colors(n_clus, palette = "Set 2"), levels(clus))

  col_pseudo <- viridis(100, option = "D")
  col_expr   <- viridis(100, option = "A")

  # Draw trajectory lines only in UMAP space
  draw_lines <- (reduction == "umap")

  # -------------------------------------------------------------------
  # Plot 1: by cluster
  # -------------------------------------------------------------------
  pdf(file.path(out_dir, "slingshot_clusters.pdf"), width = 7, height = 6)
  plot(umap_coords, col = clus_colors[clus], cex = 0.5, pch = 16,
       main = sprintf("%s - Clusters (start: %s)", label, start_cluster),
       xlab = "UMAP_1", ylab = "UMAP_2")
  for (cl in levels(clus)) {
    cx <- mean(umap_coords[clus == cl, 1])
    cy <- mean(umap_coords[clus == cl, 2])
    text(cx, cy, labels = cl, font = 2, cex = 0.8)
  }
  if (draw_lines) lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  dev.off()

  # -------------------------------------------------------------------
  # Plot 2: by sample
  # -------------------------------------------------------------------
  samples     <- sce$orig.ident
  present_lv  <- intersect(TIME_ORDER, unique(as.character(samples)))
  samples     <- factor(samples, levels = present_lv)
  sample_col  <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% present_lv]

  pdf(file.path(out_dir, "slingshot_samples.pdf"), width = 7, height = 6)
  plot(umap_coords, col = sample_col[as.character(samples)],
       cex = 0.5, pch = 16,
       main = sprintf("%s - Samples (start: %s)", label, start_cluster),
       xlab = "UMAP_1", ylab = "UMAP_2")
  if (draw_lines) lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  legend("topright", legend = names(sample_col), fill = sample_col,
         cex = 0.6, bty = "n")
  dev.off()

  # -------------------------------------------------------------------
  # Plot 3: pseudotime (mean over lineages)
  # -------------------------------------------------------------------
  combined_pseudo <- rowMeans(pseudo_mat, na.rm = TRUE)

  pdf(file.path(out_dir, "pseudotime.pdf"), width = 7, height = 6)
  par(mar = c(5, 4, 4, 7), xpd = TRUE)
  plot(umap_coords,
       col  = col_pseudo[safe_cut(combined_pseudo)],
       cex  = 0.4, pch = 16,
       main = sprintf("%s - Pseudotime (start: %s)", label, start_cluster),
       xlab = "UMAP_1", ylab = "UMAP_2")
  if (draw_lines) lines(SlingshotDataSet(sce), lwd = 2, col = "black")
  image.plot(legend.only = TRUE,
             zlim        = range(combined_pseudo, na.rm = TRUE),
             col         = col_pseudo,
             legend.lab  = "Pseudotime",
             legend.line = 3)
  dev.off()

  # -------------------------------------------------------------------
  # Plot 4: pseudotime per lineage
  # -------------------------------------------------------------------
  for (i in seq_len(ncol(pseudo_mat))) {
    pt_i  <- pseudo_mat[, i]
    valid <- !is.na(pt_i)
    pdf(file.path(out_dir, sprintf("pseudotime_lineage%d.pdf", i)),
        width = 7, height = 6)
    par(mar = c(5, 4, 4, 7), xpd = TRUE)
    plot(umap_coords, col = "grey85", cex = 0.3, pch = 16,
         main = sprintf("Lineage %d (start: %s)", i, start_cluster),
         xlab = "UMAP_1", ylab = "UMAP_2")
    points(umap_coords[valid, ],
           col = col_pseudo[safe_cut(pt_i[valid])],
           cex = 0.4, pch = 16)
    if (draw_lines) lines(SlingshotDataSet(sce), lwd = 2, col = "black")
    image.plot(legend.only = TRUE,
               zlim        = range(pt_i, na.rm = TRUE),
               col         = col_pseudo,
               legend.lab  = "Pseudotime",
               legend.line = 3)
    dev.off()
  }

  # -------------------------------------------------------------------
  # Plot 5: gene expression
  # -------------------------------------------------------------------
  gene_dir <- file.path(out_dir, "gene_expression")
  dir.create(gene_dir, showWarnings = FALSE)

  genes_ok <- GENES_TO_PLOT[GENES_TO_PLOT %in% rownames(sce)]
  missing  <- setdiff(GENES_TO_PLOT, genes_ok)
  if (length(missing) > 0)
    cat(sprintf("    [!] Genes not found: %s\n", paste(missing, collapse = ", ")))

  for (gene in genes_ok) {
    gene_expr <- assay(sce, "logcounts")[gene, ]
    upper_lim <- quantile(gene_expr, 0.99, na.rm = TRUE)
    gene_expr[gene_expr > upper_lim] <- upper_lim
    z_range <- range(gene_expr, na.rm = TRUE)
    if (diff(z_range) == 0) z_range <- z_range + c(-0.01, 0.01)
    idx <- order(gene_expr)

    pdf(file.path(gene_dir, sprintf("%s.pdf", gene)), width = 7, height = 6)
    par(mar = c(5, 4, 4, 7), xpd = TRUE)
    plot(umap_coords[idx, ],
         col  = col_expr[cut(gene_expr[idx], breaks = 100)],
         cex  = 0.4, pch = 16, main = gene,
         xlab = "UMAP_1", ylab = "UMAP_2")
    if (draw_lines) lines(SlingshotDataSet(sce), lwd = 2, col = "black")
    image.plot(legend.only = TRUE,
               zlim        = z_range,
               col         = col_expr,
               legend.lab  = "logcounts",
               legend.line = 3)
    dev.off()
  }
  cat(sprintf("    Gene expression plots done: %d genes\n", length(genes_ok)))

  # -------------------------------------------------------------------
  # Save SCE object
  # -------------------------------------------------------------------
  saveRDS(sce, file.path(out_dir, "sce_slingshot.rds"))
  cat(sprintf("    Saved: %s/\n", out_dir))
}

# =====================================================================
# Main
# =====================================================================

dir_suffix <- if (CC_CORRECTION) "_cc_corrected" else ""

for (sample_name in RUN_SAMPLES) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  config        <- SAMPLE_CONFIG[[sample_name]]
  start_cluster <- config$start_cluster
  use_subset    <- config$use_subset

  if (use_subset) {
    # ---------------------------------------------------------------
    # With subset: run all patterns under 04_subset_result
    # ---------------------------------------------------------------
    subset_base <- file.path(INPUT_BASE_SUBSET, sample_name)

    if (!dir.exists(subset_base)) {
      cat(sprintf("  Skip: %s not found\n", subset_base))
      next
    }

    subset_patterns <- list.dirs(subset_base, full.names = FALSE,
                                 recursive = FALSE)
    cat(sprintf("  Subset patterns: %s\n",
                paste(subset_patterns, collapse = ", ")))

    for (pattern_name in subset_patterns) {
      cat(sprintf("\n  -- Pattern: %s --\n", pattern_name))

      input_rds <- file.path(subset_base, pattern_name, "seurat_subset.rds")
      if (!file.exists(input_rds)) {
        cat(sprintf("    Skip: %s not found\n", input_rds))
        next
      }

      cell    <- readRDS(input_rds)
      out_dir <- file.path(OUTPUT_BASE, sample_name, pattern_name)
      label   <- sprintf("%s/%s", sample_name, pattern_name)

      run_slingshot_and_plot(cell, start_cluster, REDUCTION, out_dir, label)
    }

  } else {
    # ---------------------------------------------------------------
    # Without subset: use the processed object from 02_filter_umap_result
    # ---------------------------------------------------------------
    input_rds <- file.path(INPUT_BASE_FULL,
                            paste0(sample_name, dir_suffix),
                            "seurat_processed.rds")

    if (!file.exists(input_rds)) {
      cat(sprintf("  Skip: %s not found\n", input_rds))
      next
    }

    cell    <- readRDS(input_rds)
    out_dir <- file.path(OUTPUT_BASE, sample_name, "full")
    label   <- sprintf("%s/full", sample_name)

    run_slingshot_and_plot(cell, start_cluster, REDUCTION, out_dir, label)
  }
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
