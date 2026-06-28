# =============================================================================
# 01_qc.R
# Check QC metrics (nCount_RNA / nFeature_RNA / percent.mt)
# Output: 01_qc_result/
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# =====================================================================
# User settings
# =====================================================================

SAMPLES <- list(
  siRNA_1   = "../siRNA_1.rds",
  siRNA_2   = "../siRNA_2.rds",
  siRNA_3   = "../siRNA_3.rds",
  siRNA_all = "../siRNA_all.rds",
  control   = "../control.rds",
  siRNA_all_cont   = "../siRNA_all_control.rds"
)

TIME_ORDER <- c(
  "0h", "Cont_24h", "Cont_48h", "Cont_72h", "Cont_96h", "Cont_168h",
  "siRNA_1_24h", "siRNA_1_48h", "siRNA_1_72h", "siRNA_1_96h", "siRNA_1_168h",
  "siRNA_2_24h", "siRNA_2_48h", "siRNA_2_72h", "siRNA_2_96h", "siRNA_2_168h",
  "siRNA_3_24h", "siRNA_3_48h", "siRNA_3_72h", "siRNA_3_96h", "siRNA_3_168h"
)

OUTPUT_BASE <- "01_qc_result"

# =====================================================================
# Main
# =====================================================================

for (sample_name in names(SAMPLES)) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  out_dir <- file.path(OUTPUT_BASE, sample_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------------------------------------------
  # Load data
  # -----------------------------------------------------------------
  cell <- readRDS(SAMPLES[[sample_name]])
  DefaultAssay(cell) <- "RNA"
  cat(sprintf("  %d cells, %d features\n", ncol(cell), nrow(cell)))

  # Mitochondrial gene percentage
  cell[["percent.mt"]] <- PercentageFeatureSet(cell, pattern = "^MT-")

  # Order orig.ident by time point
  present_levels <- intersect(TIME_ORDER, unique(cell$orig.ident))
  cell$orig.ident <- factor(cell$orig.ident, levels = present_levels)

  # Dummy column to draw a single violin over all cells
  cell[["all_cells"]] <- sample_name

  # -----------------------------------------------------------------
  # Print QC summary to console
  # -----------------------------------------------------------------
  cat("  QC summary per sample:\n")
  for (s in present_levels) {
    idx   <- cell$orig.ident == s
    nfeat  <- cell$nFeature_RNA[idx]
    ncount <- cell$nCount_RNA[idx]
    mt     <- cell$percent.mt[idx]
    cat(sprintf("    %-20s  n=%4d  nFeature[med=%4.0f, max=%5.0f]  nCount[med=%5.0f]  MT[med=%.1f%%, max=%.1f%%]\n",
                s, sum(idx),
                median(nfeat), max(nfeat),
                median(ncount),
                median(mt), max(mt)))
  }

  # -----------------------------------------------------------------
  # Overview violins (all cells in one column)
  # -----------------------------------------------------------------
  p_count_all <- VlnPlot(cell,
                         features = "nCount_RNA",
                         group.by = "all_cells",
                         pt.size  = 0) +
    ggtitle("nCount_RNA") +
    xlab("") +
    theme(legend.position = "none")

  p_feature_all <- VlnPlot(cell,
                           features = "nFeature_RNA",
                           group.by = "all_cells",
                           pt.size  = 0) +
    ggtitle("nFeature_RNA") +
    xlab("") +
    theme(legend.position = "none")

  p_mt_all <- VlnPlot(cell,
                      features = "percent.mt",
                      group.by = "all_cells",
                      pt.size  = 0) +
    ggtitle("percent.mt (%)") +
    xlab("") +
    theme(legend.position = "none")

  # Three panels side by side
  pdf(file.path(out_dir, "qc_overview.pdf"), width = 10, height = 5)
  print(p_count_all | p_feature_all | p_mt_all)
  dev.off()

  # -----------------------------------------------------------------
  # Per-sample violins (separate PDFs)
  # -----------------------------------------------------------------
  n_samples <- length(present_levels)
  fig_width <- max(6, n_samples * 1.2)   # widen by number of samples

  p_count_by_sample <- VlnPlot(cell,
                               features = "nCount_RNA",
                               group.by = "orig.ident",
                               pt.size  = 0) +
    ggtitle(paste(sample_name, "- nCount_RNA by sample")) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p_feature_by_sample <- VlnPlot(cell,
                                 features = "nFeature_RNA",
                                 group.by = "orig.ident",
                                 pt.size  = 0) +
    ggtitle(paste(sample_name, "- nFeature_RNA by sample")) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  p_mt_by_sample <- VlnPlot(cell,
                             features = "percent.mt",
                             group.by = "orig.ident",
                             pt.size  = 0) +
    ggtitle(paste(sample_name, "- percent.mt (%) by sample")) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  pdf(file.path(out_dir, "qc_nCount_by_sample.pdf"),
      width = fig_width, height = 5)
  print(p_count_by_sample)
  dev.off()

  pdf(file.path(out_dir, "qc_nFeature_by_sample.pdf"),
      width = fig_width, height = 5)
  print(p_feature_by_sample)
  dev.off()

  pdf(file.path(out_dir, "qc_percent_mt_by_sample.pdf"),
      width = fig_width, height = 5)
  print(p_mt_by_sample)
  dev.off()

  cat(sprintf("  Saved: %s/\n", out_dir))
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
