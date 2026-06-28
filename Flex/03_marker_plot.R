# =============================================================================
# 03_marker_plot.R
# FeaturePlot / VlnPlot of selected genes and module-score visualization
# Output: 03_marker_plot_result/
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(scales)
set.seed(1234)

# =====================================================================
# User settings
# =====================================================================

# --- Input RDS ---
# CC_CORRECTION = TRUE  -> 02_filter_umap_result/siRNA_1_cc_corrected/
# CC_CORRECTION = FALSE -> 02_filter_umap_result/siRNA_1/
CC_CORRECTION <- TRUE

SAMPLES <- c("siRNA_1", "siRNA_2", "siRNA_3", "siRNA_all", "control", "siRNA_all_cont")

# --- Genes to plot ---
GENES_TO_PLOT <- c("ASCL1", "NEUROD1", "CALCA", "CHGA", "INSM1", "NKX2-1",
                   "YAP1", "POU2F3", "ASCL2", "CD274", "IDO1", "FGL1",
                   "NOTCH1", "DLL3", "NRARP", "DLX2")

# --- Module scores ---
MODULE_FILE <- "./module_genes_kME.txt"   # NULL to skip
MODULES <- list(
  Module17 = list(col = "module", val = "Module17", kme_col = "kME_Module17"),
  Module13 = list(col = "module", val = "Module13", kme_col = "kME_Module13"),
  Module11 = list(col = "module", val = "Module11", kme_col = "kME_Module11")
)
KME_THRESHOLD <- 0.4

# --- Time-point order ---
TIME_ORDER <- c(
  "0h", "Cont_24h", "Cont_48h", "Cont_72h", "Cont_96h", "Cont_168h",
  "siRNA_1_24h", "siRNA_1_48h", "siRNA_1_72h", "siRNA_1_96h", "siRNA_1_168h",
  "siRNA_2_24h", "siRNA_2_48h", "siRNA_2_72h", "siRNA_2_96h", "siRNA_2_168h",
  "siRNA_3_24h", "siRNA_3_48h", "siRNA_3_72h", "siRNA_3_96h", "siRNA_3_168h"
)

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

OUTPUT_BASE <- "03_marker_plot_result"
INPUT_BASE  <- "02_filter_umap_result"

# =====================================================================
# Main
# =====================================================================

for (sample_name in SAMPLES) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  # Resolve input / output paths
  dir_suffix <- if (CC_CORRECTION) "_cc_corrected" else ""
  input_rds  <- file.path(INPUT_BASE,
                           paste0(sample_name, dir_suffix),
                           "seurat_processed.rds")

  if (!file.exists(input_rds)) {
    cat(sprintf("  Skip: %s not found\n", input_rds))
    next
  }

  out_dir <- file.path(OUTPUT_BASE, paste0(sample_name, dir_suffix))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------------------------------------------
  # 1. Load
  # -----------------------------------------------------------------
  cat("  [1/3] Loading...\n")
  cell <- readRDS(input_rds)
  cat(sprintf("    %d cells, %d clusters\n",
              ncol(cell), nlevels(Idents(cell))))

  # Order orig.ident by time point
  present_levels <- intersect(TIME_ORDER, unique(cell$orig.ident))
  cell$orig.ident <- factor(cell$orig.ident, levels = present_levels)

  # Drop genes not present in the data
  missing <- GENES_TO_PLOT[!GENES_TO_PLOT %in% rownames(cell)]
  if (length(missing) > 0)
    cat(sprintf("    [!] Genes not in data: %s\n",
                paste(missing, collapse = ", ")))
  genes_ok <- GENES_TO_PLOT[GENES_TO_PLOT %in% rownames(cell)]
  cat(sprintf("    Genes to plot: %d / %d\n",
              length(genes_ok), length(GENES_TO_PLOT)))

  # -----------------------------------------------------------------
  # 2. Gene FeaturePlot & VlnPlot
  # -----------------------------------------------------------------
  cat("  [2/3] Gene plots...\n")
  present_colors <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% present_levels]

  n_cols <- 4
  n_rows <- ceiling(length(genes_ok) / n_cols)

  # --- FeaturePlot: all genes (4 columns) ---
  pdf(file.path(out_dir, "featureplot_all_genes.pdf"),
      width  = n_cols * 4,
      height = n_rows * 4)
  print(
    FeaturePlot(cell,
                features  = genes_ok,
                reduction = "umap.pca",
                ncol      = n_cols,
                order     = TRUE)
  )
  dev.off()

  # --- VlnPlot: all genes (page 1 by cluster, page 2 by sample) ---
  pdf(file.path(out_dir, "vlnplot_genes.pdf"),
      width  = n_cols * 3,
      height = n_rows * 3)

  # Page 1: by cluster
  print(
    VlnPlot(cell,
            features = genes_ok,
            ncol     = n_cols,
            pt.size  = 0) &
      theme(legend.position = "none")
  )

  # Page 2: by sample
  print(
    VlnPlot(cell,
            features = genes_ok,
            group.by = "orig.ident",
	    cols     = present_colors,
            ncol     = n_cols,
            pt.size  = 0) &
      theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 6),
            legend.position  = "none")
  )
  dev.off()

  # --- Per-gene PDF (page 1 = FeaturePlot / page 2 = VlnPlots stacked) ---
  gene_dir <- file.path(out_dir, "individual_genes")
  dir.create(gene_dir, showWarnings = FALSE)

  for (gene in genes_ok) {
    pdf(file.path(gene_dir, sprintf("%s.pdf", gene)), width = 12, height = 10)

    # Page 1: FeaturePlot only
    print(
      FeaturePlot(cell,
                  features  = gene,
                  reduction = "umap.pca",
                  order     = TRUE) +
        ggtitle(paste(gene, "- FeaturePlot"))
    )

    # Page 2: VlnPlot by cluster / by sample, stacked
    p_vl_cl <- VlnPlot(cell,
                        features = gene,
                        pt.size  = 0) +
      ggtitle(paste(gene, "- by cluster")) +
      theme(legend.position = "none")

    p_vl_sm <- VlnPlot(cell,
                        features = gene,
                        group.by = "orig.ident",
			cols     = present_colors,
                        pt.size  = 0) +
      ggtitle(paste(gene, "- by sample")) +
      theme(axis.text.x     = element_text(angle = 45, hjust = 1),
            legend.position  = "none")

    print(p_vl_cl / p_vl_sm)
    dev.off()
    cat(sprintf("    Done: %s\n", gene))
  }

  # -----------------------------------------------------------------
  # 3. Module scores (skipped if MODULE_FILE is NULL)
  # -----------------------------------------------------------------
  if (is.null(MODULE_FILE)) {
    cat("  [3/3] Module scores: skipped\n")
  } else {
    cat("  [3/3] Module scores & plots...\n")

    DATA <- read.table(MODULE_FILE, header = TRUE)

    module_names <- names(MODULES)
    for (mod_name in module_names) {
      mod_info  <- MODULES[[mod_name]]
      mod_genes <- DATA$gene_name[
        DATA[[mod_info$col]] == mod_info$val &
        DATA[[mod_info$kme_col]] > KME_THRESHOLD
      ]
      mod_genes_ok <- mod_genes[mod_genes %in% rownames(cell)]
      cat(sprintf("    %s: %d genes (in data: %d)\n",
                  mod_name, length(mod_genes), length(mod_genes_ok)))

      cell <- AddModuleScore(cell,
                             features = list(mod_genes_ok),
                             name     = mod_name)
      # AddModuleScore appends "1" to the name; strip it
      colnames(cell@meta.data) <- gsub(paste0(mod_name, "1"),
                                       mod_name,
                                       colnames(cell@meta.data))
    }

    # Check score distributions
    for (mod_name in module_names) {
      s <- cell@meta.data[[mod_name]]
      cat(sprintf("    %s score: min=%.3f, median=%.3f, max=%.3f\n",
                  mod_name, min(s), median(s), max(s)))
    }

    # FeaturePlot: all modules
    pdf(file.path(out_dir, "module_score_featureplot.pdf"),
        width  = length(module_names) * 5,
        height = 4)
    print(
      FeaturePlot(cell,
                  features  = module_names,
                  reduction = "umap.pca",
                  ncol      = length(module_names),
                  order     = TRUE) &
        scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")),
                              limits  = c(-0.1, 0.2),
                              oob     = squish)
    )
    dev.off()

    # VlnPlot: page 1 by cluster, page 2 by sample
    pdf(file.path(out_dir, "module_score_vlnplot.pdf"),
        width  = length(module_names) * 4,
        height = 5)

    # Page 1: by cluster
    print(
      VlnPlot(cell,
              features = module_names,
              ncol     = length(module_names),
              pt.size  = 0) &
        theme(legend.position = "none")
    )

    # Page 2: by sample
    print(
      VlnPlot(cell,
              features = module_names,
              group.by = "orig.ident",
	      cols     = present_colors,
              ncol     = length(module_names),
              pt.size  = 0) &
        theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
              legend.position  = "none")
    )
    dev.off()

    cat("    Module-score plots done\n")
  }

  cat(sprintf("  Done: %s/\n", out_dir))
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
