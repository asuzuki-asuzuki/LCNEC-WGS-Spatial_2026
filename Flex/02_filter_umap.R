# =============================================================================
# 02_filter_umap.R
# Filtering -> preprocessing -> UMAP
# - Filter by nFeature / percent.mt
# - Cell-cycle regression can be toggled on/off
# Output: 02_filter_umap_result/siRNA_1/ or siRNA_1_cc_corrected/
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
set.seed(1234)

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

# --- Filtering thresholds (shared across samples) ---
FILTER <- list(
  nFeature_min = 4000,    # minimum number of features
  nFeature_max = 11000,   # maximum number of features
  mt_max       = 15       # maximum mitochondrial gene percentage (%)
)

# --- Preprocessing parameters ---
N_FEATURES <- 2000
N_PCS      <- 30
RESOLUTION <- 0.5

# --- Cell-cycle regression ---
CC_CORRECTION <- TRUE   # TRUE: regress out / FALSE: no regression

# --- Time-point order (legend / color order on UMAP) ---
TIME_ORDER <- c(
  "0h", "Cont_24h", "Cont_48h", "Cont_72h", "Cont_96h", "Cont_168h",
  "siRNA_1_24h", "siRNA_1_48h", "siRNA_1_72h", "siRNA_1_96h", "siRNA_1_168h",
  "siRNA_2_24h", "siRNA_2_48h", "siRNA_2_72h", "siRNA_2_96h", "siRNA_2_168h",
  "siRNA_3_24h", "siRNA_3_48h", "siRNA_3_72h", "siRNA_3_96h", "siRNA_3_168h"
)

# --- Sample colors (0h: grey / siRNA_1: reds / siRNA_2: blues / siRNA_3: greens) ---
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

OUTPUT_BASE <- "02_filter_umap_result"

# =====================================================================
# Main
# =====================================================================

for (sample_name in names(SAMPLES)) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  # Output directory (suffix depends on cell-cycle correction)
  dir_suffix <- if (CC_CORRECTION) "_cc_corrected" else ""
  out_dir    <- file.path(OUTPUT_BASE, paste0(sample_name, dir_suffix))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------------------------------------------
  # 1. Load & filter
  # -----------------------------------------------------------------
  cat("  [1/6] Load & filter...\n")
  cell <- readRDS(SAMPLES[[sample_name]])
  DefaultAssay(cell) <- "RNA"
  cell[["percent.mt"]] <- PercentageFeatureSet(cell, pattern = "^MT-")

  n_before <- ncol(cell)
  keep <- rep(TRUE, ncol(cell))
  if (!is.null(FILTER$nFeature_min)) keep <- keep & cell$nFeature_RNA >= FILTER$nFeature_min
  if (!is.null(FILTER$nFeature_max)) keep <- keep & cell$nFeature_RNA <= FILTER$nFeature_max
  if (!is.null(FILTER$mt_max))       keep <- keep & cell$percent.mt   <= FILTER$mt_max
  cell <- cell[, keep]
  n_after <- ncol(cell)
  cat(sprintf("    %d -> %d cells (removed: %d)\n", n_before, n_after, n_before - n_after))

  # Order orig.ident by time point
  present_levels <- intersect(TIME_ORDER, unique(cell$orig.ident))
  cell$orig.ident <- factor(cell$orig.ident, levels = present_levels)

  # -----------------------------------------------------------------
  # 2. Preprocessing
  # -----------------------------------------------------------------
  cat("  [2/6] Preprocessing (Normalize -> VariableFeatures -> Scale -> PCA)...\n")
  cell <- NormalizeData(cell, verbose = FALSE)
  cell <- FindVariableFeatures(cell, nfeatures = N_FEATURES, verbose = FALSE)

  # Score cell cycle (computed regardless of regression setting)
  cell <- CellCycleScoring(cell,
                           s.features   = cc.genes.updated.2019$s.genes,
                           g2m.features = cc.genes.updated.2019$g2m.genes,
                           set.ident    = FALSE)
  phase_pct <- round(prop.table(table(cell$Phase)) * 100, 1)
  cat(sprintf("    Cell cycle: G1=%.1f%%  S=%.1f%%  G2M=%.1f%%\n",
              phase_pct["G1"], phase_pct["S"], phase_pct["G2M"]))

  # ScaleData (toggle cell-cycle regression)
  if (CC_CORRECTION) {
    cat("    ScaleData (regress out cell-cycle scores)...\n")
    cell <- ScaleData(cell,
                      vars.to.regress = c("S.Score", "G2M.Score"),
                      verbose         = FALSE)
  } else {
    cell <- ScaleData(cell, verbose = FALSE)
  }

  cell <- RunPCA(cell, npcs = N_PCS, verbose = FALSE)

  # -----------------------------------------------------------------
  # 3. UMAP & clustering
  # -----------------------------------------------------------------
  cat("  [3/6] UMAP & clustering...\n")
  cell <- RunUMAP(cell, reduction = "pca", dims = 1:N_PCS,
                  reduction.name = "umap.pca", verbose = FALSE)
  cell <- FindNeighbors(cell, reduction = "pca", dims = 1:N_PCS, verbose = FALSE)
  cell <- FindClusters(cell, resolution = RESOLUTION, verbose = FALSE)
  cat(sprintf("    Number of clusters: %d\n", nlevels(Idents(cell))))

  # -----------------------------------------------------------------
  # 4. Plots
  # -----------------------------------------------------------------
  cat("  [4/6] Plotting...\n")

  # --- UMAP by cluster ---
  p_cluster <- DimPlot(cell,
                       reduction = "umap.pca",
                       label     = TRUE,
                       pt.size   = 0.3) +
    ggtitle(paste(sample_name, "- Clusters")) +
    NoLegend()

  pdf(file.path(out_dir, "umap_clusters.pdf"), width = 6, height = 5)
  print(p_cluster)
  dev.off()

  # --- UMAP by sample ---
  present_colors <- SAMPLE_COLORS[names(SAMPLE_COLORS) %in% present_levels]

  p_sample <- DimPlot(cell,
                      reduction = "umap.pca",
                      group.by  = "orig.ident",
                      cols      = present_colors,
                      pt.size   = 0.3) +
    ggtitle(paste(sample_name, "- Samples"))

  pdf(file.path(out_dir, "umap_samples.pdf"), width = 8, height = 5)
  print(p_sample)
  dev.off()

  # --- Cluster + sample side by side ---
  pdf(file.path(out_dir, "umap_overview.pdf"), width = 14, height = 5)
  print(p_cluster | p_sample)
  dev.off()

  # --- Cell-cycle phase / score UMAP ---
  p_phase <- DimPlot(cell,
                     reduction = "umap.pca",
                     group.by  = "Phase",
                     cols      = c("G1" = "#AAAAAA", "S" = "#E69F00", "G2M" = "#56B4E9"),
                     pt.size   = 0.3) +
    ggtitle("Cell Cycle Phase")

  p_sscore <- FeaturePlot(cell,
                          features  = "S.Score",
                          reduction = "umap.pca",
                          order     = TRUE) +
    ggtitle("S.Score")

  p_g2mscore <- FeaturePlot(cell,
                            features  = "G2M.Score",
                            reduction = "umap.pca",
                            order     = TRUE) +
    ggtitle("G2M.Score")

  pdf(file.path(out_dir, "umap_cellcycle.pdf"), width = 16, height = 5)
  print(p_phase | p_sscore | p_g2mscore)
  dev.off()

  # -----------------------------------------------------------------
  # 5. DEG analysis
  # -----------------------------------------------------------------
  cat("  [5/6] DEG analysis (FindAllMarkers)...\n")

  markers <- FindAllMarkers(cell,
                            only.pos        = TRUE,
                            min.pct         = 0.25,
                            logfc.threshold = 0.25,
                            verbose         = FALSE)
  cat(sprintf("    Total DEGs: %d\n", nrow(markers)))

  top10 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))

  write.csv(markers,
            file      = file.path(out_dir, "deg_all.csv"),
            row.names = FALSE)

  write.csv(top10,
            file      = file.path(out_dir, "deg_top10.csv"),
            row.names = FALSE)

  cat(sprintf("    Saved: deg_all.csv / deg_top10.csv\n"))

  # -----------------------------------------------------------------
  # 6. Save processed object
  # -----------------------------------------------------------------
  cat("  [6/6] Save object...\n")
  out_rds <- file.path(out_dir, "seurat_processed.rds")
  saveRDS(cell, out_rds)
  cat(sprintf("    Saved: %s\n", out_rds))
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
