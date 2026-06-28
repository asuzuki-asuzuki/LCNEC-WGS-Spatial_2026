#!/usr/bin/env Rscript
# =============================================================================
# 07_pseudotime_modules.R
# Cluster pseudotime-associated genes (from tradeSeq) into modules
# Output: 07_pseudotime_modules_result/siRNA_1/group_1/lineage1/
# =============================================================================

library(tradeSeq)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(cluster)
library(patchwork)
set.seed(1234)

# =====================================================================
# User settings
# =====================================================================

SAMPLE_NAME         <- "siRNA_1"
SUBSET_PATTERN      <- "group_1"
TARGET_LINEAGE      <- 1

TOP_N_GENES         <- 1000
N_PSEUDOTIME_POINTS <- 10      # number of GAM sampling points
N_CLUSTERS_OVERRIDE <- 3       # NULL: auto by silhouette / number: manual
K_RANGE             <- 3:5
N_REPR_GENES        <- 12
PADJ_THRESHOLD      <- 0.05
EXPR_RANGE_MIN      <- 0.5

INPUT_BASE  <- "06_tradeseq_result"
OUTPUT_BASE <- "07_pseudotime_modules_result"

# =====================================================================
# Main
# =====================================================================

out_dir <- file.path(OUTPUT_BASE, SAMPLE_NAME, SUBSET_PATTERN,
                     sprintf("lineage%d", TARGET_LINEAGE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------
cat("[1/6] Loading data...\n")

sce_trade <- readRDS(file.path(INPUT_BASE, SAMPLE_NAME, SUBSET_PATTERN,
                               "sce_tradeseq.rds"))
sce_sling  <- readRDS(file.path("05_trajectory_result", SAMPLE_NAME,
                                 SUBSET_PATTERN, "sce_slingshot.rds"))
asso_res  <- read.csv(file.path(INPUT_BASE, SAMPLE_NAME, SUBSET_PATTERN,
                                sprintf("lineage%d_association.csv",
                                        TARGET_LINEAGE)))

# Get pseudotime and logcounts
pseudo_mat   <- as.matrix(slingPseudotime(sce_sling))
cell_pt      <- pseudo_mat[, TARGET_LINEAGE]      # pseudotime per cell
on_lineage   <- !is.na(cell_pt)                   # cells on this lineage
logcounts_mat <- assay(sce_sling, "logcounts")

cat(sprintf("  %d cells (on lineage%d: %d cells)\n",
            ncol(sce_sling), TARGET_LINEAGE, sum(on_lineage)))
rm(sce_sling); gc()

# -----------------------------------------------------------------
# 2. Select top genes
# -----------------------------------------------------------------
cat(sprintf("[2/6] Selecting top %d genes...\n", TOP_N_GENES))

top_genes <- asso_res %>%
  filter(padj < PADJ_THRESHOLD) %>%
  filter(expr_range >= EXPR_RANGE_MIN) %>%
  arrange(padj, desc(waldStat)) %>%
  head(TOP_N_GENES) %>%
  pull(gene)

# Keep only genes present in sce_tradeseq
top_genes <- top_genes[top_genes %in% rownames(sce_trade)]

cat(sprintf("  Significant genes: %d -> used: %d\n",
            sum(asso_res$padj < PADJ_THRESHOLD, na.rm = TRUE),
            length(top_genes)))

if (length(top_genes) < 2) stop("Too few genes. Please revisit the thresholds.")

# -----------------------------------------------------------------
# 3. Get GAM predictions & Z-score normalization
# -----------------------------------------------------------------
cat(sprintf("[3/6] Getting GAM predictions (%d points) & Z-scoring...\n",
            N_PSEUDOTIME_POINTS))

pred_mat  <- predictSmooth(sce_trade, gene = top_genes,
                           nPoints = N_PSEUDOTIME_POINTS, tidy = FALSE)
lin_cols  <- grep(sprintf("^lineage%d_", TARGET_LINEAGE),
                  colnames(pred_mat), value = TRUE)
expr_mat  <- pred_mat[, lin_cols, drop = FALSE]

# Z-score normalization (per gene / row)
expr_scaled <- t(scale(t(expr_mat)))
valid       <- apply(expr_scaled, 1,
                     function(x) !any(is.na(x) | is.infinite(x)))
expr_scaled <- expr_scaled[valid, ]
top_genes   <- top_genes[valid]
cat(sprintf("  Valid genes: %d\n", nrow(expr_scaled)))

# -----------------------------------------------------------------
# 4. Choose number of clusters & k-means
# -----------------------------------------------------------------
cat("[4/6] Clustering...\n")

if (is.null(N_CLUSTERS_OVERRIDE)) {
  cat("  Auto-selecting k by silhouette score...\n")
  sil_scores <- sapply(K_RANGE, function(k) {
    km  <- kmeans(expr_scaled, centers = k, nstart = 25, iter.max = 100)
    sil <- silhouette(km$cluster, dist(expr_scaled))
    mean(sil[, 3])
  })
  best_k <- K_RANGE[which.max(sil_scores)]

  p_sil <- ggplot(data.frame(k = K_RANGE, sil = sil_scores),
                  aes(x = k, y = sil)) +
    geom_line() + geom_point(size = 2) +
    geom_vline(xintercept = best_k, color = "red", linetype = "dashed") +
    annotate("text", x = best_k, y = max(sil_scores),
             label = sprintf("k=%d", best_k),
             color = "red", vjust = -0.5, hjust = -0.2) +
    labs(title = "Silhouette Score", x = "k", y = "Mean silhouette") +
    theme_classic()

  pdf(file.path(out_dir, "silhouette_score.pdf"), width = 6, height = 4)
  print(p_sil)
  dev.off()
  cat(sprintf("  Auto-selected: k=%d (score=%.3f)\n", best_k, max(sil_scores)))
  cat(sprintf("  You can fix it later with N_CLUSTERS_OVERRIDE=%d\n", best_k))

} else {
  best_k <- N_CLUSTERS_OVERRIDE
  cat(sprintf("  Manual: k=%d\n", best_k))
}

set.seed(1234)
km_result  <- kmeans(expr_scaled, centers = best_k, nstart = 25, iter.max = 100)
cluster_id <- km_result$cluster
cat(sprintf("  Number of clusters: %d\n", best_k))
print(table(cluster_id))

# Cluster colors
cluster_colors <- setNames(hcl.colors(best_k, palette = "Set 2"),
                           as.character(1:best_k))

# -----------------------------------------------------------------
# 5. Save result table
# -----------------------------------------------------------------
result_table <- asso_res %>%
  filter(gene %in% top_genes) %>%
  mutate(cluster = cluster_id[match(gene, top_genes)]) %>%
  arrange(cluster, padj, desc(waldStat)) %>%
  select(cluster, gene, waldStat, pvalue, padj)

write.csv(result_table,
          file.path(out_dir, "gene_cluster_table.csv"),
          row.names = FALSE)
cat("  Saved gene_cluster_table.csv\n")

# -----------------------------------------------------------------
# 6. Plots
# -----------------------------------------------------------------
cat("[5/6] Plotting...\n")

pt_axis <- seq(0, 1, length.out = N_PSEUDOTIME_POINTS)

# --- (1) PCA plot (2D view of gene patterns) ---
pca_res <- prcomp(expr_scaled, scale. = FALSE)
pca_df  <- data.frame(
  PC1     = pca_res$x[, 1],
  PC2     = pca_res$x[, 2],
  cluster = factor(cluster_id)
)
var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  labs(
    title = sprintf("PCA of gene expression patterns (k=%d)", best_k),
    x     = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y     = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_classic()

pdf(file.path(out_dir, "pca_gene_clusters.pdf"), width = 7, height = 6)
print(p_pca)
dev.off()
cat("  Saved pca_gene_clusters.pdf\n")

# --- (2) Heatmap ---
gene_order    <- order(cluster_id)
cluster_annot <- data.frame(
  Cluster   = factor(paste0("C", cluster_id[gene_order])),
  row.names = top_genes[gene_order]
)
annot_colors <- list(
  Cluster = setNames(cluster_colors, paste0("C", 1:best_k))
)

pdf(file.path(out_dir, "heatmap_all_genes.pdf"),
    width = 8, height = max(6, nrow(expr_scaled) * 0.015 + 2))
pheatmap(expr_scaled[gene_order, ],
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         annotation_row    = cluster_annot,
         annotation_colors = annot_colors,
         color             = colorRampPalette(
                               c("#053061", "white", "#67001f"))(100),
         labels_col        = c("Early",
                               rep("", N_PSEUDOTIME_POINTS - 2), "Late"),
         main              = sprintf("Gene modules along pseudotime (k=%d)",
                                     best_k))
dev.off()
cat("  Saved heatmap_all_genes.pdf\n")

# --- (3) Representative-gene summary plot (cluster mean) ---
pdf(file.path(out_dir, "representative_genes.pdf"), width = 6, height = 3.5)
for (cl in sort(unique(cluster_id))) {
  cl_genes   <- result_table %>% filter(cluster == cl) %>% pull(gene)
  repr_genes <- result_table %>%
    filter(cluster == cl) %>% arrange(desc(waldStat)) %>%
    head(N_REPR_GENES) %>% pull(gene)

  repr_df <- as.data.frame(expr_scaled[repr_genes, , drop = FALSE]) %>%
    mutate(gene = repr_genes) %>%
    pivot_longer(-gene, names_to = "tp", values_to = "zscore") %>%
    #mutate(pt = rep(pt_axis, each = length(repr_genes)))
    mutate(pt = rep(pt_axis, times = length(repr_genes)))

  mean_df <- data.frame(pt     = pt_axis,
                        zscore = colMeans(expr_scaled[cl_genes, , drop = FALSE]))

  p <- ggplot() +
    geom_line(data = repr_df,
              aes(x = pt, y = zscore, group = gene),
              color = "steelblue", alpha = 0.4, linewidth = 0.5) +
    geom_line(data = mean_df, aes(x = pt, y = zscore),
              color = "black", linewidth = 1.2) +
    labs(title = sprintf("Cluster %d  (n=%d genes, rep=%d shown)",
                         cl, length(cl_genes), length(repr_genes)),
         x = "Pseudotime", y = "Z-score") +
    scale_x_continuous(breaks = c(0, 0.5, 1),
                       labels = c("Early", "Mid", "Late")) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  print(p)
}
dev.off()
cat("  Saved representative_genes.pdf\n")

# --- Cluster mean + SD ribbon plot ---
pdf(file.path(out_dir, "cluster_mean_sd.pdf"), width = 6, height = 3.5)
for (cl in sort(unique(cluster_id))) {
  cl_genes <- result_table %>% filter(cluster == cl) %>% pull(gene)
  cl_mat   <- expr_scaled[cl_genes, , drop = FALSE]

  ribbon_df <- data.frame(
    pt   = pt_axis,
    mean = colMeans(cl_mat),
    sd   = apply(cl_mat, 2, sd)
  ) %>% mutate(ymin = mean - sd, ymax = mean + sd)

  p <- ggplot(ribbon_df, aes(x = pt)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                fill = "grey60", alpha = 0.3) +
    geom_line(aes(y = mean), color = "black", linewidth = 1.2) +
    labs(title = sprintf("Cluster %d  (n=%d genes)", cl, length(cl_genes)),
         x = "Pseudotime", y = "Z-score (mean +/- SD)") +
    scale_x_continuous(breaks = c(0, 0.5, 1),
                       labels = c("Early", "Mid", "Late")) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  print(p)
}
dev.off()
cat("  Saved cluster_mean_sd.pdf\n")

# --- (4) Per-gene dot plots (cells x pseudotime + GAM curve) ---
cat("[6/6] Per-gene dot plots for representative genes...\n")

# Cell data frame on the lineage
cells_on_lin <- which(on_lineage)
pt_cells     <- cell_pt[cells_on_lin]

pdf(file.path(out_dir, "representative_genes_individual.pdf"),
    width = 14, height = 10)

for (cl in sort(unique(cluster_id))) {

  repr_genes <- result_table %>%
    filter(cluster == cl) %>%
    arrange(desc(waldStat)) %>%
    head(N_REPR_GENES) %>%
    pull(gene)

  plot_list <- lapply(repr_genes, function(gene) {

    expr_cells <- as.numeric(logcounts_mat[gene, cells_on_lin])

    df <- data.frame(
      pseudotime = pt_cells,
      expression = expr_cells
    )

    ggplot(df, aes(x = pseudotime, y = expression)) +
      # individual cells (dots)
      geom_point(color = "grey60", alpha = 0.5, size = 0.4) +
      # GAM smoothed curve
      geom_smooth(method    = "gam",
                  formula   = y ~ s(x, bs = "cs"),
		  se        = FALSE,
                  color     = "red",
                  linewidth = 1) +
      labs(title = gene,
           x     = "Pseudotime",
           y     = "Expression (logcounts)") +
      theme_classic() +
      theme(plot.title = element_text(size = 9, face = "bold.italic"))
  })

  n_cols <- min(4, length(repr_genes))
  print(
    wrap_plots(plot_list, ncol = n_cols) +
      plot_annotation(
        title = sprintf("Cluster %d - Representative genes", cl),
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )
}
dev.off()
cat("  Saved representative_genes_individual.pdf\n")

# --- Per-gene dot plots for all genes (one directory per cluster) ---

for (cl in sort(unique(cluster_id))) {

  cl_genes <- result_table %>%
    filter(cluster == cl) %>%
    arrange(desc(waldStat)) %>%
    pull(gene)

  cl_dir <- file.path(out_dir, sprintf("cluster%d_genes", cl))
  dir.create(cl_dir, showWarnings = FALSE)

  cat(sprintf("  Cluster %d: plotting %d genes...\n", cl, length(cl_genes)))

  for (gene in cl_genes) {
    expr_cells <- as.numeric(logcounts_mat[gene, cells_on_lin])

    df <- data.frame(
      pseudotime = pt_cells,
      expression = expr_cells
    )

    p <- ggplot(df, aes(x = pseudotime, y = expression)) +
      geom_point(color = "grey60", alpha = 0.5, size = 0.3) +
      geom_smooth(method    = "gam",
                  formula   = y ~ s(x, bs = "cs"),
                  se        = FALSE,
                  color     = "red",
                  linewidth = 1) +
      labs(title = gene,
           x     = "Pseudotime",
           y     = "Expression (logcounts)") +
      theme_classic() +
      theme(plot.title = element_text(size = 9, face = "bold.italic"))

    pdf(file.path(cl_dir, sprintf("%s.pdf", gene)), width = 4, height = 3)
    print(p)
    dev.off()
  }

  cat(sprintf("  Cluster %d done: %s/\n", cl, cl_dir))
}

cat(sprintf("\n=== Done ===\nOutput: %s/\n", out_dir))
