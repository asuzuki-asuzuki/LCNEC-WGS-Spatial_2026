#!/usr/bin/env Rscript
# =============================================================================
# 06_tradeseq_pseudotime.R
# Identify pseudotime-associated genes with tradeSeq (analyzed per lineage)
# Output: 06_tradeseq_result/siRNA_1/group_1/
# =============================================================================

library(tradeSeq)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
set.seed(1234)

# =====================================================================
# User settings
# =====================================================================

SAMPLE_CONFIG <- list(
  siRNA_1 = list(subset_pattern = "group_1")
)

RUN_SAMPLES <- c("siRNA_1")

# --- Number of knots ---
USE_EVALUATE_K <- FALSE
K_FIXED        <- 6
K_RANGE        <- 3:10

# --- Filtering ---
MIN_EXPR_CELLS <- 20     # keep genes expressed in at least this many cells

# --- Thresholds ---
PADJ_THRESHOLD <- 0.05
EXPR_RANGE_MIN <- 0.5    # minimum expr_range for UMAP visualization (log1p scale)

# --- Number of genes for UMAP visualization (expr_range >= EXPR_RANGE_MIN & top waldStat) ---
TOP_N_GENES <- 12

# --- Number of pseudotime points for expr_range ---
N_POINTS_EXPR_RANGE <- 10

# --- Parallelization ---
USE_PARALLEL <- TRUE
N_CORES      <- 4

INPUT_BASE  <- "05_trajectory_result"
OUTPUT_BASE <- "06_tradeseq_result"

# =====================================================================
# Parallelization setup
# =====================================================================

library(BiocParallel)
BPPARAM <- if (USE_PARALLEL) MulticoreParam(N_CORES) else SerialParam()

# =====================================================================
# Main
# =====================================================================

for (sample_name in RUN_SAMPLES) {

  cat(sprintf("\n=== %s ===\n", sample_name))

  config         <- SAMPLE_CONFIG[[sample_name]]
  subset_pattern <- config$subset_pattern

  input_rds <- file.path(INPUT_BASE, sample_name,
                          subset_pattern, "sce_slingshot.rds")

  if (!file.exists(input_rds)) {
    cat(sprintf("  Skip: %s not found\n", input_rds))
    next
  }

  out_dir <- file.path(OUTPUT_BASE, sample_name, subset_pattern)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # -----------------------------------------------------------------
  # 1. Load SCE & inspect
  # -----------------------------------------------------------------
  cat("  [1/6] Loading SCE...\n")
  sce <- readRDS(input_rds)
  cat(sprintf("    %d cells, %d genes\n", ncol(sce), nrow(sce)))

  n_lineages  <- length(slingLineages(sce))
  pseudo_mat  <- as.matrix(slingPseudotime(sce))
  umap_coords <- reducedDims(sce)$UMAP.PCA

  cat(sprintf("    Number of lineages: %d\n", n_lineages))
  for (i in seq_len(n_lineages)) {
    cat(sprintf("    Lineage %d: %d cells\n", i,
                sum(!is.na(pseudo_mat[, i]))))
  }

  # -----------------------------------------------------------------
  # 2. Filter low-expression genes & remove mito/ribosomal genes
  # -----------------------------------------------------------------
  cat(sprintf("  [2/6] Filtering low-expression genes (>= %d cells)...\n",
              MIN_EXPR_CELLS))
  counts_mat <- assay(sce, "counts")
  expressed  <- rowSums(counts_mat > 0) >= MIN_EXPR_CELLS
  sce_filt   <- sce[expressed, ]
  cat(sprintf("    %d -> %d genes\n", nrow(sce), nrow(sce_filt)))

  # Remove mitochondrial / ribosomal genes
  exclude  <- grepl("^MT-|^RPL|^RPS", rownames(sce_filt))
  sce_filt <- sce_filt[!exclude, ]
  cat(sprintf("    After mito/ribosomal removal: %d genes\n", nrow(sce_filt)))

  # -----------------------------------------------------------------
  # 3. Choose number of knots & fitGAM
  # -----------------------------------------------------------------
  cat("  [3/6] Running fitGAM...\n")

  if (USE_EVALUATE_K) {
    cat(sprintf("    Running evaluateK (k=%d-%d)...\n",
                min(K_RANGE), max(K_RANGE)))
    set.seed(1234)
    eval_k <- evaluateK(
      counts  = assay(sce_filt, "counts"),
      sds     = SlingshotDataSet(sce_filt),
      k       = K_RANGE,
      nGenes  = 500,
      verbose = FALSE,
      BPPARAM = BPPARAM
    )
    best_k <- K_RANGE[which.min(colMeans(eval_k, na.rm = TRUE))]
    cat(sprintf("    Best number of knots: %d\n", best_k))

    pdf(file.path(out_dir, "evaluateK.pdf"), width = 7, height = 5)
    plot(K_RANGE, colMeans(eval_k, na.rm = TRUE), type = "b",
         xlab = "Number of knots", ylab = "Mean AIC",
         main = "evaluateK result")
    abline(v = best_k, col = "red", lty = 2)
    dev.off()

  } else {
    best_k <- K_FIXED
    cat(sprintf("    Number of knots: %d (fixed)\n", best_k))
  }

  cat(sprintf("    Running fitGAM (k=%d)...\n", best_k))
  sce_trade <- fitGAM(
    counts  = assay(sce_filt, "counts"),
    sds     = SlingshotDataSet(sce_filt),
    nknots  = best_k,
    verbose = FALSE,
    BPPARAM = BPPARAM
  )
  cat("    fitGAM done\n")
  saveRDS(sce_trade, file.path(out_dir, "sce_tradeseq.rds"))

  # -----------------------------------------------------------------
  # 4. Compute expr_range (all genes, log1p scale)
  # -----------------------------------------------------------------
  cat(sprintf("  [4/6] Computing expr_range (all genes x %d points, log1p)...\n",
              N_POINTS_EXPR_RANGE))

  pred_mat_all <- predictSmooth(sce_trade,
                                gene    = rownames(sce_trade),
                                nPoints = N_POINTS_EXPR_RANGE,
                                tidy    = FALSE)

  # Per lineage: max - min after log1p transform
  expr_range_list <- lapply(seq_len(n_lineages), function(i) {
    cols <- grep(sprintf("^lineage%d_", i), colnames(pred_mat_all),
                 value = TRUE)
    if (length(cols) == 0) return(NULL)
    rng <- apply(log1p(pred_mat_all[, cols, drop = FALSE]), 1,
                 function(x) max(x) - min(x))
    data.frame(gene = names(rng), expr_range = rng, lineage = i)
  })
  rm(pred_mat_all); gc()

  # -----------------------------------------------------------------
  # 5. associationTest & CSV (all genes, per lineage)
  # -----------------------------------------------------------------
  cat("  [5/6] Running associationTest...\n")
  asso_all <- associationTest(sce_trade, lineages = TRUE)

  for (i in seq_len(n_lineages)) {

    cat(sprintf("    Processing lineage %d...\n", i))

    waldStat_col <- sprintf("waldStat_%d", i)
    df_col       <- sprintf("df_%d",       i)
    pval_col     <- sprintf("pvalue_%d",   i)

    if (!waldStat_col %in% colnames(asso_all)) {
      cat(sprintf("    -> No result for lineage %d (skip)\n", i))
      next
    }

    range_df_i <- expr_range_list[[i]]

    res_i <- asso_all %>%
      select(all_of(c(waldStat_col, df_col, pval_col))) %>%
      rename(waldStat = !!waldStat_col,
             df       = !!df_col,
             pvalue   = !!pval_col) %>%
      mutate(
        gene    = rownames(asso_all),
        padj    = p.adjust(pvalue, method = "BH"),
        lineage = i
      ) %>%
      filter(!is.na(pvalue)) %>%
      left_join(range_df_i %>% select(gene, expr_range), by = "gene") %>%
      arrange(desc(waldStat)) %>%
      select(gene, expr_range, waldStat, df, pvalue, padj, lineage)

    n_sig <- sum(res_i$padj < PADJ_THRESHOLD, na.rm = TRUE)
    cat(sprintf("    -> Significant genes (padj<%.2f): %d\n",
                PADJ_THRESHOLD, n_sig))

    # Save all significant genes to CSV (sorted by expr_range)
    res_sig <- res_i %>% filter(padj < PADJ_THRESHOLD)

    write.csv(res_sig,
              file      = file.path(out_dir,
                                    sprintf("lineage%d_association.csv", i)),
              row.names = FALSE)
    cat(sprintf("    -> CSV saved: %d genes\n", nrow(res_sig)))

    # -----------------------------------------------------------------
    # 6. UMAP visualization (expr_range >= threshold & top waldStat)
    # -----------------------------------------------------------------
    cat(sprintf("    Plotting top %d genes for lineage %d...\n", TOP_N_GENES, i))

    top_genes <- res_sig %>%
      filter(expr_range >= EXPR_RANGE_MIN) %>%
      arrange(desc(waldStat)) %>%
      head(TOP_N_GENES) %>%
      pull(gene)

    if (length(top_genes) == 0) {
      cat(sprintf("    -> No significant gene with expr_range >= %.1f. Skip\n",
                  EXPR_RANGE_MIN))
      next
    }
    cat(sprintf("    -> Candidates: expr_range>=%.1f & top %d by waldStat\n",
                EXPR_RANGE_MIN, length(top_genes)))

    logcounts_mat <- assay(sce_filt, "logcounts")
    on_lineage    <- !is.na(pseudo_mat[, i])
    n_cols        <- min(4, TOP_N_GENES)
    n_rows        <- ceiling(length(top_genes) / n_cols)

    # Combined UMAP
    pdf(file.path(out_dir, sprintf("lineage%d_top%d_umap.pdf", i, TOP_N_GENES)),
        width  = n_cols * 4,
        height = n_rows * 4)

    plot_list <- lapply(top_genes, function(gene) {
      expr_vec       <- as.numeric(logcounts_mat[gene, ])
      df <- data.frame(
        UMAP1  = umap_coords[, 1],
        UMAP2  = umap_coords[, 2],
        expr   = expr_vec,
        on_lin = on_lineage
      ) %>% arrange(expr)

      ggplot() +
        geom_point(data = df[!df$on_lin, ],
                   aes(x = UMAP1, y = UMAP2),
                   color = "grey85", size = 0.4) +
        geom_point(data = df[df$on_lin, ],
                   aes(x = UMAP1, y = UMAP2, color = expr),
                   size = 0.6) +
        scale_color_viridis(option = "A", name = "logcounts") +
        ggtitle(gene) +
        theme_classic() +
        theme(
          plot.title   = element_text(size = 9, face = "bold.italic"),
          axis.title   = element_text(size = 8),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8)
        )
    })

    print(wrap_plots(plot_list, ncol = n_cols))
    dev.off()

    # Per-gene PDF
    gene_dir <- file.path(out_dir, sprintf("lineage%d_individual_genes", i))
    dir.create(gene_dir, showWarnings = FALSE)

    for (gene in top_genes) {
      expr_vec       <- as.numeric(logcounts_mat[gene, ])
      expr_range_val <- res_sig$expr_range[res_sig$gene == gene]
      df <- data.frame(
        UMAP1  = umap_coords[, 1],
        UMAP2  = umap_coords[, 2],
        expr   = expr_vec,
        on_lin = on_lineage
      ) %>% arrange(expr)

      p <- ggplot() +
        geom_point(data = df[!df$on_lin, ],
                   aes(x = UMAP1, y = UMAP2),
                   color = "grey85", size = 0.5) +
        geom_point(data = df[df$on_lin, ],
                   aes(x = UMAP1, y = UMAP2, color = expr),
                   size = 0.8) +
        scale_color_viridis(option = "A", name = "logcounts") +
        ggtitle(sprintf("%s (Lineage %d, expr_range=%.2f)",
                        gene, i, expr_range_val)) +
        theme_classic()

      pdf(file.path(gene_dir, sprintf("%s.pdf", gene)), width = 6, height = 5)
      print(p)
      dev.off()
    }

    cat(sprintf("    Lineage %d done\n", i))
  }

  cat(sprintf("  Done: %s/\n", out_dir))
}

cat("\n=== Done ===\n")
cat(sprintf("Output: %s/\n", OUTPUT_BASE))
