#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 1: MAIN GENE-LEVEL ANALYSIS (CELLS)
## =========================================================
## This script performs:
## - Quality control (density, PCA, MDS)
## - Differential expression (DESeq2)
## - Volcano plots
## - Heatmaps (top DEGs)
## - ComplexHeatmap
##
## NO PATHWAY ANALYSIS (see part2_cells_ora.R and part3_cells_fgsea.R)
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "viridis",
  "grid",
  "scales"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "DESeq2",
  "ComplexHeatmap"
)

## circlize is on CRAN (required for ComplexHeatmap)
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(viridis)
  library(grid)   ## for unit()
  library(scales)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
})

## =======================================================
## UNIFIED PLOT THEME & COLORS
## =======================================================
## Consistent styling across all plots (matches tissue analysis)

## Color palette
plot_colors <- list(
  genotype = c(CTL = "#0072B2", KAT8KD = "#D55E00"),  # colorblind-friendly blue/orange
  viridis_option = "viridis",
  ns_grey = "grey75",
  text_black = "black",
  line_grey = "grey40"
)

## Unified ggplot2 theme
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      ## Panel
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background = element_blank(),

      ## Axes
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),

      ## Title
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "grey30"),

      ## Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.position = "right",

      ## Strip (for facets)
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = base_size - 1, face = "bold")
    )
}

## Volcano label settings
volcano_label_n_fdr <- 10  # Top N genes by FDR per direction
volcano_label_n_fc <- 10   # Top N genes by |logFC| per direction

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
source(file.path(utils_dir, "save.R"))
source(file.path(utils_dir, "findsave.R"))
source(file.path(utils_dir, "purge.R"))
source(file.path(utils_dir, "dedupe.R"))

## 2.5) Run metadata ----------------------------------------
run_info <- list(
  script_name = "part1_cells_main_analysis.R",
  species     = "mouse",
  data_type   = "bulkRNAseq_3T3L1_cells",
  keywords    = c("KAT8", "3T3-L1", "adipocytes", "bulk RNA-seq", "cells", "DESeq2"),
  notes       = "Part 1: Gene-level analysis only (no pathway analysis)",
  message     = "Cells: DESeq2 + VST QC + DEGs + ComplexHeatmap"
)

## 3) init_run ----------------------------------------------
run_ctx <- init_run(
  script_name = run_info$script_name,
  species     = run_info$species,
  data_type   = run_info$data_type,
  keywords    = run_info$keywords,
  notes       = run_info$notes
)

outdir  <- run_ctx$outdir
run_tag <- run_ctx$run_tag
cat("Run directory:", normalizePath(outdir, mustWork = FALSE), "\n")

## -----------------------------
## Prefilter parameters
## -----------------------------
prefilter_min_count <- 10
prefilter_min_samples_per_group <- 3

## 4) Main logic --------------------------------------------
hero_volcano_file <- NULL

tryCatch({

  ## 4.1) Load raw counts -----------------------------------
  cat("\n=== RAW COUNTS MATRIX (FULL) ===\n")

  counts_raw <- read.delim(
    "counts.txt",
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  cat("Dimensions (rows x cols):", paste(dim(counts_raw), collapse = " x "), "\n")
  cat("First 10 column names:\n")
  print(head(colnames(counts_raw), 10))
  cat("\nHead of counts (first 5 rows, first 6 cols):\n")
  print(head(counts_raw[, 1:min(6, ncol(counts_raw))]))

  ## 4.2) Define cell samples + annotation ------------------
  ## JS_41-JS_44 = CTL, JS_45-JS_48 = KAT8KD
  cell_samples <- paste0("JS_", sprintf("%02d", 41:48))

  sample_annot <- data.frame(
    Sample   = cell_samples,
    Genotype = c(rep("CTL", 4), rep("KAT8KD", 4)),
    stringsAsFactors = FALSE
  )

  ## Check samples exist
  missing_cells <- setdiff(cell_samples, colnames(counts_raw))
  if (length(missing_cells) > 0) {
    stop("These cell sample columns are missing in counts.txt: ", paste(missing_cells, collapse = ", "))
  }

  counts_cells_raw <- counts_raw[, c("gene_id", "gene", cell_samples)]

  cat("\n=== CELL SUBSET ===\n")
  cat("Dimensions (rows x cols):", paste(dim(counts_cells_raw), collapse = " x "), "\n")
  cat("Sample annotation:\n")
  print(sample_annot)

  ## 4.3) Gene cleaning -------------------------------------
  counts_cells <- counts_cells_raw %>%
    mutate(
      ensembl_gene_id = gene_id,
      gene_name       = gene
    )

  na_or_blank <- is.na(counts_cells$gene_name) | counts_cells$gene_name == ""
  counts_cells$gene_name[na_or_blank] <- counts_cells$ensembl_gene_id[na_or_blank]

  ## Handle duplicates robustly
  if (any(duplicated(counts_cells$gene_name))) {
    dup_flag <- duplicated(counts_cells$gene_name) | duplicated(counts_cells$gene_name, fromLast = TRUE)
    counts_cells$gene_name[dup_flag] <- paste0(counts_cells$ensembl_gene_id[dup_flag], "_", counts_cells$gene_name[dup_flag])
  }
  ## As a final safety net:
  counts_cells$gene_name <- make.unique(counts_cells$gene_name)

  rownames(counts_cells) <- counts_cells$gene_name

  cat("\nAny NA in gene_name? ", any(is.na(counts_cells$gene_name)), "\n")
  cat("Any blank gene_name? ", any(counts_cells$gene_name == ""), "\n")
  cat("Any duplicated gene_name? ", any(duplicated(counts_cells$gene_name)), "\n")

  ## 4.4) Build count matrix --------------------------------
  count_matrix_cells <- as.matrix(counts_cells[, cell_samples])
  rownames(count_matrix_cells) <- counts_cells$gene_name

  ## 4.5) Attach sample annotation --------------------------
  sample_annot <- sample_annot[match(colnames(count_matrix_cells), sample_annot$Sample), ]
  if (any(is.na(sample_annot$Sample))) stop("Sample annotation failed to match count matrix columns.")

  sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
  rownames(sample_annot) <- sample_annot$Sample

  cat("\n=== BASIC QC (CELLS) ===\n")
  cat("Number of genes:   ", nrow(count_matrix_cells), "\n")
  cat("Number of samples: ", ncol(count_matrix_cells), "\n")
  cat("Library sizes (colSums raw counts):\n")
  print(colSums(count_matrix_cells))
  cat("\nSummary of library sizes:\n")
  print(summary(colSums(count_matrix_cells)))

  ## 4.6) DESeq2 object + typical prefilter ------------------
  min_count <- prefilter_min_count
  min_samples_per_group <- prefilter_min_samples_per_group

  keep <- rowSums(count_matrix_cells >= min_count) >= min_samples_per_group

  ## Comprehensive filtering sanity check
  n_before <- nrow(count_matrix_cells)
  n_after <- sum(keep)
  n_removed <- n_before - n_after
  pct_retained <- round(100 * n_after / n_before, 1)

  cat("\n==========================================================\n")
  cat("=== SANITY CHECK: GENE PREFILTERING ===\n")
  cat("==========================================================\n")
  cat("Prefilter rule: counts >= ", min_count, " in >= ", min_samples_per_group, " samples\n", sep = "")
  cat("Genes BEFORE prefiltering:  ", n_before, "\n", sep = "")
  cat("Genes AFTER prefiltering:   ", n_after, "\n", sep = "")
  cat("Genes REMOVED:              ", n_removed, " (", 100 - pct_retained, "%)\n", sep = "")
  cat("Genes RETAINED:             ", pct_retained, "%\n", sep = "")
  cat("==========================================================\n\n")

  count_matrix_filt <- count_matrix_cells[keep, , drop = FALSE]

  dds_cells <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filt),
    colData   = sample_annot,
    design    = ~ Genotype
  )

  ## Save centralized parameters (matching tissue analysis thresholds)
  logFC_cut_cells <- 1      # Same as tissue: |log2FC| > 1
  fdr_cut_cells   <- 0.05   # Same as tissue: FDR < 0.05

  params <- list(
    prefilter_rule    = list(
      min_count       = min_count,
      min_samples     = min_samples_per_group,
      description     = "Keep genes with counts >= min_count in at least min_samples samples"
    ),
    de_thresholds     = list(
      logFC_cutoff    = logFC_cut_cells,
      fdr_cutoff      = fdr_cut_cells,
      description     = "Same thresholds used for volcano, heatmap, and pathway analysis"
    ),
    vst              = list(
      blind           = TRUE,
      description     = "VST (variance stabilizing transformation) for QC visualization"
    ),
    volcano_labels   = list(
      top_n_by_fdr    = volcano_label_n_fdr,
      top_n_by_fc     = volcano_label_n_fc
    ),
    heatmap_max_genes = 300
  )

  params_file <- file.path(outdir, "logs", paste0("params_", run_tag, ".txt"))
  params_lines <- c(
    "=== RNA-seq Pipeline Parameters (Part 1: Gene-level Analysis, Cells) ===",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Prefilter Rule (count-based) ---",
    paste0("  min_count: ", params$prefilter_rule$min_count),
    paste0("  min_samples: ", params$prefilter_rule$min_samples),
    paste0("  Description: ", params$prefilter_rule$description),
    "",
    "--- DE Thresholds (unified across all analyses) ---",
    paste0("  |log2FC| cutoff: ", params$de_thresholds$logFC_cutoff),
    paste0("  FDR cutoff: ", params$de_thresholds$fdr_cutoff),
    paste0("  Note: ", params$de_thresholds$description),
    "",
    "--- VST (DESeq2 normalization for visualization) ---",
    paste0("  blind: ", params$vst$blind),
    paste0("  Description: ", params$vst$description),
    "",
    "--- Volcano Plot Labels ---",
    paste0("  Top N by FDR (per direction): ", params$volcano_labels$top_n_by_fdr),
    paste0("  Top N by |logFC| (per direction): ", params$volcano_labels$top_n_by_fc),
    "",
    "--- Heatmap ---",
    paste0("  Max genes displayed: ", params$heatmap_max_genes)
  )
  writeLines(params_lines, con = params_file)
  cat("Saved params to ", params_file, "\n")

  ## 4.7) Run DESeq2 ----------------------------------------
  cat("\n=== RUNNING DESeq2 ===\n")
  dds_cells <- DESeq(dds_cells)
  cat("[OK] DESeq2 complete.\n")

  cat("\nDESeq2 design matrix:\n")
  print(model.matrix(design(dds_cells), colData(dds_cells)))

  ## =======================================================
  ## 4.7.1) DESeq2 SANITY CHECKS
  ## =======================================================
  cat("\n--- DESeq2 SANITY CHECKS ---\n")

  ## Size factors (should be near 1, large deviations indicate library size issues)
  sf <- sizeFactors(dds_cells)
  cat("\nSize factors (should be ~1.0, range indicates library size variation):\n")
  print(round(sf, 3))
  cat("  Range: ", round(min(sf), 3), " - ", round(max(sf), 3), "\n", sep = "")
  cat("  Fold difference: ", round(max(sf)/min(sf), 2), "x\n", sep = "")

  if (max(sf)/min(sf) > 3) {
    cat("[WARN] Large variation in size factors (>3x). Check library sizes.\n")
  } else {
    cat("[OK] Size factors within normal range.\n")
  }

  ## Dispersion estimates
  cat("\nDispersion estimate summary:\n")
  disp <- dispersions(dds_cells)
  disp <- disp[!is.na(disp)]
  cat("  Median dispersion: ", round(median(disp), 4), "\n", sep = "")
  cat("  Range: ", round(min(disp), 4), " - ", round(max(disp), 4), "\n", sep = "")
  cat("  Genes with dispersion > 1: ", sum(disp > 1), " (", round(100*sum(disp > 1)/length(disp), 1), "%)\n", sep = "")

  ## Save dispersion plot
  disp_file <- paste0("Dispersion_plot_cells_", run_tag, ".png")
  png(file.path(outdir, "plots", disp_file), width = 1000, height = 800, res = 150)
  plotDispEsts(dds_cells, main = "DESeq2 Dispersion Estimates (Cells)")
  dev.off()
  ## Save PDF version
  disp_file_pdf <- paste0("Dispersion_plot_cells_", run_tag, ".pdf")
  pdf(file.path(outdir, "plots", disp_file_pdf), width = 6.67, height = 5.33)
  plotDispEsts(dds_cells, main = "DESeq2 Dispersion Estimates (Cells)")
  dev.off()
  cat("[OK] Saved dispersion plot: ", disp_file, " (PNG and PDF)\n", sep = "")

  ## Cook's distance - identify potential outlier samples
  cat("\nCook's distance (sample-level outlier detection):\n")
  cooks <- assays(dds_cells)[["cooks"]]
  if (!is.null(cooks)) {
    max_cooks_per_sample <- apply(cooks, 2, function(x) max(x, na.rm = TRUE))
    cat("  Max Cook's distance per sample:\n")
    print(round(max_cooks_per_sample, 3))
    high_cooks <- names(max_cooks_per_sample)[max_cooks_per_sample > 10]
    if (length(high_cooks) > 0) {
      cat("[WARN] Samples with high Cook's distance (>10): ", paste(high_cooks, collapse = ", "), "\n", sep = "")
    } else {
      cat("[OK] No samples with extreme Cook's distance.\n")
    }
  }

  ## 4.8) VST transform for QC ------------------------------
  cat("\n=== VST FOR QC ===\n")
  vst_cells <- vst(dds_cells, blind = TRUE)
  vst_mat <- assay(vst_cells)  ## genes x samples (log2-ish stabilized)

  cat("\n=== VST FOR HEATMAPS (blind=FALSE) ===\n")
  vst_cells_hm <- vst(dds_cells, blind = FALSE)
  vst_mat_hm <- assay(vst_cells_hm)

  ## 4.9) Density plots: BEFORE vs AFTER DESeq2 normalization ---
  cat("\n--- DENSITY PLOTS: BEFORE vs AFTER DESeq2 NORMALIZATION ---\n")

  ## Get VST for all genes (before filtering) for comparison
  dds_all_genes <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_cells),
    colData   = sample_annot,
    design    = ~ Genotype
  )
  dds_all_genes <- estimateSizeFactors(dds_all_genes)
  dds_all_genes <- estimateDispersions(dds_all_genes, fitType = "local")
  vst_all_genes <- vst(dds_all_genes, blind = TRUE)
  vst_mat_before <- assay(vst_all_genes)

  plot_density_vst_comparison <- function(vst_before, vst_after) {
    dens_before <- apply(vst_before, 2, density)
    dens_after  <- apply(vst_after, 2, density)

    xlim_before <- range(sapply(dens_before, function(d) range(d$x)))
    xlim_after  <- range(sapply(dens_after, function(d) range(d$x)))
    xlim_global <- range(c(xlim_before, xlim_after))

    max_y_before <- max(sapply(dens_before, function(d) max(d$y)))
    max_y_after  <- max(sapply(dens_after, function(d) max(d$y)))

    par(mfrow = c(1, 2), mar = c(5, 5, 5, 2), oma = c(0, 0, 2, 0))

    ## Define consistent colors
    sample_colors <- rep(c("#0072B2", "#D55E00"), each = 4)

    ## BEFORE filter
    plot(dens_before[[1]],
         main = "Before Filtering",
         sub  = paste0("n=", nrow(vst_before), " genes, ", ncol(vst_before), " samples"),
         xlab = "VST value",
         lwd  = 2,
         col  = sample_colors[1],
         xlim = xlim_global,
         ylim = c(0, max_y_before * 1.15))
    if (ncol(vst_before) > 1) {
      for (i in 2:ncol(vst_before)) {
        lines(dens_before[[i]], lwd = 2, col = sample_colors[i])
      }
    }

    ## AFTER filter
    plot(dens_after[[1]],
         main = "After Filtering",
         sub  = paste0("n=", nrow(vst_after), " genes, ", ncol(vst_after), " samples"),
         xlab = "VST value",
         lwd  = 2,
         col  = sample_colors[1],
         xlim = xlim_global,
         ylim = c(0, max_y_after * 1.15))
    if (ncol(vst_after) > 1) {
      for (i in 2:ncol(vst_after)) {
        lines(dens_after[[i]], lwd = 2, col = sample_colors[i])
      }
    }
  }

  density_file <- paste0("VST_density_cells_before_after_", run_tag, ".png")
  png(file.path(outdir, "plots", density_file), width = 1600, height = 900, res = 150)
  plot_density_vst_comparison(vst_mat_before, vst_mat)
  dev.off()
  ## Save PDF version
  density_file_pdf <- paste0("VST_density_cells_before_after_", run_tag, ".pdf")
  pdf(file.path(outdir, "plots", density_file_pdf), width = 10.67, height = 6)
  plot_density_vst_comparison(vst_mat_before, vst_mat)
  dev.off()
  cat("Saved VST density plot (before/after filtering) to ", file.path(outdir, "plots", density_file), " (PNG and PDF)\n")

  ## Sample correlation heatmap (VST)
  cat("\n--- SAMPLE CORRELATION MATRIX ---\n")
  sample_cor <- cor(vst_mat, method = "pearson")
  cat("Sample-to-sample Pearson correlations (VST):\n")
  print(round(sample_cor, 3))
  cat("\n  Min correlation: ", round(min(sample_cor[lower.tri(sample_cor)]), 3), "\n", sep = "")
  cat("  Mean correlation: ", round(mean(sample_cor[lower.tri(sample_cor)]), 3), "\n", sep = "")

  if (min(sample_cor[lower.tri(sample_cor)]) < 0.8) {
    cat("[WARN] Some samples have low correlation (<0.8). Check for outliers.\n")
  } else {
    cat("[OK] All sample correlations >= 0.8\n")
  }

  ## 4.10) PCA (VST) ----------------------------------------
  pca_cells <- prcomp(t(vst_mat), scale. = TRUE)
  pca_var <- summary(pca_cells)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    Sample   = colnames(vst_mat),
    Genotype = colData(dds_cells)$Genotype,
    PC1      = pca_cells$x[, 1],
    PC2      = pca_cells$x[, 2],
    stringsAsFactors = FALSE
  )

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Genotype, fill = Genotype, shape = Genotype, label = Sample)) +
    geom_point(size = 4.5, stroke = 1) +
    geom_text_repel(size = 3.5, max.overlaps = 20, segment.color = "grey50", color = "black") +
    scale_color_manual(values = plot_colors$genotype, name = "Genotype") +
    scale_fill_manual(values = plot_colors$genotype, name = "Genotype") +
    scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
    labs(
      title = "PCA: 3T3-L1 Cells (DESeq2 VST)",
      x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_var[2], 1), "%)")
    ) +
    theme_publication(base_size = 12) +
    guides(
      color = "none",
      fill = guide_legend(override.aes = list(shape = 21, color = unname(plot_colors$genotype), fill = unname(plot_colors$genotype))),
      shape = guide_legend(override.aes = list(fill = "white", color = "black"))
    )

  pca_file <- paste0("PCA_cells_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", pca_file), plot = p_pca, width = 7, height = 6, dpi = 300)
  ## Save PDF version
  pca_file_pdf <- paste0("PCA_cells_VST_", run_tag, ".pdf")
  ggsave(file.path(outdir, "plots", pca_file_pdf), plot = p_pca, width = 7, height = 6, device = "pdf")
  cat("Saved PCA plot to ", file.path(outdir, "plots", pca_file), " (PNG and PDF)\n")

  ## 4.11) MDS (VST) ----------------------------------------
  dist_mat <- dist(t(vst_mat))
  mds_coords <- cmdscale(dist_mat, k = 2)

  mds_df <- data.frame(
    Sample   = colnames(vst_mat),
    Genotype = colData(dds_cells)$Genotype,
    MDS1     = mds_coords[, 1],
    MDS2     = mds_coords[, 2],
    stringsAsFactors = FALSE
  )

  p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Genotype, fill = Genotype, shape = Genotype, label = Sample)) +
    geom_point(size = 4.5, stroke = 1) +
    geom_text_repel(size = 3.5, max.overlaps = 20, segment.color = "grey50", color = "black") +
    scale_color_manual(values = plot_colors$genotype, name = "Genotype") +
    scale_fill_manual(values = plot_colors$genotype, name = "Genotype") +
    scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
    labs(title = "MDS: 3T3-L1 Cells (DESeq2 VST)", x = "MDS1", y = "MDS2") +
    theme_publication(base_size = 12) +
    guides(
      color = "none",
      fill = guide_legend(override.aes = list(shape = 21, color = unname(plot_colors$genotype), fill = unname(plot_colors$genotype))),
      shape = guide_legend(override.aes = list(fill = "white", color = "black"))
    )

  mds_file <- paste0("MDS_cells_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", mds_file), plot = p_mds, width = 7, height = 6, dpi = 300)
  ## Save PDF version
  mds_file_pdf <- paste0("MDS_cells_VST_", run_tag, ".pdf")
  ggsave(file.path(outdir, "plots", mds_file_pdf), plot = p_mds, width = 7, height = 6, device = "pdf")
  cat("Saved MDS plot to ", file.path(outdir, "plots", mds_file), " (PNG and PDF)\n")

  ## 4.12) Differential expression: KAT8KD vs CTL -----------
  cat("\n=== DIFFERENTIAL EXPRESSION: KAT8KD vs CTL ===\n")

  ## DESeq2 results (KAT8KD vs CTL, since CTL is reference level)
  res <- results(dds_cells, contrast = c("Genotype", "KAT8KD", "CTL"))
  res <- res[order(res$pvalue), , drop = FALSE]

  tt <- as.data.frame(res)
  tt$logFC     <- tt$log2FoldChange
  tt$P.Value   <- tt$pvalue
  tt$adj.P.Val <- tt$padj
  tt$gene_name <- rownames(tt)

  ## Remove duplicate columns
  tt <- tt[, !duplicated(colnames(tt)), drop = FALSE]

  de_file <- paste0("DE_cells_KAT8KD_vs_CTL_", run_tag, ".csv")
  write.csv(tt, file = file.path(outdir, "tables", de_file), row.names = TRUE)

  ## Calculate comprehensive DEG statistics
  n_total <- nrow(tt)
  n_valid_pval <- sum(!is.na(tt$adj.P.Val))
  n_fdr_only <- sum(tt$adj.P.Val < fdr_cut_cells, na.rm = TRUE)
  n_up_both <- sum(tt$adj.P.Val < fdr_cut_cells & tt$logFC > logFC_cut_cells, na.rm = TRUE)
  n_down_both <- sum(tt$adj.P.Val < fdr_cut_cells & tt$logFC < -logFC_cut_cells, na.rm = TRUE)
  n_both <- n_up_both + n_down_both

  cat("\n==========================================================\n")
  cat("=== SANITY CHECK: KAT8KD_vs_CTL ===\n")
  cat("==========================================================\n")
  cat("Total genes tested by DESeq2:           ", n_total, "\n", sep = "")
  cat("Genes with valid adjusted p-value:      ", n_valid_pval, "\n", sep = "")
  cat("\n--- DEG Counts at Different Thresholds ---\n")
  cat("DEGs (FDR < ", fdr_cut_cells, " only):               ", n_fdr_only, "\n", sep = "")
  cat("DEGs (FDR < ", fdr_cut_cells, " AND |logFC| > ", logFC_cut_cells, "):  ", n_both, "\n", sep = "")
  cat("  - Up-regulated (logFC > ", logFC_cut_cells, "):       ", n_up_both, "\n", sep = "")
  cat("  - Down-regulated (logFC < -", logFC_cut_cells, "):    ", n_down_both, "\n", sep = "")
  cat("\nPercent of tested genes that are DEGs:\n")
  cat("  - FDR-only threshold:                 ", round(100 * n_fdr_only / n_total, 2), "%\n", sep = "")
  cat("  - FDR + logFC threshold:              ", round(100 * n_both / n_total, 2), "%\n", sep = "")
  cat("==========================================================\n")
  cat("DE table written: ", file.path(outdir, "tables", de_file), "\n", sep = "")
  cat("==========================================================\n\n")

  ## DEG summary
  deg_summary <- data.frame(
    Contrast          = "KAT8KD_vs_CTL",
    N_DEG_FDR_lt_0_05 = n_fdr_only,
    N_Up              = n_up_both,
    N_Down            = n_down_both,
    stringsAsFactors  = FALSE
  )

  deg_summary_file <- paste0("DEG_summary_cells_", run_tag, ".csv")
  write.csv(deg_summary, file = file.path(outdir, "tables", deg_summary_file), row.names = FALSE)
  cat("\nDEG summary written: ", file.path(outdir, "tables", deg_summary_file), "\n", sep = "")

  ## 4.13) Volcano plot -------------------------------------
  cat("\n--- VOLCANO PLOT ---\n")

  ## For volcano: only keep genes with valid padj (tested genes)
  tt_plot <- tt %>%
    dplyr::filter(is.finite(logFC), is.finite(adj.P.Val), adj.P.Val > 0) %>%
    mutate(
      negLogFDR = -log10(adj.P.Val),
      sig_cat = dplyr::case_when(
        adj.P.Val < fdr_cut_cells & logFC >=  logFC_cut_cells  ~ "Up",
        adj.P.Val < fdr_cut_cells & logFC <= -logFC_cut_cells  ~ "Down",
        TRUE                                                    ~ "NS"
      ),
      sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
    )

  cat("Genes in volcano plot:       ", nrow(tt_plot), "\n", sep = "")
  cat("  Up (FDR<0.05, logFC>", logFC_cut_cells, "):   ", sum(tt_plot$sig_cat == "Up"), "\n", sep = "")
  cat("  Down (FDR<0.05, logFC<-", logFC_cut_cells, "): ", sum(tt_plot$sig_cat == "Down"), "\n", sep = "")
  cat("  NS:                          ", sum(tt_plot$sig_cat == "NS"), "\n", sep = "")

  ## Label genes: Top N per direction by FDR and |logFC|
  tt_sig <- tt_plot %>%
    dplyr::filter(adj.P.Val < fdr_cut_cells,
                  abs(logFC) >= logFC_cut_cells,
                  sig_cat != "NS")

  up_sig <- tt_sig %>% dplyr::filter(sig_cat == "Up")
  down_sig <- tt_sig %>% dplyr::filter(sig_cat == "Down")

  up_top_fdr <- up_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = volcano_label_n_fdr)
  up_top_fc  <- up_sig %>% dplyr::arrange(dplyr::desc(logFC)) %>% dplyr::slice_head(n = volcano_label_n_fc)

  down_top_fdr <- down_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = volcano_label_n_fdr)
  down_top_fc  <- down_sig %>% dplyr::arrange(logFC) %>% dplyr::slice_head(n = volcano_label_n_fc)

  label_genes <- dplyr::bind_rows(up_top_fdr, up_top_fc, down_top_fdr, down_top_fc) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)

  cat("Genes to label: ", nrow(label_genes), "\n", sep = "")

  ## =======================================================
  ## VOLCANO PLOT: orange for up, blue for down, grey for NS
  ## =======================================================

  vol_cells <- ggplot(
    tt_plot,
    aes(x = logFC, y = negLogFDR, color = sig_cat)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("Up" = "#E69F00", "Down" = "#0072B2", "NS" = "grey70"),
      name   = "Significance"
    ) +
    geom_point(data = label_genes, aes(color = sig_cat), size = 3) +
    geom_text_repel(
      data          = label_genes,
      aes(label     = gene_name),
      color         = "black",
      size          = 3.5,
      box.padding   = unit(0.25, "lines"),
      segment.color = "gray40",
      segment.size  = 0.3,
      point.padding = unit(0.2, "lines"),
      max.overlaps  = Inf
    ) +
    geom_hline(yintercept = -log10(fdr_cut_cells), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(logFC_cut_cells, -logFC_cut_cells), linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(size = 10, face = "bold", colour = "black"),
      axis.title      = element_text(size = 14, face = "bold", colour = "black"),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    ggtitle("Volcano: KAT8KD vs CTL (3T3-L1 Cells; DESeq2)") +
    xlab(expression(log[2]~Fold~Change)) +
    ylab(expression(-log[10]~FDR))

  volcano_file <- paste0("Volcano_cells_KAT8KD_vs_CTL_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", volcano_file), plot = vol_cells, width = 8, height = 6, dpi = 300)
  ## Save PDF version
  volcano_file_pdf <- paste0("Volcano_cells_KAT8KD_vs_CTL_", run_tag, ".pdf")
  ggsave(file.path(outdir, "plots", volcano_file_pdf), plot = vol_cells, width = 8, height = 6, device = "pdf")
  cat("Volcano saved: ", file.path(outdir, "plots", volcano_file), " (PNG and PDF)\n", sep = "")

  hero_volcano_file <- volcano_file

  ## =======================================================
  ## 4.14) ComplexHeatmap for top DEGs
  ## =======================================================
  cat("\n--- ComplexHeatmap: Top DEGs ---\n")

  de_sig <- tt %>%
    dplyr::filter(adj.P.Val < fdr_cut_cells, abs(logFC) > logFC_cut_cells)

  if (nrow(de_sig) >= 2) {
    de_sig <- de_sig[order(de_sig$adj.P.Val), , drop = FALSE]
    if (nrow(de_sig) > params$heatmap_max_genes) de_sig <- de_sig[1:params$heatmap_max_genes, , drop = FALSE]
    gene_set <- de_sig$gene_name

    cat("[INFO] Number of DEGs for heatmap: ", length(gene_set),
        " (threshold: FDR<", fdr_cut_cells, ", |logFC|>", logFC_cut_cells, ")\n", sep = "")

    ## Get sample metadata
    heat_meta <- as.data.frame(colData(dds_cells))

    ## Build matrix (using blind=FALSE VST for heatmaps)
    mat <- vst_mat_hm[gene_set, , drop = FALSE]

    ## Scale rows (z-score)
    mat_scaled <- t(scale(t(mat)))
    mat_scaled <- mat_scaled[!rowSums(is.na(mat_scaled)), , drop = FALSE]

    if (nrow(mat_scaled) >= 2) {
      ## Color scale (matching tissue analysis: blue-white-orange)
      max_val <- max(abs(mat_scaled), na.rm = TRUE)
      heat_col_fun <- colorRamp2(
        c(-max_val, 0, max_val),
        c("#0072B2", "white", "#E69F00")  ## Blue (down) -> white -> orange (up)
      )

      ## Column annotation
      heat_ha <- HeatmapAnnotation(
        Genotype = heat_meta$Genotype,
        col = list(Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")),
        annotation_name_side = "left"
      )

      ## Create heatmap
      heat_deg <- Heatmap(
        mat_scaled,
        col = heat_col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = (nrow(mat_scaled) <= 50),
        column_split = heat_meta$Genotype,
        border = TRUE,
        column_gap = unit(2, "mm"),
        row_gap = unit(0, "mm"),
        top_annotation = heat_ha,
        column_title = paste0("Top DEGs: KAT8KD vs CTL (FDR<", fdr_cut_cells, ", |logFC|>", logFC_cut_cells, ")"),
        heatmap_legend_param = list(
          title = "Z-score",
          title_position = "leftcenter-rot",
          title_gp = gpar(fontsize = 10)
        ),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8)
      )

      ## Save PNG
      heat_file <- paste0("Heatmap_cells_KAT8KD_vs_CTL_", run_tag, ".png")
      png(file.path(outdir, "plots", heat_file), width = 2000, height = 3000, res = 200)
      draw(heat_deg)
      dev.off()
      ## Save PDF version
      heat_file_pdf <- paste0("Heatmap_cells_KAT8KD_vs_CTL_", run_tag, ".pdf")
      pdf(file.path(outdir, "plots", heat_file_pdf), width = 10, height = 15)
      draw(heat_deg)
      dev.off()
      cat("[OK] ComplexHeatmap saved: ", heat_file, " (PNG and PDF)\n", sep = "")

    } else {
      cat("[WARN] Not enough valid genes for heatmap after scaling\n")
    }
  } else {
    cat("[WARN] Not enough DE genes for heatmap\n")
  }

  ## 4.15) save_core bookkeeping ----------------------------
  if (!is.null(hero_volcano_file)) {
    hero_path <- file.path(outdir, "plots", hero_volcano_file)
    if (file.exists(hero_path) && exists("save_run_file")) {
      tryCatch({
        save_run_file(hero_path, tag = "hero_volcano")
      }, error = function(e) {
        cat("[WARN] save_run_file failed: ", conditionMessage(e), "\n", sep = "")
      })
    }
  }

  if (exists("dedupe")) {
    try(dedupe(outdir), silent = TRUE)
  }

  ## Save session info for reproducibility
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_", run_tag, ".txt"))
  sink(session_file)
  cat("=== R SESSION INFORMATION ===\n\n")
  print(sessionInfo())
  cat("\n=== KEY PACKAGE VERSIONS ===\n\n")
  key_pkgs <- c("DESeq2", "dplyr", "ggplot2", "ComplexHeatmap")
  print(installed.packages()[intersect(key_pkgs, rownames(installed.packages())), c("Version", "Built")])
  sink()
  cat("[INFO] Session info saved to: ", session_file, "\n", sep = "")

  cat("\n======================================\n")
  cat("=== PART 1 COMPLETE ===\n")
  cat("======================================\n")
  cat("DE results saved. Run part2_cells_ora.R and part3_cells_fgsea.R for pathway analysis.\n\n")

}, error = function(e) {

  cat("\n======================================\n")
  cat("=== ERROR IN PIPELINE ===\n")
  cat("======================================\n")
  cat(conditionMessage(e), "\n\n")

  err_file <- paste0("ERROR_", run_tag, ".txt")
  writeLines(
    c(
      "Script error:",
      conditionMessage(e),
      "",
      "Traceback:",
      capture.output(traceback())
    ),
    con = file.path(outdir, "logs", err_file)
  )

  stop(e)
})

cat("\n=== Part 1 finished ===\n")
