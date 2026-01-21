#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - COMPREHENSIVE ANALYSIS (CELLS)
## =========================================================
## This script performs complete analysis in one run:
## - Part 1: DESeq2 differential expression + QC + volcano + heatmap
## - Part 2: ORA pathway analysis (GO:BP, KEGG)
## - Part 3: FGSEA pathway analysis (GO:BP, KEGG, WikiPathways, Hallmark)
##
## Cells: 3T3-L1 adipocytes, KAT8KD vs CTL
## Samples: JS_41-JS_48 (4 CTL, 4 KAT8KD)
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
  "scales",
  "msigdbr"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "DOSE",
  "AnnotationDbi",
  "ComplexHeatmap",
  "fgsea",
  "GO.db",
  "GOSemSim"
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
  library(grid)
  library(scales)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
  library(fgsea)
  library(GO.db)
  library(msigdbr)
})

## =======================================================
## UNIFIED PLOT THEME & COLORS
## =======================================================
plot_colors <- list(
  genotype = c(CTL = "#0072B2", KAT8KD = "#D55E00"),
  viridis_option = "viridis",
  ns_grey = "grey75",
  text_black = "black",
  line_grey = "grey40"
)

theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background = element_blank(),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "grey30"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.position = "right",
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = base_size - 1, face = "bold")
    )
}

volcano_label_n_fdr <- 10
volcano_label_n_fc <- 10

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
source(file.path(utils_dir, "save.R"))
source(file.path(utils_dir, "findsave.R"))
source(file.path(utils_dir, "purge.R"))
source(file.path(utils_dir, "dedupe.R"))

## 2.5) Run metadata ----------------------------------------
run_info <- list(
  script_name = "KAT8_bulk_cells_comprehensive.R",
  species     = "mouse",
  data_type   = "bulkRNAseq_3T3L1_cells",
  keywords    = c("KAT8", "3T3-L1", "adipocytes", "bulk RNA-seq", "cells", "comprehensive"),
  notes       = "Comprehensive analysis: DESeq2 + ORA + FGSEA in one script",
  message     = "Cells: Complete DESeq2 + ORA + FGSEA analysis"
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

## =======================================================
## ANALYSIS PARAMETERS
## =======================================================
## PREFILTER: Keep genes with sufficient expression in at least one group
prefilter_min_count <- 10                # Minimum count threshold
prefilter_min_samples_per_group <- 2     # Min samples per group (was 3, now 2 = 50% of 4 samples)

## DIFFERENTIAL EXPRESSION THRESHOLDS
## NOTE: Cell culture typically shows smaller fold changes than tissues
## Using 0.5 (1.4-fold) instead of 1.0 (2-fold) to capture biologically relevant changes
logFC_cut_cells <- 0.5   # |log2FC| threshold (was 1, changed to 0.5 for cells)
fdr_cut_cells   <- 0.05  # FDR threshold (standard)

## ORA settings
use_ora_universe <- FALSE

## FGSEA settings
run_go_bp <- TRUE
run_kegg <- TRUE
run_wikipathways <- TRUE
run_hallmark <- TRUE
simplify_go_bp <- TRUE
simplify_go_cutoff <- 0.7
top_n_per_direction <- 10

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## Sanity check trackers
pathway_sanity_tracker <- data.frame(
  Contrast = character(),
  Direction = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_Mapped_Entrez = integer(),
  N_Universe_Entrez = integer(),
  N_Sig_Pathways = integer(),
  stringsAsFactors = FALSE
)

fgsea_sanity_tracker <- data.frame(
  Contrast = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_Pathways_Tested = integer(),
  N_Sig_Pathways = integer(),
  N_Sig_Up = integer(),
  N_Sig_Down = integer(),
  stringsAsFactors = FALSE
)

## 4) Main logic --------------------------------------------
hero_volcano_file <- NULL

tryCatch({

  cat("\n")
  cat("================================================================================\n")
  cat("=== PART 1: DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ===\n")
  cat("================================================================================\n")
  log_time("Starting Part 1")

  ## 4.1) Load raw counts -----------------------------------
  cat("\n=== RAW COUNTS MATRIX (FULL) ===\n")

  counts_raw <- read.delim(
    "counts.txt",
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  cat("Dimensions (rows x cols):", paste(dim(counts_raw), collapse = " x "), "\n")

  ## 4.2) Define cell samples + annotation ------------------
  cell_samples <- paste0("JS_", sprintf("%02d", 41:48))

  sample_annot <- data.frame(
    Sample   = cell_samples,
    Genotype = c(rep("CTL", 4), rep("KAT8KD", 4)),
    stringsAsFactors = FALSE
  )

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

  if (any(duplicated(counts_cells$gene_name))) {
    dup_flag <- duplicated(counts_cells$gene_name) | duplicated(counts_cells$gene_name, fromLast = TRUE)
    counts_cells$gene_name[dup_flag] <- paste0(counts_cells$ensembl_gene_id[dup_flag], "_", counts_cells$gene_name[dup_flag])
  }
  counts_cells$gene_name <- make.unique(counts_cells$gene_name)

  rownames(counts_cells) <- counts_cells$gene_name

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

  ## 4.6) DESeq2 object + prefilter -------------------------
  min_count <- prefilter_min_count
  min_samples_per_group <- prefilter_min_samples_per_group

  ## IMPROVED PREFILTER: Require min samples in AT LEAST ONE GROUP (not total)
  ## This is more appropriate for comparing two groups
  ctl_samples_idx <- sample_annot$Sample[sample_annot$Genotype == "CTL"]
  kat8kd_samples_idx <- sample_annot$Sample[sample_annot$Genotype == "KAT8KD"]

  keep_ctl <- rowSums(count_matrix_cells[, ctl_samples_idx] >= min_count) >= min_samples_per_group
  keep_kat8kd <- rowSums(count_matrix_cells[, kat8kd_samples_idx] >= min_count) >= min_samples_per_group
  keep <- keep_ctl | keep_kat8kd

  n_before <- nrow(count_matrix_cells)
  n_after <- sum(keep)
  n_ctl_only <- sum(keep_ctl & !keep_kat8kd)
  n_kat8kd_only <- sum(!keep_ctl & keep_kat8kd)
  n_both <- sum(keep_ctl & keep_kat8kd)
  pct_retained <- round(100 * n_after / n_before, 1)

  cat("\n=== GENE PREFILTERING (PER-GROUP) ===\n")
  cat("Prefilter rule: counts >= ", min_count, " in >= ", min_samples_per_group, " samples in AT LEAST ONE group\n", sep = "")
  cat("Genes BEFORE: ", n_before, "\n", sep = "")
  cat("Genes AFTER:  ", n_after, " (", pct_retained, "%)\n", sep = "")
  cat("  - Expressed in CTL only:    ", n_ctl_only, "\n", sep = "")
  cat("  - Expressed in KAT8KD only: ", n_kat8kd_only, "\n", sep = "")
  cat("  - Expressed in both:        ", n_both, "\n", sep = "")

  count_matrix_filt <- count_matrix_cells[keep, , drop = FALSE]

  dds_cells <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filt),
    colData   = sample_annot,
    design    = ~ Genotype
  )

  ## Save parameters
  params <- list(
    prefilter_rule    = list(min_count = min_count, min_samples = min_samples_per_group),
    de_thresholds     = list(logFC_cutoff = logFC_cut_cells, fdr_cutoff = fdr_cut_cells),
    vst               = list(blind = TRUE),
    volcano_labels    = list(top_n_by_fdr = volcano_label_n_fdr, top_n_by_fc = volcano_label_n_fc),
    heatmap_max_genes = 300
  )

  params_file <- file.path(outdir, "logs", paste0("params_", run_tag, ".txt"))
  params_lines <- c(
    "=== RNA-seq Pipeline Parameters (Comprehensive Analysis, Cells) ===",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Prefilter Rule ---",
    paste0("  min_count: ", params$prefilter_rule$min_count),
    paste0("  min_samples: ", params$prefilter_rule$min_samples),
    "",
    "--- DE Thresholds ---",
    paste0("  |log2FC| cutoff: ", params$de_thresholds$logFC_cutoff),
    paste0("  FDR cutoff: ", params$de_thresholds$fdr_cutoff),
    "",
    "--- ORA Settings ---",
    paste0("  Use universe: ", use_ora_universe),
    "",
    "--- FGSEA Settings ---",
    paste0("  GO:BP simplification: ", simplify_go_bp),
    paste0("  Top N per direction: ", top_n_per_direction)
  )
  writeLines(params_lines, con = params_file)
  cat("Saved params to ", params_file, "\n")

  ## 4.7) Run DESeq2 ----------------------------------------
  cat("\n=== RUNNING DESeq2 ===\n")
  dds_cells <- DESeq(dds_cells)
  cat("[OK] DESeq2 complete.\n")

  ## DESeq2 sanity checks
  sf <- sizeFactors(dds_cells)
  cat("\nSize factors range: ", round(min(sf), 3), " - ", round(max(sf), 3), "\n", sep = "")

  ## Save dispersion plot
  disp_file <- paste0("Dispersion_plot_cells_", run_tag, ".png")
  png(file.path(outdir, "plots", disp_file), width = 1000, height = 800, res = 150)
  plotDispEsts(dds_cells, main = "DESeq2 Dispersion Estimates (Cells)")
  dev.off()
  cat("[OK] Saved dispersion plot\n")

  ## 4.8) VST transform for QC ------------------------------
  cat("\n=== VST FOR QC ===\n")
  vst_cells <- vst(dds_cells, blind = TRUE)
  vst_mat <- assay(vst_cells)

  vst_cells_hm <- vst(dds_cells, blind = FALSE)
  vst_mat_hm <- assay(vst_cells_hm)

  ## 4.9) PCA (VST) -----------------------------------------
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
    theme_publication(base_size = 12)

  pca_file <- paste0("PCA_cells_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", pca_file), plot = p_pca, width = 7, height = 6, dpi = 300)
  cat("Saved PCA plot\n")

  ## 4.10) MDS (VST) ----------------------------------------
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
    theme_publication(base_size = 12)

  mds_file <- paste0("MDS_cells_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", mds_file), plot = p_mds, width = 7, height = 6, dpi = 300)
  cat("Saved MDS plot\n")

  ## 4.11) Differential expression: KAT8KD vs CTL ----------
  cat("\n=== DIFFERENTIAL EXPRESSION: KAT8KD vs CTL ===\n")

  res <- results(dds_cells, contrast = c("Genotype", "KAT8KD", "CTL"))
  res <- res[order(res$pvalue), , drop = FALSE]

  tt <- as.data.frame(res)
  tt$logFC     <- tt$log2FoldChange
  tt$P.Value   <- tt$pvalue
  tt$adj.P.Val <- tt$padj
  tt$gene_name <- rownames(tt)
  tt <- tt[, !duplicated(colnames(tt)), drop = FALSE]

  de_file <- paste0("DE_cells_KAT8KD_vs_CTL_", run_tag, ".csv")
  write.csv(tt, file = file.path(outdir, "tables", de_file), row.names = TRUE)

  n_total <- nrow(tt)
  n_fdr_only <- sum(tt$adj.P.Val < fdr_cut_cells, na.rm = TRUE)
  n_up_both <- sum(tt$adj.P.Val < fdr_cut_cells & tt$logFC > logFC_cut_cells, na.rm = TRUE)
  n_down_both <- sum(tt$adj.P.Val < fdr_cut_cells & tt$logFC < -logFC_cut_cells, na.rm = TRUE)

  cat("\n=== DEG COUNTS ===\n")
  cat("Total genes tested: ", n_total, "\n", sep = "")
  cat("DEGs (FDR < ", fdr_cut_cells, " only): ", n_fdr_only, "\n", sep = "")
  cat("DEGs (FDR < ", fdr_cut_cells, " AND |logFC| > ", logFC_cut_cells, "): ", n_up_both + n_down_both, "\n", sep = "")
  cat("  - Up-regulated:   ", n_up_both, "\n", sep = "")
  cat("  - Down-regulated: ", n_down_both, "\n", sep = "")

  ## DIAGNOSTIC: Show fold change distribution for significant genes
  cat("\n=== FOLD CHANGE DIAGNOSTICS (FDR < ", fdr_cut_cells, ") ===\n", sep = "")
  sig_genes <- tt %>% dplyr::filter(adj.P.Val < fdr_cut_cells, !is.na(logFC))
  if (nrow(sig_genes) > 0) {
    cat("Among ", nrow(sig_genes), " FDR-significant genes:\n", sep = "")
    cat("  logFC range: ", round(min(sig_genes$logFC, na.rm = TRUE), 2), " to ",
        round(max(sig_genes$logFC, na.rm = TRUE), 2), "\n", sep = "")
    cat("  |logFC| quantiles:\n")
    abs_fc_quantiles <- quantile(abs(sig_genes$logFC), probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
    cat("    25%: ", round(abs_fc_quantiles[1], 2), "\n", sep = "")
    cat("    50% (median): ", round(abs_fc_quantiles[2], 2), "\n", sep = "")
    cat("    75%: ", round(abs_fc_quantiles[3], 2), "\n", sep = "")
    cat("    90%: ", round(abs_fc_quantiles[4], 2), "\n", sep = "")
    cat("    95%: ", round(abs_fc_quantiles[5], 2), "\n", sep = "")
    cat("  Genes with |logFC| > 0.25: ", sum(abs(sig_genes$logFC) > 0.25, na.rm = TRUE), "\n", sep = "")
    cat("  Genes with |logFC| > 0.5:  ", sum(abs(sig_genes$logFC) > 0.5, na.rm = TRUE), "\n", sep = "")
    cat("  Genes with |logFC| > 1.0:  ", sum(abs(sig_genes$logFC) > 1.0, na.rm = TRUE), "\n", sep = "")
    cat("  Genes with |logFC| > 2.0:  ", sum(abs(sig_genes$logFC) > 2.0, na.rm = TRUE), "\n", sep = "")
  } else {
    cat("No FDR-significant genes found!\n")
  }

  ## 4.12) Volcano plot -------------------------------------
  cat("\n--- VOLCANO PLOT ---\n")

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

  tt_sig <- tt_plot %>%
    dplyr::filter(adj.P.Val < fdr_cut_cells, abs(logFC) >= logFC_cut_cells, sig_cat != "NS")

  up_sig <- tt_sig %>% dplyr::filter(sig_cat == "Up")
  down_sig <- tt_sig %>% dplyr::filter(sig_cat == "Down")

  up_top_fdr <- up_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = volcano_label_n_fdr)
  up_top_fc  <- up_sig %>% dplyr::arrange(dplyr::desc(logFC)) %>% dplyr::slice_head(n = volcano_label_n_fc)

  down_top_fdr <- down_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = volcano_label_n_fdr)
  down_top_fc  <- down_sig %>% dplyr::arrange(logFC) %>% dplyr::slice_head(n = volcano_label_n_fc)

  label_genes <- dplyr::bind_rows(up_top_fdr, up_top_fc, down_top_fdr, down_top_fc) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE)

  vol_cells <- ggplot(tt_plot, aes(x = logFC, y = negLogFDR, color = sig_cat)) +
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
  cat("Volcano saved\n")

  hero_volcano_file <- volcano_file

  ## 4.13) ComplexHeatmap for top DEGs ---------------------
  cat("\n--- ComplexHeatmap: Top DEGs ---\n")

  de_sig <- tt %>%
    dplyr::filter(adj.P.Val < fdr_cut_cells, abs(logFC) > logFC_cut_cells)

  if (nrow(de_sig) >= 2) {
    de_sig <- de_sig[order(de_sig$adj.P.Val), , drop = FALSE]
    if (nrow(de_sig) > params$heatmap_max_genes) de_sig <- de_sig[1:params$heatmap_max_genes, , drop = FALSE]
    gene_set <- de_sig$gene_name

    heat_meta <- as.data.frame(colData(dds_cells))
    mat <- vst_mat_hm[gene_set, , drop = FALSE]
    mat_scaled <- t(scale(t(mat)))
    mat_scaled <- mat_scaled[!rowSums(is.na(mat_scaled)), , drop = FALSE]

    if (nrow(mat_scaled) >= 2) {
      max_val <- max(abs(mat_scaled), na.rm = TRUE)
      heat_col_fun <- colorRamp2(c(-max_val, 0, max_val), c("#0072B2", "white", "#E69F00"))

      heat_ha <- HeatmapAnnotation(
        Genotype = heat_meta$Genotype,
        col = list(Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")),
        annotation_name_side = "left"
      )

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
        top_annotation = heat_ha,
        column_title = paste0("Top DEGs: KAT8KD vs CTL (FDR<", fdr_cut_cells, ", |logFC|>", logFC_cut_cells, ")"),
        heatmap_legend_param = list(title = "Z-score", title_position = "leftcenter-rot")
      )

      heat_file <- paste0("Heatmap_cells_KAT8KD_vs_CTL_", run_tag, ".png")
      png(file.path(outdir, "plots", heat_file), width = 2000, height = 3000, res = 200)
      draw(heat_deg)
      dev.off()
      cat("[OK] ComplexHeatmap saved\n")
    }
  }

  log_time("Completed Part 1")

  ## ====================================================================
  ## PART 2: ORA PATHWAY ANALYSIS
  ## ====================================================================

  cat("\n")
  cat("================================================================================\n")
  cat("=== PART 2: ORA PATHWAY ANALYSIS ===\n")
  cat("================================================================================\n")
  log_time("Starting Part 2 (ORA)")

  ## Create ORA universe if requested
  universe_entrez <- NULL
  if (use_ora_universe) {
    cat("\n=== CREATING ORA UNIVERSE ===\n")
    universe_symbols <- rownames(tt)
    universe_entrez_df <- tryCatch({
      bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
    }, error = function(e) {
      cat("[ERROR] Universe gene ID conversion failed\n")
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    })
    universe_entrez <- unique(universe_entrez_df$ENTREZID)
    universe_entrez <- universe_entrez[!is.na(universe_entrez)]
    cat("[INFO] Universe Entrez IDs: ", length(universe_entrez), "\n", sep = "")
  } else {
    cat("[INFO] NOT using background gene set (default/liberal approach)\n")
  }

  ## ORA function
  run_ora_analysis <- function(gene_list, direction, contrast_name, universe_entrez, run_tag, outdir, fdr_cutoff = 0.05) {

    cat("\n--- ORA: ", direction, " genes (", contrast_name, ") ---\n", sep = "")

    if (length(gene_list) < 5) {
      cat("[WARN] Too few genes (", length(gene_list), ")\n", sep = "")
      return(NULL)
    }

    gene_entrez <- tryCatch({
      bitr(unique(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
    }, error = function(e) {
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    })

    gene_entrez <- gene_entrez[!is.na(gene_entrez$ENTREZID), , drop = FALSE]
    gene_entrez <- gene_entrez[!duplicated(gene_entrez$ENTREZID), , drop = FALSE]

    cat("[OK] Mapped ", nrow(gene_entrez), " of ", length(unique(gene_list)), " genes to Entrez\n", sep = "")

    if (nrow(gene_entrez) < 3) {
      cat("[WARN] Too few genes mapped\n")
      return(NULL)
    }

    entrez_ids <- gene_entrez$ENTREZID
    results <- list()

    ## GO:BP
    tryCatch({
      enrich_go <- if (!is.null(universe_entrez)) {
        enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.1,
                 qvalueCutoff = 0.2, readable = TRUE, minGSSize = 5, maxGSSize = 500, universe = universe_entrez)
      } else {
        enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.1,
                 qvalueCutoff = 0.2, readable = TRUE, minGSSize = 5, maxGSSize = 500)
      }
      if (!is.null(enrich_go) && nrow(enrich_go@result) > 0) {
        enrich_go_simp <- simplify(enrich_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
        results$gobp <- enrich_go_simp
        cat("[OK] GO:BP: ", nrow(enrich_go_simp@result), " terms\n", sep = "")
      }
    }, error = function(e) cat("[WARN] GO:BP failed\n"))

    ## KEGG
    tryCatch({
      enrich_kegg <- if (!is.null(universe_entrez)) {
        enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 0.1, qvalueCutoff = 0.2,
                   minGSSize = 5, maxGSSize = 500, universe = universe_entrez)
      } else {
        enrichKEGG(gene = entrez_ids, organism = "mmu", pvalueCutoff = 0.1, qvalueCutoff = 0.2,
                   minGSSize = 5, maxGSSize = 500)
      }
      if (!is.null(enrich_kegg) && nrow(enrich_kegg@result) > 0) {
        results$kegg <- enrich_kegg
        cat("[OK] KEGG: ", nrow(enrich_kegg@result), " pathways\n", sep = "")
      }
    }, error = function(e) cat("[WARN] KEGG failed\n"))

    ## Save results
    if (length(results) > 0) {
      for (db_name in names(results)) {
        enrich_obj <- results[[db_name]]
        sig_results <- enrich_obj@result %>%
          dplyr::filter(p.adjust < fdr_cutoff) %>%
          dplyr::arrange(p.adjust)

        if (nrow(sig_results) == 0) next

        sig_results <- sig_results[, !duplicated(colnames(sig_results)), drop = FALSE]

        table_file <- paste0("ORA_", db_name, "_KAT8KD_vs_CTL_", direction, "_", run_tag, ".csv")
        write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
        cat("[OK] Saved: ", table_file, "\n", sep = "")

        ## Bar plot
        plot_data <- sig_results %>% dplyr::slice_head(n = 15)
        plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
          num <- as.numeric(x[1]); denom <- as.numeric(x[2])
          if (is.na(denom) || denom == 0) return(0)
          num / denom
        })
        plot_data <- plot_data %>%
          dplyr::arrange(GeneRatio_numeric) %>%
          dplyr::mutate(Description = factor(Description, levels = Description))

        fill_scale <- if (direction == "Up") {
          scale_fill_gradient(low = "#FDD49E", high = "#E69F00", name = "Adj.\nP-value")
        } else {
          scale_fill_gradient(low = "#9ECAE1", high = "#0072B2", name = "Adj.\nP-value")
        }

        p_pathway <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description, fill = p.adjust)) +
          geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
          fill_scale +
          labs(title = paste0(toupper(db_name), " (ORA - ", direction, ")"), x = "Gene Ratio", y = NULL) +
          theme_classic(base_size = 12) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 10, color = "black"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
          )

        plot_file <- paste0("ORA_barplot_", db_name, "_KAT8KD_vs_CTL_", direction, "_", run_tag, ".png")
        ggsave(file.path(outdir, "plots", plot_file), plot = p_pathway, width = 10,
               height = max(6, nrow(plot_data) * 0.3), dpi = 300, bg = "white")
        cat("[OK] Saved pathway bar plot\n")
      }
    }

    return(results)
  }

  ## Get up/down genes
  up_genes <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_cells, logFC > logFC_cut_cells) %>%
    rownames()

  down_genes <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_cells, logFC < -logFC_cut_cells) %>%
    rownames()

  cat("[INFO] Up genes: ", length(up_genes), "\n", sep = "")
  cat("[INFO] Down genes: ", length(down_genes), "\n", sep = "")

  ## Run ORA
  if (length(up_genes) >= 5) {
    run_ora_analysis(up_genes, "Up", "KAT8KD_vs_CTL", universe_entrez, run_tag, outdir, fdr_cut_cells)
  }

  if (length(down_genes) >= 5) {
    run_ora_analysis(down_genes, "Down", "KAT8KD_vs_CTL", universe_entrez, run_tag, outdir, fdr_cut_cells)
  }

  log_time("Completed Part 2 (ORA)")

  ## ====================================================================
  ## PART 3: FGSEA PATHWAY ANALYSIS
  ## ====================================================================

  cat("\n")
  cat("================================================================================\n")
  cat("=== PART 3: FGSEA PATHWAY ANALYSIS ===\n")
  cat("================================================================================\n")
  log_time("Starting Part 3 (FGSEA)")

  ## Build gene sets
  cat("\n=== BUILDING GENE SETS ===\n")

  ## GO:BP
  go_bp_list <- NULL
  if (run_go_bp) {
    go_bp_genes <- tryCatch({
      AnnotationDbi::select(org.Mm.eg.db, keys = keys(org.Mm.eg.db, keytype = "GOALL"),
                           columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"), keytype = "GOALL")
    }, error = function(e) NULL)

    if (!is.null(go_bp_genes)) {
      go_bp_genes <- go_bp_genes %>% dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL))
      go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)
      go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= 5 & sapply(go_bp_list, length) <= 500]
      cat("[OK] GO:BP: ", length(go_bp_list), " pathways\n", sep = "")
    }
  }

  ## KEGG
  kegg_list <- NULL
  if (run_kegg) {
    tryCatch({
      msigdb_c2 <- msigdbr(species = "Mus musculus", collection = "C2")
      kegg_msigdb <- msigdb_c2 %>% dplyr::filter(grepl("^KEGG_", gs_name))
      if (nrow(kegg_msigdb) > 0) {
        kegg_list <- split(kegg_msigdb$gene_symbol, kegg_msigdb$gs_name)
        kegg_list <- kegg_list[sapply(kegg_list, length) >= 5 & sapply(kegg_list, length) <= 500]
        cat("[OK] KEGG: ", length(kegg_list), " pathways\n", sep = "")
      }
    }, error = function(e) cat("[WARN] KEGG failed\n"))
  }

  ## WikiPathways
  wikipathways_list <- NULL
  if (run_wikipathways) {
    tryCatch({
      msigdb_wp <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
      if (nrow(msigdb_wp) > 0) {
        wikipathways_list <- split(msigdb_wp$gene_symbol, msigdb_wp$gs_name)
        wikipathways_list <- wikipathways_list[sapply(wikipathways_list, length) >= 5 & sapply(wikipathways_list, length) <= 500]
        cat("[OK] WikiPathways: ", length(wikipathways_list), " pathways\n", sep = "")
      }
    }, error = function(e) cat("[WARN] WikiPathways failed\n"))
  }

  ## Hallmark
  hallmark_list <- NULL
  if (run_hallmark) {
    tryCatch({
      msigdb_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
      if (nrow(msigdb_hallmark) > 0) {
        hallmark_list <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)
        hallmark_list <- hallmark_list[sapply(hallmark_list, length) >= 5 & sapply(hallmark_list, length) <= 500]
        cat("[OK] Hallmark: ", length(hallmark_list), " pathways\n", sep = "")
      }
    }, error = function(e) cat("[WARN] Hallmark failed\n"))
  }

  ## Create ranked gene list
  ranked_genes <- tt %>% dplyr::filter(is.finite(stat)) %>% dplyr::arrange(dplyr::desc(stat))
  ranked_vec <- setNames(ranked_genes$stat, rownames(ranked_genes))
  cat("[INFO] Ranked genes: ", length(ranked_vec), "\n", sep = "")

  ## Run FGSEA for each database
  gene_sets <- list(gobp = go_bp_list, kegg = kegg_list, wikipathways = wikipathways_list, hallmark = hallmark_list)

  for (db_name in names(gene_sets)) {
    pathways <- gene_sets[[db_name]]
    if (is.null(pathways) || length(pathways) == 0) next

    cat("\n--- FGSEA: ", toupper(db_name), " ---\n", sep = "")

    tryCatch({
      fgsea_res <- fgsea(pathways = pathways, stats = ranked_vec, minSize = 5, maxSize = 500, nPermSimple = 10000)

      if (!is.null(fgsea_res) && nrow(fgsea_res) > 0) {
        ## Add descriptions
        if (db_name == "gobp") {
          go_terms <- AnnotationDbi::select(GO.db, keys = fgsea_res$pathway, columns = "TERM", keytype = "GOID")
          fgsea_res <- fgsea_res %>% dplyr::left_join(go_terms %>% dplyr::rename(pathway = GOID, Description = TERM), by = "pathway")
        } else {
          fgsea_res <- fgsea_res %>% dplyr::mutate(Description = gsub("^KEGG_|^HALLMARK_|^WIKIPATHWAYS_|^WP|^mmu", "", pathway), Description = gsub("_", " ", Description))
        }

        fgsea_res <- fgsea_res %>% dplyr::arrange(padj)

        sig_results <- fgsea_res %>% dplyr::filter(padj < fdr_cut_cells)
        cat("[OK] ", nrow(sig_results), " significant pathways\n", sep = "")

        if (nrow(sig_results) > 0) {
          ## Convert list columns to character
          for (col_name in colnames(sig_results)) {
            if (is.list(sig_results[[col_name]])) {
              sig_results[[col_name]] <- sapply(sig_results[[col_name]], function(x) {
                if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ";")
              })
            }
          }

          table_file <- paste0("fgsea_", db_name, "_KAT8KD_vs_CTL_", run_tag, ".csv")
          write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)

          ## Leading edge
          if ("leadingEdge" %in% colnames(sig_results)) {
            leading_edge_df <- sig_results %>%
              dplyr::select(pathway, Description, NES, padj, leadingEdge) %>%
              dplyr::mutate(Direction = ifelse(NES > 0, "Up", "Down"), N_LeadingEdge = sapply(strsplit(leadingEdge, ";"), length))

            leading_edge_file <- paste0("fgsea_leadingEdge_", db_name, "_KAT8KD_vs_CTL_", run_tag, ".csv")
            write.csv(leading_edge_df, file = file.path(outdir, "tables", leading_edge_file), row.names = FALSE)
          }

          ## Plots
          top_up <- sig_results %>% dplyr::filter(NES > 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n_per_direction)
          top_down <- sig_results %>% dplyr::filter(NES < 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n_per_direction)

          plot_data <- dplyr::bind_rows(top_down, top_up) %>%
            dplyr::arrange(NES) %>%
            dplyr::mutate(
              pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
              pathway_label = factor(pathway_label, levels = pathway_label),
              Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")
            )

          if (nrow(plot_data) > 0) {
            p_fgsea <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = Direction)) +
              geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
              scale_fill_manual(values = c("Up-regulated" = "#E69F00", "Down-regulated" = "#0072B2"), name = "Direction") +
              labs(title = paste0("fgsea ", toupper(db_name)), x = "Normalized Enrichment Score (NES)", y = NULL) +
              theme_classic(base_size = 12) +
              theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
              geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

            plot_file <- paste0("fgsea_plot_", db_name, "_KAT8KD_vs_CTL_", run_tag, ".png")
            ggsave(file.path(outdir, "plots", plot_file), plot = p_fgsea, width = 12, height = max(6, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
            cat("[OK] Saved fgsea bar plot\n")
          }
        }
      }
    }, error = function(e) cat("[WARN] FGSEA ", db_name, " failed: ", conditionMessage(e), "\n", sep = ""))
  }

  log_time("Completed Part 3 (FGSEA)")

  ## ====================================================================
  ## FINALIZE
  ## ====================================================================

  ## Session info
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_", run_tag, ".txt"))
  sink(session_file)
  cat("=== R SESSION INFORMATION ===\n\n")
  print(sessionInfo())
  sink()

  ## Cleanup
  if (exists("dedupe")) try(dedupe(outdir), silent = TRUE)

  cat("\n")
  cat("================================================================================\n")
  cat("=== COMPREHENSIVE ANALYSIS COMPLETE ===\n")
  cat("================================================================================\n")
  cat("Output in: ", normalizePath(outdir), "\n", sep = "")
  cat("\n")

}, error = function(e) {
  cat("\n[ERROR] ", conditionMessage(e), "\n\n", sep = "")
  err_file <- paste0("ERROR_", run_tag, ".txt")
  writeLines(c("Script error:", conditionMessage(e), "", "Traceback:", capture.output(traceback())),
             con = file.path(outdir, "logs", err_file))
  stop(e)
})

cat("=== Analysis finished ===\n")
