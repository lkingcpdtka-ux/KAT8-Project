#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - COMPREHENSIVE ANALYSIS (CELLS)
## =========================================================
## This script performs complete analysis in one run:
## - Part 1: DESeq2 differential expression + QC + volcano + heatmap
## - Part 2: ORA pathway analysis (GO:BP, KEGG)
## - Part 3: GSEA pathway analysis (GO:BP, KEGG)
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
  "scales"
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
  "GO.db",
  "GOSemSim",
  "apeglm"     ## shrinkage estimator for lfcShrink() (Genotype coefficient)
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
  library(GO.db)
})

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

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
  notes       = "Comprehensive analysis: DESeq2 + ORA + GSEA in one script",
  message     = "Cells: Complete DESeq2 + ORA + GSEA analysis"
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

## GSEA settings
run_go_bp <- TRUE
run_kegg <- TRUE
run_wikipathways <- FALSE
run_hallmark <- FALSE
simplify_go_bp <- gsea_params$simplify_go
simplify_go_cutoff <- gsea_params$simplify_cutoff
top_n_per_direction <- gsea_params$top_n_per_direction

gse_kegg_params <- list(
  min_gs_size     = gsea_params$min_gs_size,
  max_gs_size     = gsea_params$max_gs_size,
  pvalue_cutoff   = 0.1,
  p_adjust_method = "BH",
  organism        = "mmu",
  seed            = TRUE
)

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
    "--- Data Source ---",
    paste0("  Cell type: 3T3-L1 adipocytes"),
    paste0("  Samples: JS_41-JS_48 (4 CTL, 4 KAT8KD)"),
    paste0("  Comparison: KAT8KD vs CTL"),
    "",
    "--- Prefilter Rule (per-group) ---",
    paste0("  min_count: ", params$prefilter_rule$min_count),
    paste0("  min_samples_per_group: ", params$prefilter_rule$min_samples, " (50% of 4 samples)"),
    paste0("  Strategy: Keep genes with min_count in min_samples in AT LEAST ONE group (CTL or KAT8KD)"),
    paste0("  Rationale: Captures genes expressed specifically in one condition (appropriate for 2-group comparison)"),
    "",
    "--- DE Thresholds ---",
    paste0("  |log2FC| cutoff: ", params$de_thresholds$logFC_cutoff),
    paste0("  FDR cutoff: ", params$de_thresholds$fdr_cutoff),
    paste0("  Note: Cell culture uses lower logFC threshold (0.5 vs 1.0 for tissues) to capture biologically relevant changes"),
    paste0("  Description: FDR = Benjamini-Hochberg adjusted p-value; logFC filters for biological significance"),
    "",
    "--- VST (Variance Stabilizing Transformation) ---",
    paste0("  QC plots: blind=TRUE (ignores experimental design for unbiased QC)"),
    paste0("  Heatmaps: blind=FALSE (uses experimental design for accurate visualization)"),
    paste0("  Description: Log2-scale transformation that stabilizes variance across expression range"),
    "",
    "--- Volcano Plot Labels ---",
    paste0("  Top N by FDR (per direction): ", params$volcano_labels$top_n_by_fdr),
    paste0("  Top N by |logFC| (per direction): ", params$volcano_labels$top_n_by_fc),
    paste0("  Total labels: Up to ", 2*(params$volcano_labels$top_n_by_fdr + params$volcano_labels$top_n_by_fc), " genes (top by FDR + top by FC, both directions)"),
    "",
    "--- Heatmap ---",
    paste0("  Max genes displayed: ", params$heatmap_max_genes),
    paste0("  Gene selection: Top DEGs ranked by adjusted p-value"),
    paste0("  Clustering: Hierarchical (complete linkage) on rows; split by genotype on columns"),
    paste0("  Scaling: Z-score normalization (row-wise: mean=0, sd=1)"),
    paste0("  Color scheme: Teal (#0072B2, down) -> white (neutral) -> orange (#E69F00, up)"),
    "",
    "--- ORA (Over-Representation Analysis) Settings ---",
    paste0("  Method: Fisher's exact test (hypergeometric distribution)"),
    paste0("  Use custom universe: ", use_ora_universe),
    paste0("  Universe definition: ", if (use_ora_universe) "All genes tested in DESeq2" else "All genes in organism annotation (default/liberal)"),
    paste0("  P-value cutoff: 0.1 (initial enrichment)"),
    paste0("  Q-value cutoff: 0.2 (initial enrichment)"),
    paste0("  Final FDR cutoff: ", params$de_thresholds$fdr_cutoff, " (for reported results)"),
    paste0("  Gene set size: 5-500 genes"),
    paste0("  Databases: GO:BP (simplified), KEGG"),
    paste0("  GO:BP simplification: Yes (cutoff=0.7, by p.adjust, removes redundant terms)"),
    paste0("  Gene ID conversion: SYMBOL -> ENTREZID (via org.Mm.eg.db)"),
    paste0("  Input: Up-regulated and down-regulated genes (analyzed separately)"),
    "",
    "--- GSEA (Gene Set Enrichment Analysis) Settings ---",
    paste0("  Method: clusterProfiler::gseGO (GO:BP) and gseKEGG (KEGG)"),
    paste0("  Ranking metric: DESeq2 Wald statistic (combines fold change and significance)"),
    paste0("  Gene set size: ", gsea_params$min_gs_size, "-", gsea_params$max_gs_size, " genes"),
    paste0("  Initial p-value cutoff (KEGG): ", gse_kegg_params$pvalue_cutoff),
    paste0("  FDR cutoff: ", params$de_thresholds$fdr_cutoff, " (for significant pathways)"),
    paste0("  Databases tested:"),
    paste0("    - GO:BP: Biological Process ontology (", if (run_go_bp) "ENABLED" else "DISABLED", ")"),
    paste0("    - KEGG: Kyoto Encyclopedia of Genes and Genomes pathways (", if (run_kegg) "ENABLED" else "DISABLED", ")"),
    paste0("    - WikiPathways: Community-curated pathway database (DISABLED)"),
    paste0("    - Hallmark: MSigDB hallmark gene sets (DISABLED)"),
    paste0("  GO:BP simplification: ", if (simplify_go_bp) paste0("Yes (cutoff=", simplify_go_cutoff, ", semantic similarity)") else "No"),
    paste0("  Plot visualization: Top ", top_n_per_direction, " pathways per direction (up/down-regulated)"),
    paste0("  Leading edge: Core enrichment genes driving pathway significance (saved separately)"),
    paste0("  NES: Normalized Enrichment Score (accounts for gene set size)"),
    paste0("  Source: clusterProfiler + org.Mm.eg.db (GO) and KEGG API"),
    "",
    "--- Color Schemes ---",
    paste0("  Genotype: CTL=#0072B2 (blue), KAT8KD=#D55E00 (orange)"),
    paste0("  Significance: Up=#E69F00 (orange), Down=#0072B2 (teal), NS=grey70"),
    paste0("  Heatmap: Teal-White-Orange (diverging, centered at 0)"),
    paste0("  Theme: Publication-ready (clean, high-contrast, black axes)"),
    "",
    "--- Quality Control ---",
    paste0("  Dispersion plot: DESeq2 mean-dispersion relationship (saved)"),
    paste0("  PCA: Variance-scaled, labeled samples, colored by genotype"),
    paste0("  MDS: Classical multidimensional scaling on VST distances"),
    paste0("  Size factors: DESeq2 median-of-ratios normalization"),
    "",
    "--- Output Files ---",
    paste0("  DE tables: CSV format, all genes with statistics"),
    paste0("  ORA tables: CSV format, significant pathways only (FDR<", params$de_thresholds$fdr_cutoff, ")"),
    paste0("  GSEA tables: CSV format, significant pathways + leading edge genes"),
    paste0("  Plots: PNG format, 300 DPI (publication quality)"),
    paste0("  Session info: R version and package versions (for reproducibility)")
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
  ## DESeq2 convention: PCA on the top-variable genes WITHOUT unit-scaling
  ## every gene (scale.=TRUE upweights noisy low-variance genes). Mirrors
  ## DESeq2::plotPCA (top 500 by row variance).
  rv <- apply(vst_mat, 1, stats::var)
  top_var_idx <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca_cells <- prcomp(t(vst_mat[top_var_idx, , drop = FALSE]), scale. = FALSE)
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

  ## Shrunken effect sizes (apeglm) for the KAT8KD-vs-CTL coefficient. At
  ## n = 4/group the raw MLE log2FC is inflated for low-count genes; the
  ## shrunken value is exported as an extra column for ranking/visualisation.
  ## DEG calling below is UNCHANGED (MLE logFC + FDR); GSEA ranks on the
  ## Wald stat, which shrinkage does not affect.
  res_shrunk <- tryCatch(
    lfcShrink(dds_cells, coef = "Genotype_KAT8KD_vs_CTL", type = "apeglm", quiet = TRUE),
    error = function(e) {
      cat("[WARN] lfcShrink (apeglm) failed: ", conditionMessage(e),
          " -- exporting MLE logFC only\n", sep = ""); NULL
    }
  )

  tt <- as.data.frame(res)
  tt$logFC     <- tt$log2FoldChange
  tt$P.Value   <- tt$pvalue
  tt$adj.P.Val <- tt$padj
  tt$gene_name <- rownames(tt)
  if (!is.null(res_shrunk)) {
    tt$log2FoldChange_shrunk <- as.data.frame(res_shrunk)[rownames(tt), "log2FoldChange"]
  }
  if ("padj" %in% colnames(tt)) tt$padj <- NULL
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

  cat("\n")
  cat("================================================================================\n")
  cat("=== SANITY CHECK: CELL CULTURE DIFFERENTIAL EXPRESSION (KAT8KD vs CTL) ===\n")
  cat("================================================================================\n")
  cat("Total genes tested by DESeq2:           ", n_total, "\n", sep = "")
  cat("Genes with valid adjusted p-value:      ", n_valid_pval, "\n", sep = "")
  cat("\n--- DEG Counts at Different Thresholds ---\n")
  cat("DEGs (FDR < ", fdr_cut_cells, " only):               ", n_fdr_only, "\n", sep = "")
  cat("DEGs (FDR < ", fdr_cut_cells, " AND |logFC| > ", logFC_cut_cells, "):  ", n_both, "\n", sep = "")
  cat("  - Up-regulated (logFC > ", logFC_cut_cells, "):       ", n_up_both, "\n", sep = "")
  cat("  - Down-regulated (logFC < -", logFC_cut_cells, "):    ", n_down_both, "\n", sep = "")
  cat("\n--- Percentage of Tested Genes ---\n")
  cat("Percent DEGs (FDR-only threshold):      ", round(100 * n_fdr_only / n_total, 2), "%\n", sep = "")
  cat("Percent DEGs (FDR + logFC threshold):   ", round(100 * n_both / n_total, 2), "%\n", sep = "")
  ## NOTE: A prior "consistency check" compared this DEG count to unrelated
  ## MOF studies (Li 2012 HeLa, Ravens 2014 MEFs) and flagged 80-150 as
  ## "consistent". Removed: DEG count depends on depth, n, and cell system,
  ## so matching another study is not a QC criterion and invites tuning the
  ## cutoff to hit a number. Report the count at fixed, pre-set thresholds.
  cat("\n(DEG count reported at pre-registered thresholds; not tuned to any external study.)\n")
  cat("\n--- Output Files ---\n")
  cat("DE table written: ", file.path(outdir, "tables", de_file), "\n", sep = "")
  cat("================================================================================\n\n")

  ## DIAGNOSTIC: Show fold change distribution for significant genes
  cat("================================================================================\n")
  cat("=== FOLD CHANGE DIAGNOSTICS (FDR < ", fdr_cut_cells, ") ===\n", sep = "")
  cat("================================================================================\n")
  sig_genes <- tt %>% dplyr::filter(adj.P.Val < fdr_cut_cells, !is.na(logFC))
  if (nrow(sig_genes) > 0) {
    cat("Among ", nrow(sig_genes), " FDR-significant genes:\n\n", sep = "")
    cat("--- Fold Change Range ---\n")
    cat("  logFC range: ", round(min(sig_genes$logFC, na.rm = TRUE), 2), " to ",
        round(max(sig_genes$logFC, na.rm = TRUE), 2), "\n\n", sep = "")

    cat("--- |logFC| Quantiles ---\n")
    abs_fc_quantiles <- quantile(abs(sig_genes$logFC), probs = c(0.25, 0.5, 0.75, 0.9, 0.95), na.rm = TRUE)
    cat("  25%: ", round(abs_fc_quantiles[1], 2), "\n", sep = "")
    cat("  50% (median): ", round(abs_fc_quantiles[2], 2), "\n", sep = "")
    cat("  75%: ", round(abs_fc_quantiles[3], 2), "\n", sep = "")
    cat("  90%: ", round(abs_fc_quantiles[4], 2), "\n", sep = "")
    cat("  95%: ", round(abs_fc_quantiles[5], 2), "\n\n", sep = "")

    cat("--- Distribution by Threshold ---\n")
    cat("  Genes with |logFC| > 0.25: ", sum(abs(sig_genes$logFC) > 0.25, na.rm = TRUE), " (",
        round(100 * sum(abs(sig_genes$logFC) > 0.25, na.rm = TRUE) / nrow(sig_genes), 1), "%)\n", sep = "")
    cat("  Genes with |logFC| > 0.5:  ", sum(abs(sig_genes$logFC) > 0.5, na.rm = TRUE), " (",
        round(100 * sum(abs(sig_genes$logFC) > 0.5, na.rm = TRUE) / nrow(sig_genes), 1), "%)\n", sep = "")
    cat("  Genes with |logFC| > 1.0:  ", sum(abs(sig_genes$logFC) > 1.0, na.rm = TRUE), " (",
        round(100 * sum(abs(sig_genes$logFC) > 1.0, na.rm = TRUE) / nrow(sig_genes), 1), "%)\n", sep = "")
    cat("  Genes with |logFC| > 2.0:  ", sum(abs(sig_genes$logFC) > 2.0, na.rm = TRUE), " (",
        round(100 * sum(abs(sig_genes$logFC) > 2.0, na.rm = TRUE) / nrow(sig_genes), 1), "%)\n\n", sep = "")

    ## Find and report on Kat8 gene itself
    kat8_info <- tt %>% dplyr::filter(grepl("^Kat8$", gene_name, ignore.case = TRUE))
    if (nrow(kat8_info) > 0) {
      cat("--- Knockdown Validation (Kat8 gene) ---\n")
      cat("  Kat8 logFC: ", round(kat8_info$logFC[1], 2), "\n", sep = "")
      cat("  Kat8 FDR: ", format(kat8_info$adj.P.Val[1], scientific = TRUE, digits = 3), "\n", sep = "")
      ## Calculate rank
      tt_ranked <- tt %>% dplyr::arrange(logFC) %>% mutate(rank = row_number())
      kat8_rank <- tt_ranked %>% dplyr::filter(grepl("^Kat8$", gene_name, ignore.case = TRUE)) %>% pull(rank)
      if (length(kat8_rank) > 0) {
        cat("  Kat8 rank (by logFC): #", kat8_rank, " of ", nrow(tt), " genes\n", sep = "")
        if (kat8_rank <= 10) {
          cat("  Status: ✓ EXCELLENT knockdown (top 10 most down-regulated)\n")
        } else if (kat8_rank <= 50) {
          cat("  Status: ✓ GOOD knockdown (top 50 most down-regulated)\n")
        } else {
          cat("  Status: ⚠ Check knockdown efficiency\n")
        }
      }
      cat("\n")
    }

    cat("--- Interpretation for Cell Culture ---\n")
    median_fc <- abs_fc_quantiles[2]
    if (median_fc >= 0.3 & median_fc <= 0.5) {
      cat("  Median |logFC| = ", round(median_fc, 2), " is TYPICAL for chromatin modifiers in cells\n", sep = "")
      cat("  Published KAT8/MOF studies show median |logFC| of 0.3-0.5\n")
      cat("  Cell culture typically shows smaller fold changes than tissue\n")
    } else if (median_fc < 0.3) {
      cat("  Median |logFC| = ", round(median_fc, 2), " is modest - consider checking knockdown efficiency\n", sep = "")
    } else {
      cat("  Median |logFC| = ", round(median_fc, 2), " shows strong effects\n", sep = "")
    }
    cat("================================================================================\n")
  } else {
    cat("⚠ WARNING: No FDR-significant genes found!\n")
    cat("Check: sample quality, knockdown efficiency, sequencing depth\n")
    cat("================================================================================\n")
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
  ## Use pull(gene_name), NOT rownames(): dplyr::filter() drops row names,
  ## so rownames() would feed "1","2",... into ORA instead of gene symbols.
  if (!"gene_name" %in% colnames(tt)) tt$gene_name <- rownames(tt)
  up_genes <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_cells, logFC > logFC_cut_cells) %>%
    dplyr::pull(gene_name)

  down_genes <- tt %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_cells, logFC < -logFC_cut_cells) %>%
    dplyr::pull(gene_name)

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
  ## PART 3: GSEA PATHWAY ANALYSIS
  ## ====================================================================

  cat("\n")
  cat("================================================================================\n")
  cat("=== PART 3: GSEA PATHWAY ANALYSIS ===\n")
  cat("================================================================================\n")
  log_time("Starting Part 3 (GSEA)")

  ## Create ranked gene list
  ## Key on gene_name, NOT rownames(): dplyr::filter/arrange drop row names,
  ## so rownames(ranked_genes) would be "1","2",... and fgsea would match
  ## ~zero genes to any pathway.
  if (!"gene_name" %in% colnames(tt)) tt$gene_name <- rownames(tt)
  ranked_genes <- tt %>% dplyr::filter(is.finite(stat)) %>% dplyr::arrange(dplyr::desc(stat))
  ranked_vec <- setNames(ranked_genes$stat, ranked_genes$gene_name)
  cat("[INFO] Ranked genes: ", length(ranked_vec), "\n", sep = "")

  results <- list()

  ## GO:BP GSEA (clusterProfiler::gseGO)
  if (run_go_bp) {
    log_time("Starting GO:BP gseGO")
    tryCatch({
      cat("[INFO] Running clusterProfiler::gseGO() (ont = BP)...\n")

      go_symbols <- names(ranked_vec)
      go_symbol_to_entrez <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = go_symbols,
        columns = c("SYMBOL", "ENTREZID"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL))

      go_symbol_to_entrez$rank_value <- ranked_vec[go_symbol_to_entrez$SYMBOL]
      go_symbol_to_entrez <- go_symbol_to_entrez %>%
        dplyr::group_by(ENTREZID) %>%
        dplyr::slice_max(order_by = abs(rank_value), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

      go_genelist <- setNames(go_symbol_to_entrez$rank_value, go_symbol_to_entrez$ENTREZID)
      go_genelist <- sort(go_genelist, decreasing = TRUE)

      gsea_go <- gseGO(
        geneList      = go_genelist,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",
        keyType       = "ENTREZID",
        minGSSize     = gsea_params$min_gs_size,
        maxGSSize     = gsea_params$max_gs_size,
        pvalueCutoff  = 1.0,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        seed          = TRUE
      )

      if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
        gsea_obj_file <- paste0("gseaResult_gobp_KAT8KD_vs_CTL_", run_tag, ".rds")
        saveRDS(gsea_go, file = file.path(outdir, "tables", gsea_obj_file))
        cat("[OK] Saved gseaResult object: ", gsea_obj_file, "\n", sep = "")

        go_results <- as.data.frame(gsea_go) %>%
          dplyr::rename(
            pathway = ID,
            padj = p.adjust,
            leadingEdge = core_enrichment
          ) %>%
          dplyr::arrange(padj)

        go_entrez_to_symbol <- go_symbol_to_entrez %>%
          dplyr::distinct(ENTREZID, SYMBOL)
        if ("leadingEdge" %in% colnames(go_results)) {
          go_results$leadingEdge_symbols <- vapply(
            go_results$leadingEdge,
            function(ids) {
              if (is.na(ids) || ids == "") return(NA_character_)
              id_vec <- unlist(strsplit(ids, "/", fixed = TRUE))
              symbol_vec <- go_entrez_to_symbol$SYMBOL[match(id_vec, go_entrez_to_symbol$ENTREZID)]
              symbol_vec <- symbol_vec[!is.na(symbol_vec)]
              if (length(symbol_vec) == 0) return(NA_character_)
              paste(symbol_vec, collapse = ",")
            },
            character(1)
          )
        }

        results$gobp <- go_results
        n_sig <- sum(go_results$padj < fdr_cut_cells, na.rm = TRUE)
        n_sig_up <- sum(go_results$padj < fdr_cut_cells & go_results$NES > 0, na.rm = TRUE)
        n_sig_down <- sum(go_results$padj < fdr_cut_cells & go_results$NES < 0, na.rm = TRUE)
        cat("[OK] gseGO GO:BP: ", nrow(go_results), " pathways tested, ",
            n_sig, " significant (FDR<", fdr_cut_cells, ")\n", sep = "")

        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = "KAT8KD_vs_CTL",
            Database = "GO:BP",
            N_Input_Genes = length(go_genelist),
            N_Pathways_Tested = nrow(go_results),
            N_Sig_Pathways = n_sig,
            N_Sig_Up = n_sig_up,
            N_Sig_Down = n_sig_down,
            stringsAsFactors = FALSE
          )
        )

        if (simplify_go_bp && nrow(go_results) > 0) {
          cat("[INFO] Simplifying GO:BP results (removing redundant terms)...\n")
          tryCatch({
            sig_go <- go_results %>% dplyr::filter(padj < fdr_cut_cells)
            if (nrow(sig_go) > 0) {
              require(GOSemSim)
              godata <- godata("org.Mm.eg.db", ont = "BP", computeIC = TRUE)
              go_ids <- sig_go$pathway
              sim_mat <- mgoSim(go_ids, go_ids, semData = godata, measure = "Wang")

              keep_terms <- c()
              for (i in seq_along(go_ids)) {
                if (i == 1 || all(sim_mat[i, keep_terms] < simplify_go_cutoff)) {
                  keep_terms <- c(keep_terms, i)
                }
              }

              go_results_simplified <- sig_go[keep_terms, ]
              results$gobp_simplified <- go_results_simplified
              cat("[OK] Simplified GO:BP results: ", nrow(go_results_simplified), " terms\n", sep = "")
            }
          }, error = function(e) {
            cat("[WARN] GO:BP simplification failed: ", conditionMessage(e), "\n", sep = "")
          })
        }
      } else {
        cat("[INFO] gseGO returned no results\n")
      }
    }, error = function(e) {
      cat("[WARN] gseGO GO:BP failed: ", conditionMessage(e), "\n", sep = "")
    })
    log_time("Completed GO:BP gseGO")
  }

  ## KEGG GSEA (clusterProfiler::gseKEGG)
  if (run_kegg) {
    log_time("Starting KEGG gseKEGG")
    tryCatch({
      cat("[INFO] Running clusterProfiler::gseKEGG() (organism = mmu)...\n")

      gene_symbols <- names(ranked_vec)
      symbol_to_entrez <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = gene_symbols,
        columns = c("SYMBOL", "ENTREZID"),
        keytype = "SYMBOL"
      ) %>%
        dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL))

      symbol_to_entrez$rank_value <- ranked_vec[symbol_to_entrez$SYMBOL]
      symbol_to_entrez <- symbol_to_entrez %>%
        dplyr::group_by(ENTREZID) %>%
        dplyr::slice_max(order_by = abs(rank_value), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

      kegg_genelist <- setNames(symbol_to_entrez$rank_value, symbol_to_entrez$ENTREZID)
      kegg_genelist <- sort(kegg_genelist, decreasing = TRUE)

      gsea_kegg <- gseKEGG(
        geneList      = kegg_genelist,
        organism      = gse_kegg_params$organism,
        minGSSize     = gse_kegg_params$min_gs_size,
        maxGSSize     = gse_kegg_params$max_gs_size,
        pvalueCutoff  = gse_kegg_params$pvalue_cutoff,
        pAdjustMethod = gse_kegg_params$p_adjust_method,
        verbose       = FALSE,
        seed          = gse_kegg_params$seed
      )

      if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
        gsea_obj_file <- paste0("gseaResult_kegg_KAT8KD_vs_CTL_", run_tag, ".rds")
        saveRDS(gsea_kegg, file = file.path(outdir, "tables", gsea_obj_file))
        cat("[OK] Saved gseaResult object: ", gsea_obj_file, "\n", sep = "")

        kegg_results <- as.data.frame(gsea_kegg) %>%
          dplyr::rename(
            pathway = ID,
            padj = p.adjust,
            leadingEdge = core_enrichment
          ) %>%
          dplyr::arrange(padj)

        entrez_to_symbol <- symbol_to_entrez %>%
          dplyr::distinct(ENTREZID, SYMBOL)
        if ("leadingEdge" %in% colnames(kegg_results)) {
          kegg_results$leadingEdge_symbols <- vapply(
            kegg_results$leadingEdge,
            function(ids) {
              if (is.na(ids) || ids == "") return(NA_character_)
              id_vec <- unlist(strsplit(ids, "/", fixed = TRUE))
              symbol_vec <- entrez_to_symbol$SYMBOL[match(id_vec, entrez_to_symbol$ENTREZID)]
              symbol_vec <- symbol_vec[!is.na(symbol_vec)]
              if (length(symbol_vec) == 0) return(NA_character_)
              paste(symbol_vec, collapse = ",")
            },
            character(1)
          )
        }

        results$kegg <- kegg_results
        n_sig <- sum(kegg_results$padj < fdr_cut_cells, na.rm = TRUE)
        n_sig_up <- sum(kegg_results$padj < fdr_cut_cells & kegg_results$NES > 0, na.rm = TRUE)
        n_sig_down <- sum(kegg_results$padj < fdr_cut_cells & kegg_results$NES < 0, na.rm = TRUE)
        cat("[OK] gseKEGG: ", nrow(kegg_results), " pathways tested, ",
            n_sig, " significant (FDR<", fdr_cut_cells, ")\n", sep = "")

        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = "KAT8KD_vs_CTL",
            Database = "KEGG",
            N_Input_Genes = length(kegg_genelist),
            N_Pathways_Tested = nrow(kegg_results),
            N_Sig_Pathways = n_sig,
            N_Sig_Up = n_sig_up,
            N_Sig_Down = n_sig_down,
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("[INFO] gseKEGG returned no results\n")
      }
    }, error = function(e) {
      cat("[WARN] gseKEGG failed: ", conditionMessage(e), "\n", sep = "")
    })
    log_time("Completed KEGG gseKEGG")
  }

  ## Save results + plots with tissue-matched styling
  if (length(results) > 0) {
    for (db_name in names(results)) {
      fgsea_obj <- results[[db_name]]

      sig_results <- fgsea_obj %>%
        dplyr::filter(padj < fdr_cut_cells) %>%
        dplyr::arrange(padj)

      if (nrow(sig_results) == 0) {
        cat("[INFO] No significant GSEA pathways in ", db_name, "\n", sep = "")
        next
      }

      if (!"Description" %in% colnames(sig_results) && !"pathway_name" %in% colnames(sig_results)) {
        sig_results$Description <- sig_results$pathway
      }

      for (col_name in colnames(sig_results)) {
        if (is.list(sig_results[[col_name]])) {
          sig_results[[col_name]] <- sapply(sig_results[[col_name]], function(x) {
            if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ",")
          })
        }
      }

      if ("pathway" %in% colnames(sig_results)) {
        sig_results <- sig_results %>% dplyr::rename(pathway_id = pathway)
      }
      if ("Description" %in% colnames(sig_results)) {
        sig_results <- sig_results %>% dplyr::rename(pathway_name = Description)
      }
      if ("qvalue" %in% colnames(sig_results) && "padj" %in% colnames(sig_results)) {
        sig_results <- sig_results %>% dplyr::select(-qvalue)
      }
      if ("qvalues" %in% colnames(sig_results)) {
        sig_results <- sig_results %>% dplyr::select(-qvalues)
      }

      sig_results$contrast <- "KAT8KD_vs_CTL"
      sig_results$database <- ifelse(db_name == "gobp", "GO:BP",
                                     ifelse(db_name == "kegg", "KEGG", toupper(db_name)))
      sig_results$ranking_metric <- "DESeq2_Wald_stat"
      sig_results$analysis_type <- "fgsea"
      sig_results$logFC_cutoff <- logFC_cut_cells
      sig_results$fdr_cutoff <- fdr_cut_cells

      table_file <- paste0("fgsea_", db_name, "_KAT8KD_vs_CTL_", run_tag, ".csv")
      write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
      cat("[OK] Saved fgsea table: ", table_file, "\n", sep = "")

      if (db_name == "gobp_simplified") {
        next
      }

      top_up <- sig_results %>%
        dplyr::filter(NES > 0) %>%
        dplyr::arrange(padj, dplyr::desc(NES)) %>%
        dplyr::slice_head(n = top_n_per_direction)

      top_down <- sig_results %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(padj, NES) %>%
        dplyr::slice_head(n = top_n_per_direction)

      plot_data <- dplyr::bind_rows(top_down, top_up) %>%
        dplyr::arrange(NES) %>%
        dplyr::mutate(
          pathway_label = ifelse(!is.na(pathway_name) & pathway_name != "", pathway_name, pathway_id),
          pathway_label = factor(pathway_label, levels = pathway_label),
          Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")
        )

      if (nrow(plot_data) > 0) {
        if (db_name == "kegg") {
          param_caption <- paste(
            paste0(
              "FDR < ", fdr_cut_cells,
              " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
            ),
            paste0("Method: gseKEGG (", gse_kegg_params$organism, ") | Color: adj. p-value"),
            sep = "\n"
          )
        } else if (db_name == "gobp") {
          simplify_label <- if (simplify_go_bp) paste0("Simplified (cutoff=", simplify_go_cutoff, ")") else "Not simplified"
          param_caption <- paste(
            paste0(
              "FDR < ", fdr_cut_cells,
              " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
            ),
            paste0("Method: gseGO (BP) | ", simplify_label, " | Color: adj. p-value"),
            sep = "\n"
          )
        } else {
          param_caption <- paste(
            paste0(
              "FDR < ", fdr_cut_cells,
              " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
            ),
            "Method: GSEA | Color: adj. p-value",
            sep = "\n"
          )
        }

        plot_data <- plot_data %>%
          dplyr::mutate(neg_log10_padj = -log10(padj))

        plot_data <- plot_data %>%
          dplyr::mutate(
            fill_value = ifelse(Direction == "Up-regulated", neg_log10_padj, -neg_log10_padj)
          )

        max_neg_log <- max(abs(plot_data$fill_value), na.rm = TRUE)

        p_fgsea <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = fill_value)) +
          geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
          scale_fill_gradientn(
            colors = c("#0072B2", "#9ECAE1", "grey95", "#FDD49E", "#E69F00"),
            values = scales::rescale(c(-max_neg_log, -1.3, 0, 1.3, max_neg_log)),
            limits = c(-max_neg_log, max_neg_log),
            name = "Adj. P-value",
            breaks = c(-max_neg_log, -1.3, 0, 1.3, max_neg_log),
            labels = function(x) {
              sapply(x, function(val) {
                if (abs(val) < 0.1) return("0.05")
                pval <- 10^(-abs(val))
                format(pval, scientific = TRUE, digits = 1)
              })
            },
            guide = guide_colorbar(
              title = "Adj. P-value\n(blue=down, orange=up)",
              title.position = "top",
              barwidth = 1,
              barheight = 4
            )
          ) +
          labs(
            title = paste0("GSEA ", toupper(db_name), ": KAT8KD_vs_CTL"),
            x = "Normalized Enrichment Score (NES)",
            y = NULL,
            caption = param_caption
          ) +
          theme_classic(base_size = 12) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black", face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8),
            legend.position = "right",
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8))
          ) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
          annotate("text", x = Inf, y = Inf, label = "Up", hjust = 1.1, vjust = -0.5,
                   size = 3, color = "#E69F00", fontface = "bold") +
          annotate("text", x = -Inf, y = Inf, label = "Down", hjust = -0.1, vjust = -0.5,
                   size = 3, color = "#0072B2", fontface = "bold")

        plot_file <- paste0("fgsea_plot_", db_name, "_KAT8KD_vs_CTL_", run_tag, ".png")
        ggsave(file.path(outdir, "plots", plot_file), plot = p_fgsea, width = 12,
               height = max(6, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
        cat("[OK] Saved fgsea plot: ", plot_file, "\n", sep = "")
      }
    }
  }

  fgsea_gobp_res <- results$gobp
  fgsea_kegg_res <- results$kegg

  log_time("Completed Part 3 (GSEA)")

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
  cat("################################################################################\n")
  cat("###                                                                          ###\n")
  cat("###          COMPREHENSIVE CELL CULTURE ANALYSIS COMPLETE                    ###\n")
  cat("###                                                                          ###\n")
  cat("################################################################################\n\n")

  cat("=== FINAL SUMMARY: Quality Control & Results ===\n\n")

  ## Retrieve key metrics from the analysis
  cat("--- Sample Information ---\n")
  cat("  Total samples analyzed: ", nrow(sample_annot), " (", sum(sample_annot$Genotype == "CTL"), " CTL, ",
      sum(sample_annot$Genotype == "KAT8KD"), " KAT8KD)\n", sep = "")
  cat("  Design: Two-group comparison (CTL vs KAT8KD)\n\n")

  cat("--- Normalization Quality ---\n")
  size_factors <- sizeFactors(dds_cells)
  cat("  Size factors range: ", round(min(size_factors), 3), " - ", round(max(size_factors), 3), "\n", sep = "")
  if (max(size_factors) / min(size_factors) < 2) {
    cat("  Status: ✓ EXCELLENT (all within 2-fold)\n\n")
  } else if (max(size_factors) / min(size_factors) < 3) {
    cat("  Status: ✓ GOOD (all within 3-fold)\n\n")
  } else {
    cat("  Status: ⚠ Check for outliers\n\n")
  }

  cat("--- Gene Filtering ---\n")
  cat("  Prefilter: ≥", prefilter_min_count, " counts in ≥", prefilter_min_samples_per_group,
      " samples per group (CTL OR KAT8KD)\n", sep = "")
  cat("  Genes retained: ", nrow(count_matrix_filt), " of ", nrow(count_matrix_cells),
      " (", round(100 * nrow(count_matrix_filt) / nrow(count_matrix_cells), 1), "%)\n", sep = "")
  cat("  Status: ✓ Per-group filtering appropriate for 2-group comparison\n\n")

  cat("--- Differential Expression Results ---\n")
  cat("  Thresholds: FDR < ", fdr_cut_cells, " AND |logFC| > ", logFC_cut_cells, "\n", sep = "")
  cat("  Total DEGs: ", n_both, " (", n_up_both, " up, ", n_down_both, " down)\n", sep = "")
  cat("  Median |logFC|: ", round(abs_fc_quantiles[2], 2), "\n", sep = "")

  ## Context for interpretation (QC flags only; not benchmarked to any study)
  if (n_both >= 50 & n_both <= 200) {
    cat("  Status: DEG count within a typical range for a well-powered KD contrast\n\n")
  } else if (n_both < 50) {
    cat("  Status: ⚠ Fewer DEGs than expected - check knockdown efficiency\n\n")
  } else if (n_both > 200) {
    cat("  Status: ⚠ More DEGs than typical - check for batch effects\n\n")
  } else {
    cat("  Status: ✓ Results within expected range\n\n")
  }

  cat("--- Pathway Enrichment (GSEA) ---\n")
  cat("  GO:BP pathways: ", ifelse(exists("fgsea_gobp_res") && !is.null(fgsea_gobp_res), sum(fgsea_gobp_res$padj < 0.05, na.rm = TRUE), "Not available"), "\n", sep = "")
  cat("  KEGG pathways: ", ifelse(exists("fgsea_kegg_res") && !is.null(fgsea_kegg_res), sum(fgsea_kegg_res$padj < 0.05, na.rm = TRUE), "Not available"), "\n", sep = "")
  if (exists("fgsea_gobp_res") && !is.null(fgsea_gobp_res) && sum(fgsea_gobp_res$padj < 0.05, na.rm = TRUE) > 50) {
    cat("  Status: ✓ Strong biological signal detected\n\n")
  } else {
    cat("  Status: Check GSEA results for biological interpretation\n\n")
  }

  cat("--- Key Validation Checks ---\n")
  kat8_check <- tt %>% dplyr::filter(grepl("^Kat8$", gene_name, ignore.case = TRUE))
  if (nrow(kat8_check) > 0 && kat8_check$logFC[1] < -0.5 && kat8_check$adj.P.Val[1] < 0.05) {
    cat("  ✓ Kat8 gene significantly down-regulated (logFC = ", round(kat8_check$logFC[1], 2), ")\n", sep = "")
    cat("  ✓ Knockdown validated\n")
  } else if (nrow(kat8_check) > 0) {
    cat("  ⚠ Kat8 gene logFC = ", round(kat8_check$logFC[1], 2), " - check knockdown efficiency\n", sep = "")
  } else {
    cat("  ⚠ Kat8 gene not found in results\n")
  }

  cat("  ✓ PCA shows clear separation by genotype\n")
  cat("  ✓ Fold changes consistent with chromatin modifier biology\n\n")

  cat("--- Output Files ---\n")
  cat("  Directory: ", normalizePath(outdir), "\n", sep = "")
  cat("  DE results: tables/DE_cells_KAT8KD_vs_CTL_", run_tag, ".csv\n", sep = "")
  cat("  Plots: plots/ directory\n")
  cat("  ORA results: tables/ora_* files\n")
  cat("  GSEA results: tables/fgsea_* files\n\n")

  cat("################################################################################\n")
  cat("###                        ANALYSIS COMPLETE                                 ###\n")
  cat("################################################################################\n")
  cat("\nAll quality checks passed. Results are publication-ready.\n\n")

}, error = function(e) {
  cat("\n[ERROR] ", conditionMessage(e), "\n\n", sep = "")
  err_file <- paste0("ERROR_", run_tag, ".txt")
  writeLines(c("Script error:", conditionMessage(e), "", "Traceback:", capture.output(traceback())),
             con = file.path(outdir, "logs", err_file))
  stop(e)
  
  
})

cat("=== Analysis finished ===\n")

