#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 1: MAIN GENE-LEVEL ANALYSIS
## =========================================================
## This script performs:
## - Quality control (density, PCA, MDS)
## - Differential expression (DESeq2)
## - Volcano plots
## - Heatmaps (top DEGs and genes of interest)
## - ComplexHeatmap and EnhancedVolcano for genes of interest
##
## NO PATHWAY ANALYSIS (see part2_ora.R and part3_fgsea.R)
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c(
  "dplyr",
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "viridis",
  "pheatmap",
  "grid",
  "scales"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "DESeq2",
  "ComplexHeatmap",
  "EnhancedVolcano"
)

## circlize is on CRAN
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(viridis)
  library(pheatmap)
  library(grid)   ## for unit()
  library(scales)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  library(EnhancedVolcano)
})

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
source(file.path(utils_dir, "save.R"))
source(file.path(utils_dir, "findsave.R"))
source(file.path(utils_dir, "purge.R"))
source(file.path(utils_dir, "dedupe.R"))

## 2.5) Run metadata ----------------------------------------
run_info <- list(
  script_name = "part1_main_analysis.R",
  species     = "mouse",
  data_type   = "bulkRNAseq_WAT_tissue",
  keywords    = c("KAT8", "WAT", "iWAT", "gWAT", "bulk RNA-seq", "tissue", "DESeq2"),
  notes       = "Part 1: Gene-level analysis only (no pathway analysis)",
  message     = "DESeq2 + VST QC + DEGs + ComplexHeatmap + EnhancedVolcano"
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
## Genes of interest for focused heatmap/volcano
## -----------------------------
genes_of_interest <- c(
  ## Collagens / ECM
  "Col4a1", "Col4a2", "Col15a1", "Col5a3", "Col6a6", "Col13a1",
  "Col4a5", "Col14a1", "Col27a1", "Col16a1", "Col6a5",
  ## Adipocyte markers & metabolism
  "Lep", "Ppargc1a", "Sorbs1", "Srebf1", "Ppara", "Pparg",
  "Fabp4", "Cd36", "Plin1", "Fabp5", "Angptl4", "Cav1",
  ## Matrix remodeling
  "Fn1", "Mmp3", "Timp4", "Mmp12", "Mmp14", "Mmp16"
)
cat("[INFO] Genes of interest for focused plots: ", length(genes_of_interest), "\n")

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

  ## 4.2) Define tissue samples + annotation ----------------
  sample_ids <- paste0("JS_", sprintf("%02d", 1:40))

  sample_annot_full <- data.frame(
    Sample   = sample_ids,
    Depot    = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
    Sex      = c(
      rep("F", 5),  rep("F", 5),
      rep("M", 5),  rep("M", 5),
      rep("F", 5),  rep("F", 5),
      rep("M", 5),  rep("M", 5)
    ),
    Genotype = c(
      rep("CTL", 5),    rep("KAT8KD", 5),
      rep("CTL", 5),    rep("KAT8KD", 5),
      rep("CTL", 5),    rep("KAT8KD", 5),
      rep("CTL", 5),    rep("KAT8KD", 5)
    ),
    stringsAsFactors = FALSE
  )

  ## Outliers (removed before analysis)
  outlier_samples <- c("JS_08", "JS_28")

  ## 4.2.1) QC evidence log BEFORE outlier removal ----------
  cat("\n=== QC BEFORE OUTLIER REMOVAL (log only) ===\n")

  all_tissue_samples <- sample_annot_full$Sample
  missing_all <- setdiff(all_tissue_samples, colnames(counts_raw))
  if (length(missing_all) > 0) stop("Missing columns for QC: ", paste(missing_all, collapse = ", "))

  counts_all_tissue <- as.matrix(counts_raw[, all_tissue_samples])
  rownames(counts_all_tissue) <- counts_raw$gene_id

  ## Library sizes from raw counts (simple + transparent)
  lib_sizes <- colSums(counts_all_tissue, na.rm = TRUE)
  lib_summary <- summary(lib_sizes)
  outlier_lib_sizes <- lib_sizes[names(lib_sizes) %in% outlier_samples]

  outlier_log_file <- file.path(outdir, "logs", paste0("outliers_", run_tag, ".txt"))
  outlier_log_lines <- c(
    "=== Outlier Removal Rationale ===",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Outlier Sample IDs ---",
    paste0("  ", paste(outlier_samples, collapse = ", ")),
    "",
    "--- Library Sizes (All Samples; colSums raw counts) ---",
    paste0("  Min:     ", format(lib_summary["Min."], big.mark = ",")),
    paste0("  1st Qu:  ", format(lib_summary["1st Qu."], big.mark = ",")),
    paste0("  Median:  ", format(lib_summary["Median"], big.mark = ",")),
    paste0("  Mean:    ", format(lib_summary["Mean"], big.mark = ",")),
    paste0("  3rd Qu:  ", format(lib_summary["3rd Qu."], big.mark = ",")),
    paste0("  Max:     ", format(lib_summary["Max."], big.mark = ",")),
    "",
    "--- Outlier Library Sizes ---",
    if (length(outlier_lib_sizes) == 0) "  (none matched?)" else paste0("  ", names(outlier_lib_sizes), ": ", format(outlier_lib_sizes, big.mark = ",")),
    "",
    "--- Rationale ---",
    "  Flagged based on prior QC review and consistency across samples."
  )
  writeLines(outlier_log_lines, con = outlier_log_file)
  cat("Saved outlier rationale log to ", outlier_log_file, "\n")
  cat("=== END QC BEFORE OUTLIER REMOVAL ===\n\n")

  ## Continue with outlier removal
  sample_annot <- sample_annot_full %>% dplyr::filter(!(Sample %in% outlier_samples))
  tissue_samples <- sample_annot$Sample

  missing_tissue <- setdiff(tissue_samples, colnames(counts_raw))
  if (length(missing_tissue) > 0) {
    stop("These tissue sample columns are missing in counts.txt: ", paste(missing_tissue, collapse = ", "))
  }

  counts_tissue_raw <- counts_raw[, c("gene_id", "gene", tissue_samples)]

  cat("\n=== TISSUE SUBSET (outliers removed) ===\n")
  cat("Dimensions (rows x cols):", paste(dim(counts_tissue_raw), collapse = " x "), "\n")

  ## 4.3) Gene cleaning -------------------------------------
  counts_tissue <- counts_tissue_raw %>%
    mutate(
      ensembl_gene_id = gene_id,
      gene_name       = gene
    )

  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]

  ## Handle duplicates robustly
  if (any(duplicated(counts_tissue$gene_name))) {
    dup_flag <- duplicated(counts_tissue$gene_name) | duplicated(counts_tissue$gene_name, fromLast = TRUE)
    counts_tissue$gene_name[dup_flag] <- paste0(counts_tissue$ensembl_gene_id[dup_flag], "_", counts_tissue$gene_name[dup_flag])
  }
  ## As a final safety net:
  counts_tissue$gene_name <- make.unique(counts_tissue$gene_name)

  rownames(counts_tissue) <- counts_tissue$gene_name

  ## 4.4) Build count matrix --------------------------------
  count_matrix_tissue <- as.matrix(counts_tissue[, tissue_samples])
  rownames(count_matrix_tissue) <- counts_tissue$gene_name

  ## 4.5) Attach sample annotation --------------------------
  sample_annot <- sample_annot[match(colnames(count_matrix_tissue), sample_annot$Sample), ]
  if (any(is.na(sample_annot$Sample))) stop("Sample annotation failed to match count matrix columns.")

  sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
  sample_annot$DepotSex <- factor(paste(sample_annot$Depot, sample_annot$Sex, sep = "_"),
                                  levels = c("iWAT_F", "iWAT_M", "gWAT_F", "gWAT_M"))
  sample_annot$GroupFull <- factor(
    paste(sample_annot$Depot, sample_annot$Sex, sample_annot$Genotype, sep = "_"),
    levels = c(
      "iWAT_F_CTL",   "iWAT_F_KAT8KD",
      "iWAT_M_CTL",   "iWAT_M_KAT8KD",
      "gWAT_F_CTL",   "gWAT_F_KAT8KD",
      "gWAT_M_CTL",   "gWAT_M_KAT8KD"
    )
  )
  rownames(sample_annot) <- sample_annot$Sample

  cat("\n=== BASIC QC (TISSUE) ===\n")
  cat("Number of genes:   ", nrow(count_matrix_tissue), "\n")
  cat("Number of samples: ", ncol(count_matrix_tissue), "\n")
  cat("Library sizes (colSums raw counts):\n")
  print(colSums(count_matrix_tissue))

  ## 4.6) DESeq2 object + typical prefilter ------------------
  min_count <- 10
  min_samples <- min(table(sample_annot$GroupFull))

  keep <- rowSums(count_matrix_tissue >= min_count) >= min_samples

  ## Comprehensive filtering sanity check
  n_before <- nrow(count_matrix_tissue)
  n_after <- sum(keep)
  n_removed <- n_before - n_after
  pct_retained <- round(100 * n_after / n_before, 1)

  cat("\n==========================================================\n")
  cat("=== SANITY CHECK: GENE PREFILTERING ===\n")
  cat("==========================================================\n")
  cat("Prefilter rule: counts >= ", min_count, " in >= ", min_samples, " samples\n", sep = "")
  cat("  (min_samples = smallest group size for per-stratum contrasts)\n\n")
  cat("Genes BEFORE prefiltering:  ", n_before, "\n", sep = "")
  cat("Genes AFTER prefiltering:   ", n_after, "\n", sep = "")
  cat("Genes REMOVED:              ", n_removed, " (", 100 - pct_retained, "%)\n", sep = "")
  cat("Genes RETAINED:             ", pct_retained, "%\n", sep = "")
  cat("==========================================================\n\n")

  count_matrix_filt <- count_matrix_tissue[keep, , drop = FALSE]

  dds_tissue <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filt),
    colData   = sample_annot,
    design    = ~ 0 + GroupFull
  )

  ## Save centralized parameters
  logFC_cut_tissue <- 1
  fdr_cut_tissue   <- 0.05

  params <- list(
    outlier_samples   = outlier_samples,
    prefilter_rule    = list(
      min_count       = min_count,
      min_samples     = min_samples,
      description     = "Keep genes with counts >= min_count in at least min_samples (smallest group size)"
    ),
    de_thresholds     = list(
      logFC_cutoff    = logFC_cut_tissue,
      fdr_cutoff      = fdr_cut_tissue
    ),
    vst              = list(
      blind           = TRUE,
      description     = "VST used for QC (blind=TRUE); heatmaps use blind=FALSE"
    ),
    volcano_labels   = list(
      top_n_by_fdr    = 10,
      top_n_by_fc     = 10
    ),
    heatmap_max_genes = 300
  )

  params_file <- file.path(outdir, "logs", paste0("params_", run_tag, ".txt"))
  params_lines <- c(
    "=== RNA-seq Pipeline Parameters (Part 1: Gene-level Analysis) ===",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Outlier Samples ---",
    paste0("  IDs: ", paste(params$outlier_samples, collapse = ", ")),
    "",
    "--- Prefilter Rule (count-based) ---",
    paste0("  min_count: ", params$prefilter_rule$min_count),
    paste0("  min_samples: ", params$prefilter_rule$min_samples),
    paste0("  Description: ", params$prefilter_rule$description),
    "",
    "--- DE Thresholds ---",
    paste0("  |log2FC| cutoff: ", params$de_thresholds$logFC_cutoff),
    paste0("  FDR cutoff: ", params$de_thresholds$fdr_cutoff),
    "",
    "--- VST ---",
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
  dds_tissue <- DESeq(dds_tissue)
  cat("[OK] DESeq2 complete.\n")

  cat("\nDESeq2 design matrix (first rows):\n")
  print(head(model.matrix(design(dds_tissue), colData(dds_tissue))))

  ## 4.8) VST transform for QC ------------------------------
  cat("\n=== VST FOR QC ===\n")
  vst_tissue <- vst(dds_tissue, blind = TRUE)
  vst_mat <- assay(vst_tissue)  ## genes x samples (log2-ish stabilized)

  cat("\n=== VST FOR HEATMAPS (blind=FALSE) ===\n")
  vst_tissue_hm <- vst(dds_tissue, blind = FALSE)
  vst_mat_hm <- assay(vst_tissue_hm)

  genotype_colors <- c("CTL" = "#0072B2", "KAT8KD" = "#E69F00")
  depot_sex_fill <- c(
    "iWAT_F" = "#56B4E9",
    "iWAT_M" = "#0072B2",
    "gWAT_F" = "#F0E442",
    "gWAT_M" = "#E69F00"
  )

  ## 4.9) Density plots (VST before/after filtering) --------
  cat("\n=== Creating before/after filtering density plots ===\n")

  ## Get VST for all genes (before filtering) for comparison
  dds_all_genes <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_tissue),
    colData   = sample_annot,
    design    = ~ 0 + GroupFull
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

    ## Define consistent colors for each sample across both plots
    depot_sex <- colData(dds_tissue)$DepotSex
    sample_colors <- depot_sex_fill[depot_sex]

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

  density_file <- paste0("VST_density_tissue_before_after_", run_tag, ".png")
  png(file.path(outdir, "plots", density_file), width = 1600, height = 900, res = 150)
  plot_density_vst_comparison(vst_mat_before, vst_mat)
  dev.off()
  cat("Saved VST density plot (before/after filtering) to ", file.path(outdir, "plots", density_file), "\n")

  ## 4.10) PCA (VST) ----------------------------------------
  pca_tissue <- prcomp(t(vst_mat), scale. = TRUE)
  pca_var <- summary(pca_tissue)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    Sample   = colnames(vst_mat),
    Genotype = colData(dds_tissue)$Genotype,
    DepotSex = colData(dds_tissue)$DepotSex,
    Depot    = colData(dds_tissue)$Depot,
    Sex      = colData(dds_tissue)$Sex,
    PC1      = pca_tissue$x[, 1],
    PC2      = pca_tissue$x[, 2],
    stringsAsFactors = FALSE
  )

  p_pca <- ggplot(
    pca_df,
    aes(x = PC1, y = PC2, color = DepotSex, fill = DepotSex, shape = Genotype, label = Sample)
  ) +
    geom_point(size = 3, stroke = 1) +
    geom_text_repel(size = 3, max.overlaps = 60, color = "black") +
    scale_color_manual(values = depot_sex_fill, name = "Depot/Sex") +
    scale_fill_manual(values = depot_sex_fill, name = "Depot/Sex") +
    scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
    labs(
      title = "Tissue PCA (VST, DESeq2)",
      x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_var[2], 1), "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right") +
    guides(
      color = "none",
      fill = guide_legend(override.aes = list(shape = 21, color = depot_sex_fill, fill = depot_sex_fill)),
      shape = guide_legend(override.aes = list(fill = "white", color = "black"))
    )

  pca_file <- paste0("PCA_tissue_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", pca_file), plot = p_pca, width = 8, height = 6, dpi = 300)
  cat("Saved PCA plot to ", file.path(outdir, "plots", pca_file), "\n")

  ## 4.11) MDS (VST) ----------------------------------------
  dist_mat <- dist(t(vst_mat))
  mds_coords <- cmdscale(dist_mat, k = 2)

  mds_df <- data.frame(
    Sample   = colnames(vst_mat),
    Genotype = colData(dds_tissue)$Genotype,
    DepotSex = colData(dds_tissue)$DepotSex,
    Depot    = colData(dds_tissue)$Depot,
    Sex      = colData(dds_tissue)$Sex,
    MDS1     = mds_coords[, 1],
    MDS2     = mds_coords[, 2],
    stringsAsFactors = FALSE
  )

  p_mds <- ggplot(
    mds_df,
    aes(x = MDS1, y = MDS2, color = DepotSex, fill = DepotSex, shape = Genotype, label = Sample)
  ) +
    geom_point(size = 3, stroke = 1) +
    geom_text_repel(size = 3, max.overlaps = 60, color = "black") +
    scale_color_manual(values = depot_sex_fill, name = "Depot/Sex") +
    scale_fill_manual(values = depot_sex_fill, name = "Depot/Sex") +
    scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
    labs(title = "Tissue MDS (VST, DESeq2)", x = "MDS1", y = "MDS2") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right") +
    guides(
      color = "none",
      fill = guide_legend(override.aes = list(shape = 21, color = depot_sex_fill, fill = depot_sex_fill)),
      shape = guide_legend(override.aes = list(fill = "white", color = "black"))
    )

  mds_file <- paste0("MDS_tissue_VST_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", mds_file), plot = p_mds, width = 8, height = 6, dpi = 300)
  cat("Saved MDS plot to ", file.path(outdir, "plots", mds_file), "\n")

  ## 4.12) Contrasts ----------------------------------------
  contrast_definitions <- list(
    iWAT_F_KD_vs_CTL = c("GroupFull", "iWAT_F_KAT8KD", "iWAT_F_CTL"),
    iWAT_M_KD_vs_CTL = c("GroupFull", "iWAT_M_KAT8KD", "iWAT_M_CTL"),
    gWAT_F_KD_vs_CTL = c("GroupFull", "gWAT_F_KAT8KD", "gWAT_F_CTL"),
    gWAT_M_KD_vs_CTL = c("GroupFull", "gWAT_M_KAT8KD", "gWAT_M_CTL")
  )
  contrast_names <- names(contrast_definitions)

  ## 4.13) Containers ---------------------------------------
  deg_summary <- data.frame(
    Contrast          = contrast_names,
    N_DEG_FDR_lt_0_05 = NA_integer_,
    N_Up              = NA_integer_,
    N_Down            = NA_integer_,
    stringsAsFactors  = FALSE
  )

  tt_list <- vector("list", length(contrast_names))
  names(tt_list) <- contrast_names

  ## 4.14) Per-contrast DE + volcano + heatmap + genes of interest ----
  for (cn in contrast_names) {

    cat("\n=== Contrast: ", cn, " ===\n", sep = "")

    res <- results(dds_tissue, contrast = contrast_definitions[[cn]])
    res <- res[order(res$pvalue), , drop = FALSE]

    tt <- as.data.frame(res)
    tt$logFC     <- tt$log2FoldChange
    tt$P.Value   <- tt$pvalue
    tt$adj.P.Val <- tt$padj
    tt$gene_name <- rownames(tt)

    ## Remove duplicate columns
    tt <- tt[, !duplicated(colnames(tt)), drop = FALSE]

    tt_list[[cn]] <- tt

    de_file <- paste0("DE_tissue_", cn, "_", run_tag, ".csv")
    write.csv(tt, file = file.path(outdir, "tables", de_file), row.names = TRUE)

    ## Calculate comprehensive DEG statistics
    n_total <- nrow(tt)
    n_valid_pval <- sum(!is.na(tt$adj.P.Val))
    n_fdr_only <- sum(tt$adj.P.Val < fdr_cut_tissue, na.rm = TRUE)
    n_up_both <- sum(tt$adj.P.Val < fdr_cut_tissue & tt$logFC > logFC_cut_tissue, na.rm = TRUE)
    n_down_both <- sum(tt$adj.P.Val < fdr_cut_tissue & tt$logFC < -logFC_cut_tissue, na.rm = TRUE)
    n_both <- n_up_both + n_down_both

    cat("\n==========================================================\n")
    cat("=== SANITY CHECK: ", cn, " ===\n", sep = "")
    cat("==========================================================\n")
    cat("Total genes tested by DESeq2:           ", n_total, "\n", sep = "")
    cat("Genes with valid adjusted p-value:      ", n_valid_pval, "\n", sep = "")
    cat("\n--- DEG Counts at Different Thresholds ---\n")
    cat("DEGs (FDR < ", fdr_cut_tissue, " only):               ", n_fdr_only, "\n", sep = "")
    cat("DEGs (FDR < ", fdr_cut_tissue, " AND |logFC| > ", logFC_cut_tissue, "):  ", n_both, "\n", sep = "")
    cat("  - Up-regulated (logFC > ", logFC_cut_tissue, "):       ", n_up_both, "\n", sep = "")
    cat("  - Down-regulated (logFC < -", logFC_cut_tissue, "):    ", n_down_both, "\n", sep = "")
    cat("\nPercent of tested genes that are DEGs:\n")
    cat("  - FDR-only threshold:                 ", round(100 * n_fdr_only / n_total, 2), "%\n", sep = "")
    cat("  - FDR + logFC threshold:              ", round(100 * n_both / n_total, 2), "%\n", sep = "")
    cat("==========================================================\n")
    cat("DE table written: ", file.path(outdir, "tables", de_file), "\n", sep = "")
    cat("==========================================================\n\n")

    deg_summary$N_DEG_FDR_lt_0_05[deg_summary$Contrast == cn] <- n_fdr_only
    deg_summary$N_Up[deg_summary$Contrast == cn] <- n_up_both
    deg_summary$N_Down[deg_summary$Contrast == cn] <- n_down_both

    ## Volcano
    tt_plot <- tt %>%
      dplyr::filter(is.finite(logFC), is.finite(adj.P.Val), adj.P.Val > 0) %>%
      mutate(negLogFDR = -log10(adj.P.Val))

    tt_plot <- tt_plot %>%
      mutate(
        sig_cat = dplyr::case_when(
          adj.P.Val < fdr_cut_tissue & logFC >=  logFC_cut_tissue ~ "Up",
          adj.P.Val < fdr_cut_tissue & logFC <= -logFC_cut_tissue ~ "Down",
          TRUE                                                    ~ "NS"
        ),
        sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
      )

    tt_sig <- tt_plot %>%
      dplyr::filter(adj.P.Val < fdr_cut_tissue,
                    abs(logFC) >= logFC_cut_tissue,
                    sig_cat != "NS")

    top_n_fdr <- params$volcano_labels$top_n_by_fdr
    top_n_fc  <- params$volcano_labels$top_n_by_fc

    up_sig <- tt_sig %>% dplyr::filter(sig_cat == "Up")
    down_sig <- tt_sig %>% dplyr::filter(sig_cat == "Down")

    up_top_fdr <- up_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = top_n_fdr)
    up_top_fc  <- up_sig %>% dplyr::arrange(dplyr::desc(logFC)) %>% dplyr::slice_head(n = top_n_fc)

    down_top_fdr <- down_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = top_n_fdr)
    down_top_fc  <- down_sig %>% dplyr::arrange(logFC) %>% dplyr::slice_head(n = top_n_fc)

    label_genes <- dplyr::bind_rows(up_top_fdr, up_top_fc, down_top_fdr, down_top_fc) %>%
      dplyr::distinct(gene_name, .keep_all = TRUE)

    vol_tissue <- ggplot(
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
      geom_hline(yintercept = -log10(fdr_cut_tissue), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(logFC_cut_tissue, -logFC_cut_tissue), linetype = "dashed", color = "grey40") +
      theme_bw(base_size = 14) +
      theme(
        panel.grid      = element_blank(),
        axis.text       = element_text(size = 10, face = "bold", colour = "black"),
        axis.title      = element_text(size = 14, face = "bold", colour = "black"),
        plot.title      = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
      ) +
      ggtitle(paste0("Volcano: ", cn, " (WAT tissue; DESeq2)")) +
      xlab(expression(log[2]~Fold~Change)) +
      ylab(expression(-log[10]~FDR))

    volcano_file <- paste0("Volcano_tissue_", cn, "_", run_tag, ".png")
    ggsave(file.path(outdir, "plots", volcano_file), plot = vol_tissue, width = 8, height = 6, dpi = 300)
    cat("Volcano saved: ", file.path(outdir, "plots", volcano_file), "\n", sep = "")

    if (is.null(hero_volcano_file)) hero_volcano_file <- volcano_file

    ## Heatmap (VST) of strongest DEGs using ComplexHeatmap
    de_sig <- tt %>%
      dplyr::filter(adj.P.Val < fdr_cut_tissue, abs(logFC) > logFC_cut_tissue)

    if (nrow(de_sig) >= 2) {
      de_sig <- de_sig[order(de_sig$adj.P.Val), , drop = FALSE]
      if (nrow(de_sig) > params$heatmap_max_genes) de_sig <- de_sig[1:params$heatmap_max_genes, , drop = FALSE]
      gene_set <- de_sig$gene_name

      depot <- NA_character_
      sex   <- NA_character_
      if (grepl("^iWAT_F", cn))      { depot <- "iWAT"; sex <- "F" }
      else if (grepl("^iWAT_M", cn)) { depot <- "iWAT"; sex <- "M" }
      else if (grepl("^gWAT_F", cn)) { depot <- "gWAT"; sex <- "F" }
      else if (grepl("^gWAT_M", cn)) { depot <- "gWAT"; sex <- "M" }

      if (!is.na(depot) && !is.na(sex)) {

        idx_samples <- with(as.data.frame(colData(dds_tissue)), Depot == depot & Sex == sex)
        samples_heat <- rownames(colData(dds_tissue))[idx_samples]
        heat_meta_deg <- as.data.frame(colData(dds_tissue))[idx_samples, , drop = FALSE]

        mat <- vst_mat_hm[gene_set, samples_heat, drop = FALSE]

        ## Scale rows (z-score)
        mat_scaled <- t(scale(t(mat)))
        mat_scaled <- mat_scaled[!rowSums(is.na(mat_scaled)), , drop = FALSE]

        if (nrow(mat_scaled) >= 2) {
          ## Use volcano colors (teal to white to orange)
          max_val_deg <- max(abs(mat_scaled), na.rm = TRUE)
          deg_col_fun <- colorRamp2(
            c(-max_val_deg, 0, max_val_deg),
            c("#0072B2", "white", "#E69F00")  ## Teal (down) -> white -> orange (up)
          )

          ## Column annotation
          deg_ha <- HeatmapAnnotation(
            Genotype = heat_meta_deg$Genotype,
            col = list(Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")),
            annotation_name_side = "left"
          )

          ## Create heatmap
          heat_deg <- Heatmap(
            mat_scaled,
            col = deg_col_fun,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_column_names = TRUE,
            show_row_names = (nrow(mat_scaled) <= 50),
            column_split = heat_meta_deg$Genotype,
            border = TRUE,
            column_gap = unit(2, "mm"),
            row_gap = unit(0, "mm"),
            top_annotation = deg_ha,
            column_title = paste0("Top DEGs: ", cn, " (FDR<0.05, |logFC|>", logFC_cut_tissue, ")"),
            heatmap_legend_param = list(
              title = "Z-score",
              title_position = "leftcenter-rot",
              title_gp = gpar(fontsize = 10)
            ),
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 8)
          )

          ## Save
          heat_file <- paste0("Heatmap_tissue_", cn, "_", run_tag, ".png")
          png(file.path(outdir, "plots", heat_file), width = 2000, height = 3000, res = 200)
          draw(heat_deg)
          dev.off()
          cat("Heatmap saved: ", file.path(outdir, "plots", heat_file), "\n", sep = "")
        } else {
          cat("Not enough valid genes for heatmap after scaling (", cn, ").\n", sep = "")
        }

      } else {
        cat("Could not map contrast ", cn, " to depot/sex for heatmap.\n", sep = "")
      }
    } else {
      cat("Not enough DE genes for heatmap (", cn, ").\n", sep = "")
    }

    ## -------------------------------------------------------
    ## ComplexHeatmap for genes of interest
    ## -------------------------------------------------------
    cat("\n--- ComplexHeatmap: Genes of Interest (", cn, ") ---\n", sep = "")

    ## Find which genes of interest are in our data
    goi_in_data <- genes_of_interest[genes_of_interest %in% rownames(vst_mat)]
    cat("[INFO] Genes of interest found in data: ", length(goi_in_data), "/", length(genes_of_interest), "\n", sep = "")

    if (length(goi_in_data) >= 2) {

      ## Get depot/sex for this contrast
      depot_h <- NA_character_
      sex_h   <- NA_character_
      if (grepl("^iWAT_F", cn))      { depot_h <- "iWAT"; sex_h <- "F" }
      else if (grepl("^iWAT_M", cn)) { depot_h <- "iWAT"; sex_h <- "M" }
      else if (grepl("^gWAT_F", cn)) { depot_h <- "gWAT"; sex_h <- "F" }
      else if (grepl("^gWAT_M", cn)) { depot_h <- "gWAT"; sex_h <- "M" }

      if (!is.na(depot_h) && !is.na(sex_h)) {

        ## Get sample indices for this contrast
        heat_idx <- with(as.data.frame(colData(dds_tissue)), Depot == depot_h & Sex == sex_h)
        heat_samples <- rownames(colData(dds_tissue))[heat_idx]
        heat_meta <- as.data.frame(colData(dds_tissue))[heat_idx, , drop = FALSE]

        ## Get control samples for centering
        heat_ctrls <- heat_meta %>% dplyr::filter(Genotype == "CTL") %>% rownames()

        ## Build matrix
        heat_matrix <- vst_mat_hm[goi_in_data, heat_samples, drop = FALSE]

        ## Center by control means (z-score relative to controls)
        heat_ctrl_means <- rowMeans(heat_matrix[, heat_ctrls, drop = FALSE])
        heat_matrix_scaled <- t(scale(t(heat_matrix), center = heat_ctrl_means, scale = TRUE))
        heat_matrix_scaled <- heat_matrix_scaled[!rowSums(is.na(heat_matrix_scaled)), , drop = FALSE]

        if (nrow(heat_matrix_scaled) >= 2) {

          ## Use volcano colors (teal to white to orange)
          max_val <- max(abs(heat_matrix_scaled), na.rm = TRUE)
          heat_col_fun <- colorRamp2(
            c(-max_val, 0, max_val),
            c("#0072B2", "white", "#E69F00")  ## Teal (down) -> white -> orange (up)
          )

          ## Column annotation
          heat_ha <- HeatmapAnnotation(
            Genotype = heat_meta$Genotype,
            col = list(Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")),
            annotation_name_side = "left"
          )

          ## Create heatmap
          heat_goi <- Heatmap(
            heat_matrix_scaled,
            col = heat_col_fun,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            show_column_names = TRUE,
            show_row_names = TRUE,
            column_split = heat_meta$Genotype,
            border = TRUE,
            column_gap = unit(2, "mm"),
            row_gap = unit(0, "mm"),
            top_annotation = heat_ha,
            column_title = paste0("Genes of Interest: ", cn),
            heatmap_legend_param = list(
              title = "Z-score\n(vs CTL)",
              title_position = "leftcenter-rot",
              title_gp = gpar(fontsize = 10)
            ),
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8)
          )

          ## Save
          heat_goi_file <- paste0("Heatmap_GOI_", cn, "_", run_tag, ".png")
          png(file.path(outdir, "plots", heat_goi_file), width = 2000, height = 2500, res = 200)
          draw(heat_goi)
          dev.off()
          cat("[OK] ComplexHeatmap saved: ", heat_goi_file, "\n", sep = "")

        } else {
          cat("[WARN] Not enough genes of interest with valid values for heatmap\n")
        }
      }
    } else {
      cat("[WARN] Too few genes of interest found in data for heatmap\n")
    }

    ## -------------------------------------------------------
    ## EnhancedVolcano labeling genes of interest
    ## -------------------------------------------------------
    cat("\n--- EnhancedVolcano: Genes of Interest (", cn, ") ---\n", sep = "")

    ## Prepare data for EnhancedVolcano
    volcano_df <- tt %>%
      dplyr::filter(is.finite(logFC), is.finite(adj.P.Val), adj.P.Val > 0)

    ## Only label genes of interest that are SIGNIFICANT (pass both thresholds)
    goi_significant <- volcano_df %>%
      dplyr::filter(gene_name %in% genes_of_interest,
                    adj.P.Val < fdr_cut_tissue,
                    abs(logFC) >= logFC_cut_tissue) %>%
      dplyr::pull(gene_name)

    cat("[INFO] Genes of interest (total): ", length(genes_of_interest), "\n", sep = "")
    cat("[INFO] Genes of interest (significant): ", length(goi_significant), "\n", sep = "")

    if (nrow(volcano_df) > 0) {

      ## Create custom color scheme using keyvals
      keyvals_color <- rep("grey70", nrow(volcano_df))
      names(keyvals_color) <- rep("NS", nrow(volcano_df))

      ## Significant UP (orange)
      sig_up_idx <- which(volcano_df$adj.P.Val < fdr_cut_tissue &
                            volcano_df$logFC >= logFC_cut_tissue)
      keyvals_color[sig_up_idx] <- "#E69F00"
      names(keyvals_color)[sig_up_idx] <- "Significant Up"

      ## Significant DOWN (teal/blue)
      sig_down_idx <- which(volcano_df$adj.P.Val < fdr_cut_tissue &
                              volcano_df$logFC <= -logFC_cut_tissue)
      keyvals_color[sig_down_idx] <- "#0072B2"
      names(keyvals_color)[sig_down_idx] <- "Significant Down"

      ev_plot <- EnhancedVolcano(
        volcano_df,
        lab = volcano_df$gene_name,
        selectLab = goi_significant,
        x = "logFC",
        y = "adj.P.Val",
        title = paste0("Volcano: ", cn),
        subtitle = "ECM, metabolism, and remodeling genes labeled (significant only)",
        pCutoff = fdr_cut_tissue,
        FCcutoff = logFC_cut_tissue,
        pointSize = 2.0,
        labSize = 4.0,
        labCol = "black",
        labFace = "bold",
        boxedLabels = TRUE,
        drawConnectors = TRUE,
        widthConnectors = 0.5,
        colConnectors = "gray40",
        colCustom = keyvals_color,
        colAlpha = 0.8,
        legendPosition = "right",
        legendLabSize = 10,
        legendIconSize = 3.0,
        max.overlaps = Inf,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      )

      ev_file <- paste0("EnhancedVolcano_GOI_", cn, "_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", ev_file), plot = ev_plot, width = 10, height = 8, dpi = 300)
      cat("[OK] EnhancedVolcano saved: ", ev_file, "\n", sep = "")

    } else {
      cat("[WARN] No valid data for EnhancedVolcano\n")
    }
  }


  ## 4.15) OVERALL CTL vs KD contrast (all depots/sexes combined) ----
  cat("\n\n=================================================================\n")
  cat("=== OVERALL CONTRAST: CTL vs KD (all depots/sexes combined) ===\n")
  cat("=================================================================\n\n")

  ## Rebuild DESeq2 object with proper design
  dds_overall <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filt),
    colData   = sample_annot,
    design    = ~ Depot + Sex + Genotype
  )

  cat("[INFO] Running DESeq2 for overall contrast...\n")
  dds_overall <- DESeq(dds_overall)
  cat("[OK] DESeq2 complete for overall contrast.\n")

  ## Extract results: KAT8KD vs CTL
  res_overall <- results(dds_overall, contrast = c("Genotype", "KAT8KD", "CTL"))
  res_overall <- res_overall[order(res_overall$pvalue), , drop = FALSE]

  tt_overall <- as.data.frame(res_overall)
  tt_overall$logFC     <- tt_overall$log2FoldChange
  tt_overall$P.Value   <- tt_overall$pvalue
  tt_overall$adj.P.Val <- tt_overall$padj
  tt_overall$gene_name <- rownames(tt_overall)

  ## Remove duplicate columns
  tt_overall <- tt_overall[, !duplicated(colnames(tt_overall)), drop = FALSE]

  ## Save DE table
  de_overall_file <- paste0("DE_tissue_OVERALL_KD_vs_CTL_", run_tag, ".csv")
  write.csv(tt_overall, file = file.path(outdir, "tables", de_overall_file), row.names = TRUE)

  ## Calculate comprehensive DEG statistics for overall contrast
  n_total_overall <- nrow(tt_overall)
  n_valid_pval_overall <- sum(!is.na(tt_overall$adj.P.Val))
  n_sig_overall <- sum(tt_overall$adj.P.Val < fdr_cut_tissue, na.rm = TRUE)
  n_up_overall <- sum(tt_overall$adj.P.Val < fdr_cut_tissue & tt_overall$logFC > logFC_cut_tissue, na.rm = TRUE)
  n_down_overall <- sum(tt_overall$adj.P.Val < fdr_cut_tissue & tt_overall$logFC < -logFC_cut_tissue, na.rm = TRUE)
  n_both_overall <- n_up_overall + n_down_overall

  cat("\n==========================================================\n")
  cat("=== SANITY CHECK: OVERALL_KD_vs_CTL ===\n")
  cat("==========================================================\n")
  cat("Total genes tested by DESeq2:           ", n_total_overall, "\n", sep = "")
  cat("Genes with valid adjusted p-value:      ", n_valid_pval_overall, "\n", sep = "")
  cat("\n--- DEG Counts at Different Thresholds ---\n")
  cat("DEGs (FDR < ", fdr_cut_tissue, " only):               ", n_sig_overall, "\n", sep = "")
  cat("DEGs (FDR < ", fdr_cut_tissue, " AND |logFC| > ", logFC_cut_tissue, "):  ", n_both_overall, "\n", sep = "")
  cat("  - Up-regulated (logFC > ", logFC_cut_tissue, "):       ", n_up_overall, "\n", sep = "")
  cat("  - Down-regulated (logFC < -", logFC_cut_tissue, "):    ", n_down_overall, "\n", sep = "")
  cat("\nPercent of tested genes that are DEGs:\n")
  cat("  - FDR-only threshold:                 ", round(100 * n_sig_overall / n_total_overall, 2), "%\n", sep = "")
  cat("  - FDR + logFC threshold:              ", round(100 * n_both_overall / n_total_overall, 2), "%\n", sep = "")
  cat("==========================================================\n")
  cat("DE table written: ", file.path(outdir, "tables", de_overall_file), "\n", sep = "")
  cat("==========================================================\n\n")

  ## Add to summary table
  deg_summary <- rbind(
    deg_summary,
    data.frame(
      Contrast          = "OVERALL_KD_vs_CTL",
      N_DEG_FDR_lt_0_05 = n_sig_overall,
      N_Up              = n_up_overall,
      N_Down            = n_down_overall,
      stringsAsFactors  = FALSE
    )
  )

  ## Volcano plot for overall contrast
  tt_overall_plot <- tt_overall %>%
    dplyr::filter(is.finite(logFC), is.finite(adj.P.Val), adj.P.Val > 0) %>%
    mutate(negLogFDR = -log10(adj.P.Val))

  tt_overall_plot <- tt_overall_plot %>%
    mutate(
      sig_cat = dplyr::case_when(
        adj.P.Val < fdr_cut_tissue & logFC >=  logFC_cut_tissue ~ "Up",
        adj.P.Val < fdr_cut_tissue & logFC <= -logFC_cut_tissue ~ "Down",
        TRUE                                                    ~ "NS"
      ),
      sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
    )

  ## Label top genes
  tt_overall_sig <- tt_overall_plot %>%
    dplyr::filter(adj.P.Val < fdr_cut_tissue, abs(logFC) >= logFC_cut_tissue, sig_cat != "NS")

  overall_up_sig <- tt_overall_sig %>% dplyr::filter(sig_cat == "Up")
  overall_down_sig <- tt_overall_sig %>% dplyr::filter(sig_cat == "Down")

  overall_up_top_fdr <- overall_up_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = 10)
  overall_up_top_fc  <- overall_up_sig %>% dplyr::arrange(dplyr::desc(logFC)) %>% dplyr::slice_head(n = 10)

  overall_down_top_fdr <- overall_down_sig %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice_head(n = 10)
  overall_down_top_fc  <- overall_down_sig %>% dplyr::arrange(logFC) %>% dplyr::slice_head(n = 10)

  overall_label_genes <- dplyr::bind_rows(
    overall_up_top_fdr, overall_up_top_fc, overall_down_top_fdr, overall_down_top_fc
  ) %>% dplyr::distinct(gene_name, .keep_all = TRUE)

  vol_overall <- ggplot(
    tt_overall_plot,
    aes(x = logFC, y = negLogFDR, color = sig_cat)
  ) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c("Up" = "#E69F00", "Down" = "#0072B2", "NS" = "grey70"),
      name   = "Significance"
    ) +
    geom_point(data = overall_label_genes, aes(color = sig_cat), size = 3) +
    geom_text_repel(
      data          = overall_label_genes,
      aes(label     = gene_name),
      color         = "black",
      size          = 3.5,
      box.padding   = unit(0.25, "lines"),
      segment.color = "gray40",
      segment.size  = 0.3,
      point.padding = unit(0.2, "lines"),
      max.overlaps  = Inf
    ) +
    geom_hline(yintercept = -log10(fdr_cut_tissue), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(logFC_cut_tissue, -logFC_cut_tissue), linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid      = element_blank(),
      axis.text       = element_text(size = 10, face = "bold", colour = "black"),
      axis.title      = element_text(size = 14, face = "bold", colour = "black"),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    ) +
    ggtitle("Volcano: OVERALL KD vs CTL (all depots/sexes combined)") +
    xlab(expression(log[2]~Fold~Change)) +
    ylab(expression(-log[10]~FDR))

  volcano_overall_file <- paste0("Volcano_tissue_OVERALL_KD_vs_CTL_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", volcano_overall_file), plot = vol_overall, width = 8, height = 6, dpi = 300)
  cat("Volcano saved: ", file.path(outdir, "plots", volcano_overall_file), "\n", sep = "")

  ## Heatmap for overall contrast
  cat("\n--- Heatmap: OVERALL contrast ---\n")
  de_sig_overall <- tt_overall %>%
    dplyr::filter(adj.P.Val < fdr_cut_tissue, abs(logFC) > logFC_cut_tissue)

  if (nrow(de_sig_overall) >= 2) {
    de_sig_overall <- de_sig_overall[order(de_sig_overall$adj.P.Val), , drop = FALSE]
    if (nrow(de_sig_overall) > params$heatmap_max_genes) {
      de_sig_overall <- de_sig_overall[1:params$heatmap_max_genes, , drop = FALSE]
    }
    gene_set_overall <- de_sig_overall$gene_name

    ## Get VST matrix for overall (using blind=FALSE VST from the overall dds object)
    vst_overall_hm <- vst(dds_overall, blind = FALSE)
    vst_mat_overall_hm <- assay(vst_overall_hm)

    mat_overall <- vst_mat_overall_hm[gene_set_overall, , drop = FALSE]

    ## Scale rows (z-score)
    mat_overall_scaled <- t(scale(t(mat_overall)))
    mat_overall_scaled <- mat_overall_scaled[!rowSums(is.na(mat_overall_scaled)), , drop = FALSE]

    if (nrow(mat_overall_scaled) >= 2) {
      ## Color function
      max_val_overall <- max(abs(mat_overall_scaled), na.rm = TRUE)
      overall_col_fun <- colorRamp2(
        c(-max_val_overall, 0, max_val_overall),
        c("#0072B2", "white", "#E69F00")
      )

      ## Column annotation (Depot, Sex, Genotype)
      heat_meta_overall <- as.data.frame(colData(dds_overall))
      overall_ha <- HeatmapAnnotation(
        Depot = heat_meta_overall$Depot,
        Sex = heat_meta_overall$Sex,
        Genotype = heat_meta_overall$Genotype,
        col = list(
          Depot = c(iWAT = "#1B9E77", gWAT = "#D95F02"),
          Sex = c(F = "#E7298A", M = "#7570B3"),
          Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")
        ),
        annotation_name_side = "left"
      )

      ## Create heatmap
      heat_overall <- Heatmap(
        mat_overall_scaled,
        col = overall_col_fun,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        show_row_names = (nrow(mat_overall_scaled) <= 50),
        column_split = heat_meta_overall$Genotype,
        border = TRUE,
        column_gap = unit(2, "mm"),
        row_gap = unit(0, "mm"),
        top_annotation = overall_ha,
        column_title = paste0("Top DEGs: OVERALL KD vs CTL (FDR<0.05, |logFC|>", logFC_cut_tissue, ")"),
        heatmap_legend_param = list(
          title = "Z-score",
          title_position = "leftcenter-rot",
          title_gp = gpar(fontsize = 10)
        ),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 6)
      )

      ## Save (increased height to prevent sample name cutoff)
      heat_overall_file <- paste0("Heatmap_tissue_OVERALL_KD_vs_CTL_", run_tag, ".png")
      png(file.path(outdir, "plots", heat_overall_file), width = 2400, height = 4200, res = 200)
      draw(heat_overall, padding = unit(c(2, 2, 40, 2), "mm"))  ## Extra bottom padding for sample names
      dev.off()
      cat("Heatmap saved: ", file.path(outdir, "plots", heat_overall_file), "\n", sep = "")
    } else {
      cat("Not enough valid genes for overall heatmap after scaling.\n")
    }
  } else {
    cat("Not enough DE genes for overall heatmap.\n")
  }

  ## 4.16) DEG summary table --------------------------------
  deg_summary_file <- paste0("DEG_summary_tissue_", run_tag, ".csv")
  write.csv(deg_summary, file = file.path(outdir, "tables", deg_summary_file), row.names = FALSE)
  cat("\nDEG summary written: ", file.path(outdir, "tables", deg_summary_file), "\n", sep = "")

  ## 4.17) save_core bookkeeping ----------------------------
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
  cat("\n======================================\n")
  cat("=== PART 1 COMPLETE ===\n")
  cat("======================================\n")
  cat("DE results saved. Run part2_ora.R and part3_fgsea.R for pathway analysis.\n\n")

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
