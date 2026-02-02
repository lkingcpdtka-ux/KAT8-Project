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

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  source(file.path(utils_dir, "purge.R"))
  source(file.path(utils_dir, "dedupe.R"))
  use_save_core <- TRUE
  cat("[OK] save_core utilities loaded\n")
} else {
  use_save_core <- FALSE
  cat("[WARN] save_core not found at: ", utils_dir, "\n", sep = "")
  cat("[WARN] Creating minimal output directory structure manually...\n")

  ## Create minimal init_run function if save_core is missing
  init_run <- function(script_name, species, data_type, keywords, notes) {
    run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
    outdir <- file.path(getwd(), "savepoints", paste0("RUN_", run_tag))
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
    return(list(outdir = outdir, run_tag = run_tag))
  }
}

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

## -----------------------------
## Prefilter parameters (tune as needed)
## -----------------------------
## Updated to global filtering (2026-01-21)
## Rationale: KAT8 is a chromatin modifier causing context-dependent chromatin remodeling
##   - Global filtering captures genes with condition-specific accessibility changes
##   - Per-group filtering removes genes showing KAT8-dependent ectopic expression
##   - More appropriate for studying genome-wide chromatin effects
prefilter_min_count <- 10
prefilter_min_samples <- 4  # total across experiment, not per group

## 4) Main logic --------------------------------------------
hero_volcano_file <- NULL

tryCatch({
  
  ## 4.1) Load raw counts -----------------------------------
  cat("\n=== RAW COUNTS MATRIX (FULL) ===\n")

  ## Check if counts file exists
  counts_file <- "counts.txt"
  if (!file.exists(counts_file)) {
    stop("counts.txt not found in working directory: ", getwd(), "\n",
         "Please ensure counts.txt is in the project root directory.")
  }

  counts_raw <- read.delim(
    counts_file,
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
  sample_annot$Depot <- factor(sample_annot$Depot, levels = c("iWAT", "gWAT"))
  sample_annot$Sex <- factor(sample_annot$Sex, levels = c("F", "M"))
  sample_annot$DepotSex <- factor(paste(sample_annot$Depot, sample_annot$Sex, sep = "_"),
                                  levels = c("iWAT_F", "iWAT_M", "gWAT_F", "gWAT_M"))
  ## Depot-only grouping (combining sexes within each depot)
  sample_annot$GroupDepot <- factor(
    paste(sample_annot$Depot, sample_annot$Genotype, sep = "_"),
    levels = c(
      "iWAT_CTL",   "iWAT_KAT8KD",
      "gWAT_CTL",   "gWAT_KAT8KD"
    )
  )
  rownames(sample_annot) <- sample_annot$Sample
  
  cat("\n=== BASIC QC (TISSUE) ===\n")
  cat("Number of genes:   ", nrow(count_matrix_tissue), "\n")
  cat("Number of samples: ", ncol(count_matrix_tissue), "\n")
  cat("Library sizes (colSums raw counts):\n")
  print(colSums(count_matrix_tissue))
  
  ## 4.6) DESeq2 object + typical prefilter ------------------
  ## DESIGN NOTE (2026-01-23): Depot-focused analysis
  ## - Analyzing by depot only (combining males and females within each depot)
  ## - Design ~ 0 + GroupDepot creates separate coefficients for each depot/genotype combination
  ## - Male and female pathways are simple enough; substantial differences are between depots
  ## - This approach captures the primary depot-specific effects of KAT8 knockdown
  min_count <- prefilter_min_count
  min_samples <- prefilter_min_samples

  ## Global filtering: keep genes with >= min_count in >= min_samples anywhere
  keep <- rowSums(count_matrix_tissue >= min_count) >= min_samples
  
  ## Comprehensive filtering sanity check
  n_before <- nrow(count_matrix_tissue)
  n_after <- sum(keep)
  n_removed <- n_before - n_after
  pct_retained <- round(100 * n_after / n_before, 1)
  
  cat("\n==========================================================\n")
  cat("=== SANITY CHECK: GENE PREFILTERING ===\n")
  cat("==========================================================\n")
  cat("Prefilter rule: counts >= ", min_count, " in >= ", min_samples,
      " samples (global filtering)\n", sep = "")
  cat("  (keeps genes expressed anywhere in the experiment)\n")
  cat("  Rationale: Chromatin modifier study - captures context-dependent accessibility\n\n")
  cat("Genes BEFORE prefiltering:  ", n_before, "\n", sep = "")
  cat("Genes AFTER prefiltering:   ", n_after, "\n", sep = "")
  cat("Genes REMOVED:              ", n_removed, " (", 100 - pct_retained, "%)\n", sep = "")
  cat("Genes RETAINED:             ", pct_retained, "%\n", sep = "")
  cat("==========================================================\n\n")
  
  count_matrix_filt <- count_matrix_tissue[keep, , drop = FALSE]
  
  dds_tissue <- DESeqDataSetFromMatrix(
    countData = round(count_matrix_filt),
    colData   = sample_annot,
    design    = ~ 0 + GroupDepot
  )
  
  ## Save centralized parameters
  ## ============================================================
  ## DEPOT-SPECIFIC CUTOFFS (from central parameters.R)
  ## ============================================================
  ## All cutoffs are now defined in parameters.R
  ## iWAT_logFC_cut, iWAT_fdr_cut, gWAT_logFC_cut, gWAT_fdr_cut,
  ## default_logFC_cut, default_fdr_cut are loaded from there
  ## ============================================================

  ## Helper function is now defined in parameters.R as get_tissue_thresholds()
  ## Wrapper for backward compatibility
  get_cutoffs_for_contrast <- function(contrast_name) {
    thresholds <- get_tissue_thresholds(contrast_name)
    return(list(logFC = thresholds$logFC_cut, fdr = thresholds$fdr_cut))
  }

  ## Legacy variables for backward compatibility (use iWAT as default)
  logFC_cut_tissue <- iWAT_logFC_cut
  fdr_cut_tissue   <- iWAT_fdr_cut
  
  params <- list(
    outlier_samples   = outlier_samples,
    prefilter_rule    = list(
      min_count       = min_count,
      min_samples     = min_samples,
      description     = "Keep genes with counts >= min_count in at least min_samples (global filtering for chromatin modifier study)"
    ),
    de_thresholds     = list(
      iWAT = list(logFC_cutoff = iWAT_logFC_cut, fdr_cutoff = iWAT_fdr_cut),
      gWAT = list(logFC_cutoff = gWAT_logFC_cut, fdr_cutoff = gWAT_fdr_cut),
      default = list(logFC_cutoff = default_logFC_cut, fdr_cutoff = default_fdr_cut),
      description = "Depot-specific cutoffs; change iWAT_logFC_cut, gWAT_logFC_cut, etc. to adjust"
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
    paste0("  Description: FDR-corrected p-value threshold for differential expression; logFC threshold filters for biological significance"),
    "",
    "--- VST ---",
    paste0("  blind: ", params$vst$blind),
    paste0("  Description: ", params$vst$description),
    "",
    "--- Volcano Plot Labels ---",
    paste0("  Top N by FDR (per direction): ", params$volcano_labels$top_n_by_fdr),
    paste0("  Top N by |logFC| (per direction): ", params$volcano_labels$top_n_by_fc),
    paste0("  Description: Number of top genes labeled in volcano plots, selected by smallest FDR and largest fold change per direction"),
    "",
    "--- Heatmap ---",
    paste0("  Max genes displayed: ", params$heatmap_max_genes),
    paste0("  Clustering: Hierarchical clustering on rows (genes), split by genotype on columns"),
    paste0("  Scaling: Z-score normalization (row-wise)"),
    paste0("  Color scheme: Teal (down) -> white (neutral) -> orange (up)"),
    "",
    "--- ComplexHeatmap (Genes of Interest) ---",
    paste0("  Gene order: First contrast (iWAT) establishes reference clustering order"),
    paste0("  Description: All subsequent contrasts use same gene order for direct comparison"),
    paste0("  Normalization: Z-score relative to control group means"),
    "",
    "NOTE: This is Part 1 only (gene-level analysis).",
    "For ORA and FGSEA parameters, see part2_ora.R and part3_fgsea.R or comprehensive scripts."
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

  ## Save VST matrix for use by part4 (individual sample heatmaps)
  vst_file <- paste0("VST_matrix_heatmap_", run_tag, ".rds")
  saveRDS(
    list(
      vst_mat = vst_mat_hm,
      sample_info = sample_annot
    ),
    file = file.path(outdir, "tables", vst_file)
  )
  cat("[OK] Saved VST matrix for heatmaps: ", vst_file, "\n", sep = "")

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
    design    = ~ 0 + GroupDepot
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
  ## Depot-specific contrasts (combining sexes within each depot)
  contrast_definitions <- list(
    iWAT_KD_vs_CTL = c("GroupDepot", "iWAT_KAT8KD", "iWAT_CTL"),
    gWAT_KD_vs_CTL = c("GroupDepot", "gWAT_KAT8KD", "gWAT_CTL")
  )
  contrast_names <- names(contrast_definitions)
  
  ## 4.13) Containers ---------------------------------------
  ## DEG_summary: Uses AUTHORITATIVE DEG definition
  ## DEG = gene with padj < fdr_cutoff AND |log2FoldChange| > logFC_cutoff
  ## N_DEG = N_Up + N_Down (no FDR-only counts - that would be misleading)
  deg_summary <- data.frame(
    Contrast     = contrast_names,
    N_DEG        = NA_integer_,
    N_Up         = NA_integer_,
    N_Down       = NA_integer_,
    logFC_cutoff = logFC_cut_tissue,
    fdr_cutoff   = fdr_cut_tissue,
    stringsAsFactors = FALSE
  )

  ## Container for authoritative DEG lists (for downstream ORA/fgsea)
  authoritative_deg_lists <- list()
  
  tt_list <- vector("list", length(contrast_names))
  names(tt_list) <- contrast_names
  
  ## 4.14) Per-contrast DE + volcano + heatmap + genes of interest ----
  ## Initialize reference gene order for genes of interest (will be set from iWAT)
  goi_reference_order <- NULL

  ## Ensure contrasts are processed in order (iWAT must be first for reference)
  contrast_names <- c("iWAT_KD_vs_CTL", "gWAT_KD_vs_CTL")

  for (cn in contrast_names) {

    cat("\n=== Contrast: ", cn, " ===\n", sep = "")

    ## Get depot-specific cutoffs for this contrast
    contrast_cutoffs <- get_cutoffs_for_contrast(cn)
    cn_logFC_cut <- contrast_cutoffs$logFC
    cn_fdr_cut <- contrast_cutoffs$fdr
    cat("[INFO] Using cutoffs: FDR < ", cn_fdr_cut, ", |logFC| > ", cn_logFC_cut, "\n", sep = "")

    res <- results(dds_tissue, contrast = contrast_definitions[[cn]])
    res <- res[order(res$pvalue), , drop = FALSE]

    tt <- as.data.frame(res)
    tt$logFC     <- tt$log2FoldChange
    tt$gene_name <- rownames(tt)

    ## Remove duplicate columns
    tt <- tt[, !duplicated(colnames(tt)), drop = FALSE]

    tt_list[[cn]] <- tt

    de_file <- paste0("DE_tissue_", cn, "_", run_tag, ".csv")
    write.csv(tt, file = file.path(outdir, "tables", de_file), row.names = TRUE)

    ## Calculate comprehensive DEG statistics using depot-specific cutoffs
    n_total <- nrow(tt)
    n_valid_pval <- sum(!is.na(tt$padj))
    n_fdr_only <- sum(tt$padj < cn_fdr_cut, na.rm = TRUE)
    n_up_both <- sum(tt$padj < cn_fdr_cut & tt$logFC > cn_logFC_cut, na.rm = TRUE)
    n_down_both <- sum(tt$padj < cn_fdr_cut & tt$logFC < -cn_logFC_cut, na.rm = TRUE)
    n_both <- n_up_both + n_down_both

    cat("\n==========================================================\n")
    cat("=== SANITY CHECK: ", cn, " ===\n", sep = "")
    cat("==========================================================\n")
    cat("Total genes tested by DESeq2:           ", n_total, "\n", sep = "")
    cat("Genes with valid adjusted p-value:      ", n_valid_pval, "\n", sep = "")
    cat("\n--- DEG Counts at Different Thresholds ---\n")
    cat("DEGs (FDR < ", cn_fdr_cut, " only):               ", n_fdr_only, "\n", sep = "")
    cat("DEGs (FDR < ", cn_fdr_cut, " AND |logFC| > ", cn_logFC_cut, "):  ", n_both, "\n", sep = "")
    cat("  - Up-regulated (logFC > ", cn_logFC_cut, "):       ", n_up_both, "\n", sep = "")
    cat("  - Down-regulated (logFC < -", cn_logFC_cut, "):    ", n_down_both, "\n", sep = "")
    cat("\nPercent of tested genes that are DEGs:\n")
    cat("  - FDR-only threshold:                 ", round(100 * n_fdr_only / n_total, 2), "%\n", sep = "")
    cat("  - FDR + logFC threshold:              ", round(100 * n_both / n_total, 2), "%\n", sep = "")
    cat("==========================================================\n")
    cat("DE table written: ", file.path(outdir, "tables", de_file), "\n", sep = "")
    cat("==========================================================\n\n")

    ## Record AUTHORITATIVE DEG counts (padj + logFC thresholds)
    ## N_DEG = N_Up + N_Down (no FDR-only counts to avoid confusion)
    deg_summary$N_DEG[deg_summary$Contrast == cn] <- n_both
    deg_summary$N_Up[deg_summary$Contrast == cn] <- n_up_both
    deg_summary$N_Down[deg_summary$Contrast == cn] <- n_down_both
    deg_summary$logFC_cutoff[deg_summary$Contrast == cn] <- cn_logFC_cut
    deg_summary$fdr_cutoff[deg_summary$Contrast == cn] <- cn_fdr_cut

    ## Store authoritative DEG lists for downstream ORA/fgsea scripts
    ## These lists include the cutoffs used, so downstream scripts
    ## will automatically use the same parameters and display them in captions
    up_genes <- tt %>%
      dplyr::filter(!is.na(padj), padj < cn_fdr_cut, logFC > cn_logFC_cut) %>%
      rownames()
    down_genes <- tt %>%
      dplyr::filter(!is.na(padj), padj < cn_fdr_cut, logFC < -cn_logFC_cut) %>%
      rownames()

    authoritative_deg_lists[[cn]] <- list(
      up = up_genes,
      down = down_genes,
      logFC_cutoff = cn_logFC_cut,
      fdr_cutoff = cn_fdr_cut
    )
    
    ## Volcano
    tt_plot <- tt %>%
      dplyr::filter(is.finite(logFC), is.finite(padj), padj > 0) %>%
      mutate(negLogFDR = -log10(padj))
    
    tt_plot <- tt_plot %>%
      mutate(
        sig_cat = dplyr::case_when(
          padj < cn_fdr_cut & logFC >=  cn_logFC_cut ~ "Up",
          padj < cn_fdr_cut & logFC <= -cn_logFC_cut ~ "Down",
          TRUE                                            ~ "NS"
        ),
        sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
      )

    tt_sig <- tt_plot %>%
      dplyr::filter(padj < cn_fdr_cut,
                    abs(logFC) >= cn_logFC_cut,
                    sig_cat != "NS")
    
    top_n_fdr <- params$volcano_labels$top_n_by_fdr
    top_n_fc  <- params$volcano_labels$top_n_by_fc
    
    up_sig <- tt_sig %>% dplyr::filter(sig_cat == "Up")
    down_sig <- tt_sig %>% dplyr::filter(sig_cat == "Down")
    
    up_top_fdr <- up_sig %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n_fdr)
    up_top_fc  <- up_sig %>% dplyr::arrange(dplyr::desc(logFC)) %>% dplyr::slice_head(n = top_n_fc)
    
    down_top_fdr <- down_sig %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n_fdr)
    down_top_fc  <- down_sig %>% dplyr::arrange(logFC) %>% dplyr::slice_head(n = top_n_fc)
    
    label_genes <- dplyr::bind_rows(up_top_fdr, up_top_fc, down_top_fdr, down_top_fc) %>%
      dplyr::distinct(gene_name, .keep_all = TRUE)
    
    vol_tissue <- ggplot(
      tt_plot,
      aes(x = logFC, y = negLogFDR, color = sig_cat)
    ) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(
        values = c("Up" = "#FDE725", "Down" = "#440154", "NS" = "grey70"),  ## Viridis yellow/purple
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
      geom_hline(yintercept = -log10(cn_fdr_cut), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(cn_logFC_cut, -cn_logFC_cut), linetype = "dashed", color = "grey40") +
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
      dplyr::filter(padj < cn_fdr_cut, abs(logFC) > cn_logFC_cut)
    
    if (nrow(de_sig) >= 2) {
      de_sig <- de_sig[order(de_sig$padj), , drop = FALSE]
      if (nrow(de_sig) > params$heatmap_max_genes) de_sig <- de_sig[1:params$heatmap_max_genes, , drop = FALSE]
      gene_set <- de_sig$gene_name
      
      depot <- NA_character_
      if (grepl("^iWAT", cn))      { depot <- "iWAT" }
      else if (grepl("^gWAT", cn)) { depot <- "gWAT" }

      if (!is.na(depot)) {

        idx_samples <- with(as.data.frame(colData(dds_tissue)), Depot == depot)
        samples_heat <- rownames(colData(dds_tissue))[idx_samples]
        heat_meta_deg <- as.data.frame(colData(dds_tissue))[idx_samples, , drop = FALSE]
        
        mat <- vst_mat_hm[gene_set, samples_heat, drop = FALSE]
        
        ## Scale rows (z-score)
        mat_scaled <- t(scale(t(mat)))
        mat_scaled <- mat_scaled[!rowSums(is.na(mat_scaled)), , drop = FALSE]
        
        if (nrow(mat_scaled) >= 2) {
          ## Use viridis-inspired diverging colors (purple -> white -> yellow)
          max_val_deg <- max(abs(mat_scaled), na.rm = TRUE)
          deg_col_fun <- colorRamp2(
            c(-max_val_deg, 0, max_val_deg),
            c("#440154", "white", "#FDE725")  ## Purple (down) -> white -> yellow (up)
          )
          
          ## Column annotation (including Sex since depots combine both sexes)
          deg_ha <- HeatmapAnnotation(
            Sex = heat_meta_deg$Sex,
            Genotype = heat_meta_deg$Genotype,
            col = list(
              Sex = c(F = "#E7298A", M = "#7570B3"),
              Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")
            ),
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
            column_title = paste0("Top DEGs: ", cn, " (FDR<", cn_fdr_cut, ", |logFC|>", cn_logFC_cut, ")"),
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
      
      ## Get depot for this contrast
      depot_h <- NA_character_
      if (grepl("^iWAT", cn))      { depot_h <- "iWAT" }
      else if (grepl("^gWAT", cn)) { depot_h <- "gWAT" }

      if (!is.na(depot_h)) {

        ## Get sample indices for this contrast
        heat_idx <- with(as.data.frame(colData(dds_tissue)), Depot == depot_h)
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

          ## Use viridis-inspired diverging colors (purple -> white -> yellow)
          max_val <- max(abs(heat_matrix_scaled), na.rm = TRUE)
          heat_col_fun <- colorRamp2(
            c(-max_val, 0, max_val),
            c("#440154", "white", "#FDE725")  ## Purple (down) -> white -> yellow (up)
          )

          ## Column annotation (including Sex since depots combine both sexes)
          heat_ha <- HeatmapAnnotation(
            Sex = heat_meta$Sex,
            Genotype = heat_meta$Genotype,
            col = list(
              Sex = c(F = "#E7298A", M = "#7570B3"),
              Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02")
            ),
            annotation_name_side = "left"
          )

          ## Establish or apply reference gene order (iWAT is reference)
          if (cn == "iWAT_KD_vs_CTL" && is.null(goi_reference_order)) {
            ## First contrast: cluster rows and save order as gene names
            cat("[INFO] Establishing reference gene order from iWAT...\n")
            ## Perform hierarchical clustering manually to get order
            row_dist <- dist(heat_matrix_scaled)
            row_hclust <- hclust(row_dist, method = "complete")
            goi_reference_order <- rownames(heat_matrix_scaled)[row_hclust$order]
            cat("[INFO] Reference order established: ", length(goi_reference_order), " genes\n", sep = "")

            ## Reorder matrix to match clustered order
            heat_matrix_scaled <- heat_matrix_scaled[goi_reference_order, , drop = FALSE]

            ## Print gene order for verification
            cat("[INFO] Gene order (first 10): ", paste(head(rownames(heat_matrix_scaled), 10), collapse=", "), "\n", sep="")

            ## Create heatmap WITHOUT clustering (using pre-clustered order)
            heat_goi <- Heatmap(
              heat_matrix_scaled,
              col = heat_col_fun,
              cluster_rows = FALSE,  ## No clustering - already in clustered order
              cluster_columns = FALSE,
              show_column_names = TRUE,
              show_row_names = TRUE,
              column_split = heat_meta$Genotype,
              border = TRUE,
              column_gap = unit(2, "mm"),
              row_gap = unit(0, "mm"),
              top_annotation = heat_ha,
              column_title = paste0("Genes of Interest: ", cn, " (REFERENCE order)"),
              heatmap_legend_param = list(
                title = "Z-score\n(vs CTL)",
                title_position = "leftcenter-rot",
                title_gp = gpar(fontsize = 10)
              ),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
            )
          } else if (!is.null(goi_reference_order)) {
            ## Subsequent contrasts: use reference order (no clustering)
            cat("[INFO] Applying reference gene order from iWAT...\n")
            ## Reorder matrix to match reference (only use genes present in both)
            common_genes <- intersect(goi_reference_order, rownames(heat_matrix_scaled))
            heat_matrix_scaled <- heat_matrix_scaled[common_genes, , drop = FALSE]
            cat("[INFO] Using ", length(common_genes), " genes in common with reference\n", sep = "")
            cat("[INFO] Gene order (first 10): ", paste(head(rownames(heat_matrix_scaled), 10), collapse=", "), "\n", sep="")

            heat_goi <- Heatmap(
              heat_matrix_scaled,
              col = heat_col_fun,
              cluster_rows = FALSE,  ## No clustering - use reference order
              cluster_columns = FALSE,
              show_column_names = TRUE,
              show_row_names = TRUE,
              column_split = heat_meta$Genotype,
              border = TRUE,
              column_gap = unit(2, "mm"),
              row_gap = unit(0, "mm"),
              top_annotation = heat_ha,
              column_title = paste0("Genes of Interest: ", cn, " (ordered by iWAT)"),
              heatmap_legend_param = list(
                title = "Z-score\n(vs CTL)",
                title_position = "leftcenter-rot",
                title_gp = gpar(fontsize = 10)
              ),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8)
            )
          } else {
            ## Fallback: cluster normally if reference not established yet
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
          }
          
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
      dplyr::filter(is.finite(logFC), is.finite(padj), padj > 0)
    
    ## Only label genes of interest that are SIGNIFICANT (pass both thresholds)
    goi_significant <- volcano_df %>%
      dplyr::filter(gene_name %in% genes_of_interest,
                    padj < cn_fdr_cut,
                    abs(logFC) >= cn_logFC_cut) %>%
      dplyr::pull(gene_name)
    
    cat("[INFO] Genes of interest (total): ", length(genes_of_interest), "\n", sep = "")
    cat("[INFO] Genes of interest (significant): ", length(goi_significant), "\n", sep = "")
    
    if (nrow(volcano_df) > 0) {
      
      ## Create custom color scheme using keyvals (viridis-inspired)
      keyvals_color <- rep("grey70", nrow(volcano_df))
      names(keyvals_color) <- rep("NS", nrow(volcano_df))

      ## Significant UP (viridis yellow)
      sig_up_idx <- which(volcano_df$padj < cn_fdr_cut &
                            volcano_df$logFC >= cn_logFC_cut)
      keyvals_color[sig_up_idx] <- "#FDE725"
      names(keyvals_color)[sig_up_idx] <- "Significant Up"

      ## Significant DOWN (viridis purple)
      sig_down_idx <- which(volcano_df$padj < cn_fdr_cut &
                              volcano_df$logFC <= -cn_logFC_cut)
      keyvals_color[sig_down_idx] <- "#440154"
      names(keyvals_color)[sig_down_idx] <- "Significant Down"
      
      ev_plot <- EnhancedVolcano(
        volcano_df,
        lab = volcano_df$gene_name,
        selectLab = goi_significant,
        x = "logFC",
        y = "padj",
        title = paste0("Volcano: ", cn),
        subtitle = "ECM, metabolism, and remodeling genes labeled (significant only)",
        pCutoff = cn_fdr_cut,
        FCcutoff = cn_logFC_cut,
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
  

  ## 4.15) DEG summary table --------------------------------
  deg_summary_file <- paste0("DEG_summary_tissue_", run_tag, ".csv")
  write.csv(deg_summary, file = file.path(outdir, "tables", deg_summary_file), row.names = FALSE)
  cat("\nDEG summary written: ", file.path(outdir, "tables", deg_summary_file), "\n", sep = "")

  ## 4.15b) Export authoritative DEG lists for downstream ORA/fgsea ----
  ## This ensures ORA and fgsea use IDENTICAL gene sets (no drift)
  deg_lists_file <- paste0("DEG_lists_authoritative_", run_tag, ".rds")
  saveRDS(authoritative_deg_lists, file = file.path(outdir, "tables", deg_lists_file))
  cat("Authoritative DEG lists saved: ", file.path(outdir, "tables", deg_lists_file), "\n", sep = "")

  ## Also print validation summary
  cat("\n==========================================================\n")
  cat("=== VALIDATION: AUTHORITATIVE DEG DEFINITION ===\n")
  cat("==========================================================\n")
  cat("DEG definition: padj < fdr_cutoff AND |log2FC| > logFC_cutoff\n")
  cat("(Cutoffs are depot-specific - see per-contrast details below)\n")
  cat("\nSUMMARY CHECK: N_DEG = N_Up + N_Down\n")
  for (cn in contrast_names) {
    n_deg <- deg_summary$N_DEG[deg_summary$Contrast == cn]
    n_up <- deg_summary$N_Up[deg_summary$Contrast == cn]
    n_down <- deg_summary$N_Down[deg_summary$Contrast == cn]
    cn_logfc <- deg_summary$logFC_cutoff[deg_summary$Contrast == cn]
    cn_fdr <- deg_summary$fdr_cutoff[deg_summary$Contrast == cn]
    check <- ifelse(n_deg == n_up + n_down, "PASS", "FAIL")
    cat("  ", cn, ":\n", sep = "")
    cat("    Cutoffs: FDR < ", cn_fdr, ", |logFC| > ", cn_logfc, "\n", sep = "")
    cat("    N_DEG=", n_deg, " (Up=", n_up, " + Down=", n_down, ") [", check, "]\n", sep = "")
  }
  cat("==========================================================\n")
  
  ## 4.16) save_core bookkeeping ----------------------------
  if (use_save_core && !is.null(hero_volcano_file)) {
    hero_path <- file.path(outdir, "plots", hero_volcano_file)
    if (file.exists(hero_path) && exists("save_run_file")) {
      tryCatch({
        save_run_file(hero_path, tag = "hero_volcano")
      }, error = function(e) {
        cat("[WARN] save_run_file failed: ", conditionMessage(e), "\n", sep = "")
      })
    }
  }

  if (use_save_core && exists("dedupe")) {
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
