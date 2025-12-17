#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq: Tissue KAT8KD vs CTL
## iWAT & gWAT, F/M; FL vs AKO/MAKO
##
## COMPLETE PUBLICATION-READY PIPELINE:
## - Quality control: Density plots, PCA, MDS
## - Differential expression: limma-voom
## - Visualizations: Volcano plots, heatmaps
## - Pathway analysis: Hallmark, KEGG, GO:BP
## - Publication-ready figures (300 DPI)
## - Comprehensive pathway bar plots (up/down regulated)
## - Session info for reproducibility
##
## WITH save_core / run management support
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c(
  "dplyr",
  "limma",
  "edgeR",
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "viridis",
  "pheatmap",
  "grid",
  "scales"
)

## Pathway analysis packages (Bioconductor)
required_bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "DOSE",
  "pathview",
  "msigdbr"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

## Install Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(RColorBrewer)
  library(ggrepel)
  library(viridis)
  library(pheatmap)
  library(grid)
  library(scales)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(msigdbr)
})

cat("[OK] All packages loaded successfully\n")

## 2) save_core utilities -----------------------------------

utils_dir <- file.path(getwd(), "save_core")
if (dir.exists(utils_dir)) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))

  source(file.path(utils_dir, "purge.R"))
  source(file.path(utils_dir, "dedupe.R"))
  cat("[OK] save_core utilities loaded\n")
}

## Always define fallback functions if they don't exist after sourcing
if (!exists("init_run")) {
  cat("[INFO] Using basic output structure (init_run not found)\n")
  init_run <- function(script_name, species, data_type, keywords, notes) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_tag <- timestamp
    outdir <- file.path(getwd(), "output", paste0("RUN_", run_tag))
    dir.create(file.path(outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(outdir, "logs"), recursive = TRUE, showWarnings = FALSE)
    list(outdir = outdir, run_tag = run_tag)
  }
}

if (!exists("save_run_file")) {
  save_run_file <- function(file_path, tag = NULL) {
    invisible(NULL)
  }
}

if (!exists("dedupe")) {
  dedupe <- function(outdir) {
    invisible(NULL)
  }
}

## 2.5) Run metadata ----------------------------------------

run_info <- list(
  script_name = "KAT8_bulk_tissue_KAT8KD.R",
  species     = "mouse",
  data_type   = "bulkRNAseq_WAT_tissue",
  keywords    = c("KAT8", "WAT", "iWAT", "gWAT", "bulk RNA-seq", "tissue"),
  notes       = "KAT8KD vs CTL in iWAT/gWAT, F/M; per-depot/sex contrasts",
  message     = "Tissue KAT8KD vs CTL: limma-voom + PCA/MDS + DEGs + Pathway Analysis"
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

## 3.5) Load pathway gene sets -----------------------------

cat("\n=== LOADING PATHWAY GENE SETS ===\n")

## MSigDB Hallmark gene sets for mouse
msigdb_hallmark <- msigdbr(species = "Mus musculus", category = "H")

## MSigDB C2 KEGG gene sets - use robust loading with error handling
msigdb_kegg <- tryCatch({
  msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
}, error = function(e) {
  cat("[WARN] CP:KEGG subcategory not found, trying KEGG legacy format...\n")
  tryCatch({
    msigdbr(species = "Mus musculus", category = "C2", subcategory = "KEGG")
  }, error = function(e2) {
    cat("[WARN] Loading all C2 and filtering for KEGG...\n")
    all_c2 <- msigdbr(species = "Mus musculus", category = "C2")
    all_c2 %>% dplyr::filter(grepl("KEGG", gs_name, ignore.case = TRUE))
  })
})

## MSigDB C5 GO Biological Process - use robust loading
msigdb_gobp <- tryCatch({
  msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
}, error = function(e) {
  cat("[WARN] GO:BP subcategory syntax issue, trying alternative...\n")
  tryCatch({
    msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
  }, error = function(e2) {
    cat("[WARN] Loading all C5 and filtering for GO:BP...\n")
    all_c5 <- msigdbr(species = "Mus musculus", category = "C5")
    all_c5 %>% dplyr::filter(gs_subcat == "GO:BP" | grepl("^GO.*BIOLOGICAL.*PROCESS", gs_name, ignore.case = TRUE))
  })
})

cat("[OK] Loaded", length(unique(msigdb_hallmark$gs_name)), "Hallmark pathways\n")
cat("[OK] Loaded", length(unique(msigdb_kegg$gs_name)), "KEGG pathways\n")
cat("[OK] Loaded", length(unique(msigdb_gobp$gs_name)), "GO:BP terms\n")

## Validate we have gene sets
if (length(unique(msigdb_hallmark$gs_name)) == 0) {
  warning("No Hallmark pathways loaded - pathway analysis may be incomplete")
}
if (length(unique(msigdb_kegg$gs_name)) == 0) {
  warning("No KEGG pathways loaded - pathway analysis may be incomplete")
}
if (length(unique(msigdb_gobp$gs_name)) == 0) {
  warning("No GO:BP terms loaded - pathway analysis may be incomplete")
}

## 4) Main logic --------------------------------------------

hero_volcano_file <- NULL

tryCatch({

  ## 4.1) Load raw counts -----------------------------------

  cat("\n=== RAW COUNTS MATRIX (FULL) ===\n")

  ## Check if counts file exists
  counts_file <- "counts.txt"
  if (!file.exists(counts_file)) {
    stop("ERROR: counts.txt file not found in working directory: ", getwd(),
         "\nPlease ensure your counts file is present.")
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

  ## ---- SANITY CHECK 1: Validate counts file structure ----
  cat("\n--- SANITY CHECK 1: Counts file structure ---\n")

  required_cols <- c("gene_id", "gene")
  missing_required <- setdiff(required_cols, colnames(counts_raw))
  if (length(missing_required) > 0) {
    stop("ERROR: Missing required columns in counts.txt: ",
         paste(missing_required, collapse = ", "))
  }
  cat("[PASS] Required columns 'gene_id' and 'gene' found\n")

  ## Check for numeric count columns
  sample_cols <- setdiff(colnames(counts_raw), c("gene_id", "gene"))
  cat("[INFO] Found", length(sample_cols), "sample columns\n")

  ## Check for negative values
  count_matrix_check <- as.matrix(counts_raw[, sample_cols])
  if (any(count_matrix_check < 0, na.rm = TRUE)) {
    warning("WARNING: Negative values found in counts matrix!")
    cat("[WARN] Negative values detected - this is unusual for RNA-seq counts\n")
  } else {
    cat("[PASS] All count values are non-negative\n")
  }

  ## Check for NA values
  na_count <- sum(is.na(count_matrix_check))
  if (na_count > 0) {
    warning("WARNING: ", na_count, " NA values found in counts matrix!")
    cat("[WARN]", na_count, "NA values detected in counts\n")
  } else {
    cat("[PASS] No NA values in counts matrix\n")
  }

  ## Check gene_id uniqueness
  if (any(duplicated(counts_raw$gene_id))) {
    n_dup <- sum(duplicated(counts_raw$gene_id))
    cat("[WARN]", n_dup, "duplicated gene_id entries found\n")
  } else {
    cat("[PASS] All gene_id entries are unique\n")
  }

  ## 4.2) Define tissue samples + annotation ----------------

  sample_ids <- paste0("JS_", sprintf("%02d", 1:40))

  sample_annot_full <- data.frame(
    Sample   = sample_ids,
    Depot    = c(
      rep("iWAT", 10),
      rep("iWAT", 10),
      rep("gWAT", 10),
      rep("gWAT", 10)
    ),
    Sex      = c(
      rep("F", 5),  rep("F", 5),
      rep("M", 5),  rep("M", 5),
      rep("F", 5),  rep("F", 5),
      rep("M", 5),  rep("M", 5)
    ),
    Genotype = c(
      rep("CTL",    5), rep("KAT8KD", 5),
      rep("CTL",    5), rep("KAT8KD", 5),
      rep("CTL",    5), rep("KAT8KD", 5),
      rep("CTL",    5), rep("KAT8KD", 5)
    ),
    stringsAsFactors = FALSE
  )

  ## ---- SANITY CHECK 2: Sample annotation validation ----
  cat("\n--- SANITY CHECK 2: Sample annotation ---\n")

  ## Check expected vs found samples
  expected_samples <- sample_annot_full$Sample
  found_samples <- intersect(expected_samples, colnames(counts_raw))
  missing_samples <- setdiff(expected_samples, colnames(counts_raw))

  cat("[INFO] Expected samples:", length(expected_samples), "\n")
  cat("[INFO] Found in counts:", length(found_samples), "\n")

  if (length(missing_samples) > 0) {
    cat("[WARN] Missing samples:", paste(missing_samples, collapse = ", "), "\n")
  } else {
    cat("[PASS] All expected samples found in counts file\n")
  }

  ## Verify annotation structure
  cat("[INFO] Annotation structure:\n")
  cat("  - Depots:", paste(unique(sample_annot_full$Depot), collapse = ", "), "\n")
  cat("  - Sexes:", paste(unique(sample_annot_full$Sex), collapse = ", "), "\n")
  cat("  - Genotypes:", paste(unique(sample_annot_full$Genotype), collapse = ", "), "\n")

  ## Check sample counts per group
  cat("[INFO] Samples per group (before outlier removal):\n")
  group_counts <- sample_annot_full %>%
    dplyr::group_by(Depot, Sex, Genotype) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  print(as.data.frame(group_counts))

  ## =========================================================
  ## 4.2.1) QC PLOTS BEFORE OUTLIER REMOVAL
  ## Evidence supporting outlier removal decision
  ## =========================================================

  cat("\n=== QC BEFORE OUTLIER REMOVAL ===\n")

  ## Get counts for ALL tissue samples (including outliers)
  all_tissue_samples <- sample_annot_full$Sample
  missing_all <- setdiff(all_tissue_samples, colnames(counts_raw))
  if (length(missing_all) > 0) {
    stop("Missing columns for QC: ", paste(missing_all, collapse = ", "))
  }

  counts_all_tissue <- as.matrix(counts_raw[, all_tissue_samples])
  rownames(counts_all_tissue) <- counts_raw$gene_id

  ## Create temp DGEList with ALL samples
  DGE_all <- DGEList(counts = counts_all_tissue)
  DGE_all <- calcNormFactors(DGE_all)

  ## Attach sample info
  DGE_all$samples$Sample   <- sample_annot_full$Sample
  DGE_all$samples$Depot    <- sample_annot_full$Depot
  DGE_all$samples$Sex      <- sample_annot_full$Sex
  DGE_all$samples$Genotype <- factor(sample_annot_full$Genotype,
                                      levels = c("CTL", "KAT8KD"))
  DGE_all$samples$DepotSex <- factor(
    paste(sample_annot_full$Depot, sample_annot_full$Sex, sep = "_"),
    levels = c("iWAT_F", "iWAT_M", "gWAT_F", "gWAT_M")
  )

  ## Define outlier samples
  outlier_samples <- c("JS_08", "JS_28")
  DGE_all$samples$IsOutlier <- DGE_all$samples$Sample %in% outlier_samples

  ## Compute logCPM for all samples
  logCPM_all <- cpm(DGE_all, log = TRUE, prior.count = 3)

  ## --- (a) Library Size Plot ---
  lib_sizes_all <- data.frame(
    Sample     = colnames(DGE_all$counts),
    LibSize    = DGE_all$samples$lib.size,
    DepotSex   = DGE_all$samples$DepotSex,
    Genotype   = DGE_all$samples$Genotype,
    IsOutlier  = DGE_all$samples$IsOutlier,
    stringsAsFactors = FALSE
  )
  lib_sizes_all$Sample <- factor(lib_sizes_all$Sample,
                                  levels = lib_sizes_all$Sample[order(lib_sizes_all$LibSize)])

  qc_libsize_plot <- ggplot(lib_sizes_all,
                             aes(x = Sample, y = LibSize / 1e6,
                                 fill = DepotSex, alpha = IsOutlier)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_text(
      data = lib_sizes_all[lib_sizes_all$IsOutlier, ],
      aes(label = Sample),
      vjust = -0.5, hjust = 0.5, size = 3, fontface = "bold", color = "red"
    ) +
    scale_fill_manual(
      values = c("iWAT_F" = "#1b9e77", "iWAT_M" = "#d95f02",
                 "gWAT_F" = "#7570b3", "gWAT_M" = "#e7298a"),
      name = "Depot/Sex"
    ) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.5),
                       name = "Outlier", labels = c("No", "Yes")) +
    labs(
      title    = "Library Sizes (Before Outlier Removal)",
      subtitle = paste0("Outliers highlighted: ", paste(outlier_samples, collapse = ", ")),
      x        = "Sample",
      y        = "Library Size (millions)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5, color = "red"),
      legend.position = "right"
    )

  libsize_file <- paste0("QC_library_size_before_outlier_removal_", run_tag, ".png")
  ggsave(
    file.path(outdir, "plots", libsize_file),
    plot   = qc_libsize_plot,
    width  = 12,
    height = 6,
    dpi    = 300
  )
  cat("Saved QC library size plot to",
      file.path(outdir, "plots", libsize_file), "\n")

  ## --- (b) Sample-Sample Correlation Heatmap ---
  cor_mat <- cor(logCPM_all, method = "pearson")

  annot_row_corr <- data.frame(
    DepotSex = DGE_all$samples$DepotSex,
    Genotype = DGE_all$samples$Genotype,
    Outlier  = ifelse(DGE_all$samples$IsOutlier, "Yes", "No")
  )
  rownames(annot_row_corr) <- colnames(DGE_all$counts)

  annot_colors_corr <- list(
    DepotSex = c("iWAT_F" = "#1b9e77", "iWAT_M" = "#d95f02",
                 "gWAT_F" = "#7570b3", "gWAT_M" = "#e7298a"),
    Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02"),
    Outlier  = c(Yes = "red", No = "grey80")
  )

  corr_heatmap_file <- paste0("QC_sample_correlation_before_outlier_removal_", run_tag, ".png")
  png(file.path(outdir, "plots", corr_heatmap_file),
      width = 2400, height = 2200, res = 200)
  pheatmap(
    cor_mat,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks            = seq(0.7, 1, length.out = 101),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 8,
    annotation_row    = annot_row_corr,
    annotation_col    = annot_row_corr,
    annotation_colors = annot_colors_corr,
    main              = paste0("Sample Correlation (logCPM, before outlier removal)\n",
                               "Outliers: ", paste(outlier_samples, collapse = ", "))
  )
  dev.off()
  cat("Saved QC sample correlation heatmap to",
      file.path(outdir, "plots", corr_heatmap_file), "\n")

  ## --- (c) PCA Before Outlier Removal ---
  pca_all <- prcomp(t(logCPM_all), scale. = TRUE)
  pca_var_all <- summary(pca_all)$importance[2, 1:2] * 100

  pca_all_df <- data.frame(
    Sample    = colnames(DGE_all$counts),
    Genotype  = DGE_all$samples$Genotype,
    DepotSex  = DGE_all$samples$DepotSex,
    Depot     = DGE_all$samples$Depot,
    Sex       = DGE_all$samples$Sex,
    IsOutlier = DGE_all$samples$IsOutlier,
    PC1       = pca_all$x[, 1],
    PC2       = pca_all$x[, 2],
    stringsAsFactors = FALSE
  )

  qc_pca_plot <- ggplot(
    pca_all_df,
    aes(x = PC1, y = PC2,
        color = DepotSex,
        shape = Genotype,
        label = Sample)
  ) +
    geom_point(aes(size = IsOutlier), alpha = 0.85) +
    geom_point(
      data = pca_all_df[pca_all_df$IsOutlier, ],
      color = "red", shape = 1, size = 6, stroke = 2
    ) +
    geom_text_repel(
      data          = pca_all_df[pca_all_df$IsOutlier, ],
      aes(label     = Sample),
      color         = "red",
      size          = 4,
      fontface      = "bold",
      box.padding   = unit(0.5, "lines"),
      max.overlaps  = Inf
    ) +
    geom_text_repel(
      data          = pca_all_df[!pca_all_df$IsOutlier, ],
      aes(label     = Sample),
      size          = 2.5,
      color         = "gray30",
      max.overlaps  = 30
    ) +
    scale_color_manual(
      values = c("iWAT_F" = "#1b9e77", "iWAT_M" = "#d95f02",
                 "gWAT_F" = "#7570b3", "gWAT_M" = "#e7298a"),
      name = "Depot/Sex"
    ) +
    scale_shape_manual(
      values = c("CTL" = 16, "KAT8KD" = 17),
      name = "Genotype"
    ) +
    scale_size_manual(
      values = c("FALSE" = 3, "TRUE" = 4),
      guide = "none"
    ) +
    labs(
      title    = "PCA Before Outlier Removal",
      subtitle = paste0("Outliers circled in red: ", paste(outlier_samples, collapse = ", ")),
      x        = paste0("PC1 (", round(pca_var_all[1], 1), "%)"),
      y        = paste0("PC2 (", round(pca_var_all[2], 1), "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5, color = "red"),
      legend.position = "right"
    )

  pca_qc_file <- paste0("QC_PCA_before_outlier_removal_", run_tag, ".png")
  ggsave(
    file.path(outdir, "plots", pca_qc_file),
    plot   = qc_pca_plot,
    width  = 10,
    height = 8,
    dpi    = 300
  )
  cat("Saved QC PCA plot to",
      file.path(outdir, "plots", pca_qc_file), "\n")

  ## --- (d) Outlier Rationale Log File ---
  outlier_log_file <- file.path(outdir, "logs", paste0("outliers_", run_tag, ".txt"))

  lib_summary <- summary(DGE_all$samples$lib.size)
  outlier_lib_sizes <- DGE_all$samples$lib.size[DGE_all$samples$IsOutlier]
  names(outlier_lib_sizes) <- DGE_all$samples$Sample[DGE_all$samples$IsOutlier]

  outlier_log_lines <- c(
    "=== Outlier Removal Rationale ===",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Outlier Sample IDs ---",
    paste0("  ", paste(outlier_samples, collapse = ", ")),
    "",
    "--- Library Sizes (All Samples) ---",
    paste0("  Min:     ", format(lib_summary["Min."], big.mark = ",")),
    paste0("  1st Qu:  ", format(lib_summary["1st Qu."], big.mark = ",")),
    paste0("  Median:  ", format(lib_summary["Median"], big.mark = ",")),
    paste0("  Mean:    ", format(lib_summary["Mean"], big.mark = ",")),
    paste0("  3rd Qu:  ", format(lib_summary["3rd Qu."], big.mark = ",")),
    paste0("  Max:     ", format(lib_summary["Max."], big.mark = ",")),
    "",
    "--- Outlier Library Sizes ---",
    sapply(names(outlier_lib_sizes), function(s) {
      paste0("  ", s, ": ", format(outlier_lib_sizes[s], big.mark = ","))
    }),
    "",
    "--- Rationale ---",
    "  Flagged based on separation in PCA/MDS and/or low library size in QC plots.",
    "  See QC plots for visual evidence:",
    paste0("    - ", libsize_file),
    paste0("    - ", corr_heatmap_file),
    paste0("    - ", pca_qc_file)
  )
  writeLines(outlier_log_lines, con = outlier_log_file)
  cat("Saved outlier rationale log to", outlier_log_file, "\n")

  cat("=== END QC BEFORE OUTLIER REMOVAL ===\n\n")

  ## Clean up QC objects to free memory
  rm(DGE_all, logCPM_all, cor_mat, pca_all)
  gc()

  ## =========================================================
  ## Continue with outlier removal
  ## =========================================================

  ## Remove outlier samples at the annotation level
  ## NOTE: JS_08 and JS_28 identified as outliers from preliminary QC
  ## (extreme deviation in PCA/MDS, low library size, or sample quality issues)
  cat("\n[INFO] Removing outlier samples:", paste(outlier_samples, collapse = ", "), "\n")
  sample_annot <- sample_annot_full %>%
    dplyr::filter(!(Sample %in% outlier_samples))

  tissue_samples <- sample_annot$Sample

  ## Sanity check
  missing_tissue <- setdiff(tissue_samples, colnames(counts_raw))
  if (length(missing_tissue) > 0) {
    stop(
      "These tissue sample columns are missing in counts.txt: ",
      paste(missing_tissue, collapse = ", ")
    )
  }

  ## Check sample counts per group after outlier removal
  cat("[INFO] Samples per group (after outlier removal):\n")
  group_counts_after <- sample_annot %>%
    dplyr::group_by(Depot, Sex, Genotype) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
  print(as.data.frame(group_counts_after))

  ## Warn if any group has fewer than 3 samples
  min_group_size <- min(group_counts_after$n)
  if (min_group_size < 3) {
    cat("[WARN] Some groups have fewer than 3 samples - statistical power may be limited\n")
  } else {
    cat("[PASS] All groups have at least 3 samples (min =", min_group_size, ")\n")
  }

  counts_tissue_raw <- counts_raw[, c("gene_id", "gene", tissue_samples)]

  cat("\n=== TISSUE SUBSET (with outliers removed) ===\n")
  cat("Dimensions (rows x cols):", paste(dim(counts_tissue_raw), collapse = " x "), "\n")
  cat("Column names:\n")
  print(colnames(counts_tissue_raw))

  ## 4.3) Gene cleaning -------------------------------------

  counts_tissue <- counts_tissue_raw %>%
    mutate(
      ensembl_gene_id = gene_id,
      gene_name       = gene
    )

  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]

  cat("\nAny NA in gene_name (tissue)?  ", any(is.na(counts_tissue$gene_name)), "\n")
  cat("Any blank gene_name (tissue)?  ", any(counts_tissue$gene_name == ""), "\n")
  cat("Any duplicated gene_name (tissue)?  ", any(duplicated(counts_tissue$gene_name)), "\n")

  dup_count_tissue <- sum(duplicated(counts_tissue$gene_name))
  cat("Number of duplicated gene_name entries (tissue):", dup_count_tissue, "\n")

  if (dup_count_tissue > 0) {
    cat("Renaming duplicated gene_name entries as ensembl_gene_id_gene_name (tissue).\n")

    dup_flag <- duplicated(counts_tissue$gene_name) |
      duplicated(counts_tissue$gene_name, fromLast = TRUE)

    counts_tissue$gene_name[dup_flag] <- paste0(
      counts_tissue$ensembl_gene_id[dup_flag],
      "_",
      counts_tissue$gene_name[dup_flag]
    )
  }

  rownames(counts_tissue) <- counts_tissue$gene_name

  cat("\n=== AFTER GENE CLEANING (TISSUE) ===\n")
  cat("Dimensions (rows x cols):", paste(dim(counts_tissue), collapse = " x "), "\n")

  ## 4.4) DGEList construction ------------------------------

  count_matrix_tissue <- as.matrix(counts_tissue[, tissue_samples])
  rownames(count_matrix_tissue) <- counts_tissue$gene_name

  DGE_tissue <- DGEList(counts = count_matrix_tissue)

  DGE_tissue$genes <- data.frame(
    gene_name       = rownames(DGE_tissue$counts),
    ensembl_gene_id = counts_tissue$ensembl_gene_id[
      match(rownames(DGE_tissue$counts), counts_tissue$gene_name)
    ],
    stringsAsFactors = FALSE
  )

  ## 4.5) Attach sample annotation --------------------------

  sample_annot <- sample_annot[match(colnames(DGE_tissue$counts),
                                     sample_annot$Sample), ]

  DGE_tissue$samples$Sample   <- sample_annot$Sample
  DGE_tissue$samples$Depot    <- sample_annot$Depot
  DGE_tissue$samples$Sex      <- sample_annot$Sex
  DGE_tissue$samples$Genotype <- factor(sample_annot$Genotype,
                                        levels = c("CTL", "KAT8KD"))

  DGE_tissue$samples$DepotSex <- factor(
    paste(DGE_tissue$samples$Depot, DGE_tissue$samples$Sex, sep = "_"),
    levels = c("iWAT_F", "iWAT_M", "gWAT_F", "gWAT_M")
  )

  DGE_tissue$samples$GroupFull <- factor(
    paste(DGE_tissue$samples$Depot,
          DGE_tissue$samples$Sex,
          DGE_tissue$samples$Genotype,
          sep = "_"),
    levels = c(
      "iWAT_F_CTL",   "iWAT_F_KAT8KD",
      "iWAT_M_CTL",   "iWAT_M_KAT8KD",
      "gWAT_F_CTL",   "gWAT_F_KAT8KD",
      "gWAT_M_CTL",   "gWAT_M_KAT8KD"
    )
  )

  rownames(DGE_tissue$samples) <- colnames(DGE_tissue$counts)

  cat("\nSample annotation head:\n")
  print(head(DGE_tissue$samples))

  cat("\n=== DGE BASIC QC (TISSUE) ===\n")
  cat("Number of genes:  ", nrow(DGE_tissue$counts), "\n")
  cat("Number of samples:", ncol(DGE_tissue$counts), "\n")
  cat("Library sizes:\n")
  print(DGE_tissue$samples$lib.size)
  cat("\nSummary of library sizes:\n")
  print(summary(DGE_tissue$samples$lib.size))

  ## 4.6) TMM normalization + filtering ---------------------

  DGE_tissue <- calcNormFactors(DGE_tissue)

  ## ---- SANITY CHECK 3: Library sizes and normalization factors ----
  cat("\n--- SANITY CHECK 3: Library sizes and normalization ---\n")

  lib_sizes <- DGE_tissue$samples$lib.size
  norm_factors <- DGE_tissue$samples$norm.factors

  cat("[INFO] Library size summary:\n")
  print(summary(lib_sizes))

  cat("\n[INFO] Normalization factors:\n")
  print(data.frame(
    Sample = rownames(DGE_tissue$samples),
    LibSize = format(lib_sizes, big.mark = ","),
    NormFactor = round(norm_factors, 4)
  ))

  ## Check for abnormally low library sizes (< 1 million reads)
  low_lib_samples <- rownames(DGE_tissue$samples)[lib_sizes < 1e6]
  if (length(low_lib_samples) > 0) {
    cat("[WARN] Samples with low library size (<1M reads):",
        paste(low_lib_samples, collapse = ", "), "\n")
  } else {
    cat("[PASS] All samples have library size >= 1M reads\n")
  }

  ## Check normalization factors (should be close to 1, typically 0.8-1.2)
  extreme_norm <- rownames(DGE_tissue$samples)[norm_factors < 0.7 | norm_factors > 1.3]
  if (length(extreme_norm) > 0) {
    cat("[WARN] Samples with extreme normalization factors (outside 0.7-1.3):",
        paste(extreme_norm, collapse = ", "), "\n")
  } else {
    cat("[PASS] All normalization factors within expected range (0.7-1.3)\n")
  }

  ## Calculate coefficient of variation for library sizes
  lib_cv <- sd(lib_sizes) / mean(lib_sizes) * 100
  cat("[INFO] Library size CV:", round(lib_cv, 1), "%\n")
  if (lib_cv > 50) {
    cat("[WARN] High variability in library sizes - consider investigating\n")
  }

  n_tissue <- ncol(DGE_tissue$counts)

  logCPM_tissue <- cpm(DGE_tissue, log = TRUE, prior.count = 3)

  keep_tissue <- rowSums(cpm(DGE_tissue) > 1) >= (n_tissue / 2)

  ## ---- SANITY CHECK 4: Gene filtering statistics ----
  cat("\n--- SANITY CHECK 4: Gene filtering ---\n")

  cat("[INFO] Total genes before filtering:", length(keep_tissue), "\n")
  cat("[INFO] Genes passing CPM > 1 in >= 50% samples:", sum(keep_tissue), "\n")
  cat("[INFO] Genes removed:", sum(!keep_tissue), "\n")
  cat("[INFO] Retention rate:", round(sum(keep_tissue)/length(keep_tissue)*100, 1), "%\n")

  ## Typical retention is 15-50% for tissue RNA-seq with CPM > 1 in 50% samples
  ## Lower retention is normal when many genes are tissue-specific or lowly expressed
  retention_rate <- sum(keep_tissue) / length(keep_tissue) * 100
  if (retention_rate < 10) {
    cat("[WARN] Very low gene retention rate (<10%) - check if filter is too stringent\n")
  } else if (retention_rate > 60) {
    cat("[WARN] Very high gene retention rate (>60%) - filter may be too lenient\n")
  } else {
    cat("[PASS] Gene retention rate within expected range (10-60%)\n")
  }

  DGE_tissue_filt <- DGE_tissue[keep_tissue, , keep.lib.sizes = FALSE]

  cat("[INFO] Final gene count:", nrow(DGE_tissue_filt$counts), "\n")

  logCPM_tissue_filt <- cpm(DGE_tissue_filt, log = TRUE, prior.count = 3)

  ## 4.7) Density plots (before/after filter) ---------------

  ## Define density plot function with explicit parameters to avoid global variable issues
  plot_density_comparison <- function(logCPM_before, logCPM_after,
                                      sample_names_before, sample_names_after) {

    n_before <- ncol(logCPM_before)
    n_after  <- ncol(logCPM_after)

    ## Compute density objects for all samples
    dens_before <- lapply(seq_len(n_before), function(i) density(logCPM_before[, i]))
    dens_after  <- lapply(seq_len(n_after),  function(i) density(logCPM_after[, i]))

    ## Compute global xlim from density objects
    xlim_before <- range(sapply(dens_before, function(d) range(d$x)))
    xlim_after  <- range(sapply(dens_after,  function(d) range(d$x)))
    xlim_global <- range(c(xlim_before, xlim_after))

    ## Compute ylim from density objects
    max_y_before <- max(sapply(dens_before, function(d) max(d$y)))
    max_y_after  <- max(sapply(dens_after,  function(d) max(d$y)))

    par(mfrow = c(1, 2), mar = c(5, 5, 5, 2), oma = c(0, 0, 2, 0))
    par(xpd = FALSE)

    cols_before <- viridis(n_before, option = "D", end = 0.95)
    cols_after  <- viridis(n_after, option = "D", end = 0.95)

    ## BEFORE filter
    plot(dens_before[[1]],
         main = "Tissue: log2 CPM (before filter)",
         sub  = paste0("Each line = one sample (n=", n_before, ")"),
         xlab = "log2 CPM",
         lwd  = 2,
         col  = cols_before[1],
         xlim = xlim_global,
         ylim = c(0, max_y_before * 1.15))

    ## FIX: Only loop if more than 1 sample (avoid 2:1 issue)
    if (n_before > 1) {
      for (i in 2:n_before) {
        lines(dens_before[[i]], lwd = 2, col = cols_before[i])
      }
    }

    ## AFTER filter
    plot(dens_after[[1]],
         main = "Tissue: log2 CPM (after filter)",
         sub  = paste0("Each line = one sample (n=", n_after, ")"),
         xlab = "log2 CPM",
         lwd  = 2,
         col  = cols_after[1],
         xlim = xlim_global,
         ylim = c(0, max_y_after * 1.15))

    ## FIX: Only loop if more than 1 sample
    if (n_after > 1) {
      for (i in 2:n_after) {
        lines(dens_after[[i]], lwd = 2, col = cols_after[i])
      }
    }
  }

  density_file_tissue <- paste0("logCPM_density_KAT8KD_tissue_before_after_", run_tag, ".png")
  png(file.path(outdir, "plots", density_file_tissue),
      width = 1600, height = 900, res = 150)
  plot_density_comparison(
    logCPM_before        = logCPM_tissue,
    logCPM_after         = logCPM_tissue_filt,
    sample_names_before  = colnames(DGE_tissue$counts),
    sample_names_after   = colnames(DGE_tissue_filt$counts)
  )
  dev.off()

  cat("\nSaved density plots to",
      file.path(outdir, "plots", density_file_tissue), "\n")

  ## 4.8) PCA -----------------------------------------------

  pca_tissue <- prcomp(t(logCPM_tissue_filt), scale. = TRUE)
  pca_var_tissue <- summary(pca_tissue)$importance[2, 1:2] * 100

  pca_tissue_df <- data.frame(
    Sample   = colnames(DGE_tissue_filt$counts),
    Genotype = DGE_tissue_filt$samples$Genotype,
    DepotSex = DGE_tissue_filt$samples$DepotSex,
    Depot    = DGE_tissue_filt$samples$Depot,
    Sex      = DGE_tissue_filt$samples$Sex,
    PC1      = pca_tissue$x[, 1],
    PC2      = pca_tissue$x[, 2],
    stringsAsFactors = FALSE
  )

  p_tissue <- ggplot(
    pca_tissue_df,
    aes(x = PC1, y = PC2,
        color = DepotSex,
        shape = Genotype,
        label = Sample)
  ) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 60) +
    scale_color_manual(
      values = c(
        "iWAT_F" = "#1b9e77",
        "iWAT_M" = "#d95f02",
        "gWAT_F" = "#7570b3",
        "gWAT_M" = "#e7298a"
      ),
      name = "Depot/Sex"
    ) +
    scale_shape_manual(
      values = c("CTL" = 16, "KAT8KD" = 17),
      name = "Genotype"
    ) +
    labs(
      title = "Tissue PCA (log2 CPM, filtered)",
      x = paste0("PC1 (", round(pca_var_tissue[1], 1), "%)"),
      y = paste0("PC2 (", round(pca_var_tissue[2], 1), "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  pca_file_tissue <- paste0("PCA_KAT8KD_tissue_PC1_PC2_", run_tag, ".png")
  ggsave(
    file.path(outdir, "plots", pca_file_tissue),
    plot   = p_tissue,
    width  = 8,
    height = 6,
    dpi    = 300
  )

  cat("Saved PCA plot to",
      file.path(outdir, "plots", pca_file_tissue), "\n")

  ## 4.9) MDS -----------------------------------------------

  dist_mat <- dist(t(logCPM_tissue_filt))
  mds_coords <- cmdscale(dist_mat, k = 2)
  mds_df <- data.frame(
    Sample   = colnames(DGE_tissue_filt$counts),
    Genotype = DGE_tissue_filt$samples$Genotype,
    DepotSex = DGE_tissue_filt$samples$DepotSex,
    Depot    = DGE_tissue_filt$samples$Depot,
    Sex      = DGE_tissue_filt$samples$Sex,
    MDS1     = mds_coords[, 1],
    MDS2     = mds_coords[, 2],
    stringsAsFactors = FALSE
  )

  mds_plot <- ggplot(
    mds_df,
    aes(x = MDS1, y = MDS2,
        color = DepotSex,
        shape = Genotype,
        label = Sample)
  ) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 60) +
    scale_color_manual(
      values = c(
        "iWAT_F" = "#1b9e77",
        "iWAT_M" = "#d95f02",
        "gWAT_F" = "#7570b3",
        "gWAT_M" = "#e7298a"
      ),
      name = "Depot/Sex"
    ) +
    scale_shape_manual(
      values = c("CTL" = 16, "KAT8KD" = 17),
      name = "Genotype"
    ) +
    labs(
      title = "Tissue MDS (log2 CPM, filtered)",
      x = "MDS1",
      y = "MDS2"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  mds_file_tissue <- paste0("MDS_KAT8KD_tissue_", run_tag, ".png")
  ggsave(
    file.path(outdir, "plots", mds_file_tissue),
    plot   = mds_plot,
    width  = 8,
    height = 6,
    dpi    = 300
  )

  cat("Saved MDS plot to",
      file.path(outdir, "plots", mds_file_tissue), "\n")

  ## 4.10) Design matrix + voom -----------------------------

  design_tissue <- model.matrix(~ 0 + DGE_tissue_filt$samples$GroupFull)
  colnames(design_tissue) <- levels(DGE_tissue_filt$samples$GroupFull)

  ## ---- SANITY CHECK 5: Design matrix validation ----
  cat("\n--- SANITY CHECK 5: Design matrix ---\n")

  cat("[INFO] Design matrix dimensions:", nrow(design_tissue), "samples x",
      ncol(design_tissue), "groups\n")
  cat("[INFO] Group names:", paste(colnames(design_tissue), collapse = ", "), "\n")

  ## Check for full rank
  design_rank <- qr(design_tissue)$rank
  expected_rank <- ncol(design_tissue)
  cat("[INFO] Design matrix rank:", design_rank, "of", expected_rank, "\n")

  if (design_rank < expected_rank) {
    cat("[WARN] Design matrix is not full rank - some coefficients may not be estimable\n")
    cat("[WARN] This could indicate empty groups or perfect collinearity\n")
  } else {
    cat("[PASS] Design matrix is full rank\n")
  }

  ## Check samples per group in design
  samples_per_group <- colSums(design_tissue)
  cat("[INFO] Samples per group in design:\n")
  print(samples_per_group)

  ## Check for empty groups
  empty_groups <- names(samples_per_group)[samples_per_group == 0]
  if (length(empty_groups) > 0) {
    cat("[WARN] Empty groups detected:", paste(empty_groups, collapse = ", "), "\n")
  } else {
    cat("[PASS] All groups have at least one sample\n")
  }

  cat("\nDesign matrix (first rows):\n")
  print(head(design_tissue))

  v_tissue   <- voom(DGE_tissue_filt, design_tissue, plot = FALSE)
  fit_tissue <- lmFit(v_tissue, design_tissue)

  ## ---- SANITY CHECK 6: Voom transformation ----
  cat("\n--- SANITY CHECK 6: Voom transformation ---\n")

  cat("[INFO] Voom E matrix dimensions:", nrow(v_tissue$E), "genes x",
      ncol(v_tissue$E), "samples\n")

  ## Check for extreme values in voom-transformed data
  voom_range <- range(v_tissue$E)
  cat("[INFO] Voom log2-expression range:", round(voom_range[1], 2), "to",
      round(voom_range[2], 2), "\n")

  ## Typical range is roughly -5 to 15
  if (voom_range[1] < -10 || voom_range[2] > 20) {
    cat("[WARN] Voom expression values outside typical range\n")
  } else {
    cat("[PASS] Voom expression values within expected range\n")
  }

  ## Check weights
  cat("[INFO] Voom weights summary:\n")
  print(summary(as.vector(v_tissue$weights)))

  ## 4.11) Define contrasts ---------------------------------

  cont_tissue <- makeContrasts(
    iWAT_F_KD_vs_CTL = iWAT_F_KAT8KD - iWAT_F_CTL,
    iWAT_M_KD_vs_CTL = iWAT_M_KAT8KD - iWAT_M_CTL,
    gWAT_F_KD_vs_CTL = gWAT_F_KAT8KD - gWAT_F_CTL,
    gWAT_M_KD_vs_CTL = gWAT_M_KAT8KD - gWAT_M_CTL,
    levels = design_tissue
  )

  contrast_names <- colnames(cont_tissue)

  ## 4.12) DE thresholds ------------------------------------

  logFC_cut_tissue <- 1
  fdr_cut_tissue   <- 0.05

  ## 4.13) Container for DE summary -------------------------

  deg_summary <- data.frame(
    Contrast          = contrast_names,
    N_DEG_FDR_lt_0_05 = NA_integer_,
    N_Up              = NA_integer_,
    N_Down            = NA_integer_,
    stringsAsFactors  = FALSE
  )

  ## Store per-contrast DE tables
  tt_list <- vector("list", length(contrast_names))
  names(tt_list) <- contrast_names

  ## Store pathway analysis results
  pathway_results <- list()

  ## 4.13.5) Pathway Analysis Function ----------------------

  run_pathway_analysis <- function(gene_list, direction, contrast_name,
                                   run_tag, outdir, fdr_cutoff = 0.05,
                                   msigdb_hallmark, msigdb_kegg, msigdb_gobp) {

    cat("\n--- Pathway analysis:", direction, "genes (", contrast_name, ") ---\n")

    if (length(gene_list) < 5) {
      cat("[WARN] Too few genes (", length(gene_list), ") for pathway analysis\n")
      return(NULL)
    }

    ## Convert gene symbols to Entrez IDs
    gene_entrez <- tryCatch({
      bitr(
        gene_list,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db
      )
    }, error = function(e) {
      cat("[WARN] Gene ID conversion failed:", conditionMessage(e), "\n")
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    })

    if (nrow(gene_entrez) < 3) {
      cat("[WARN] Too few genes mapped to Entrez IDs\n")
      return(NULL)
    }

    cat("[OK] Mapped", nrow(gene_entrez), "out of", length(gene_list), "genes to Entrez IDs\n")

    entrez_ids <- gene_entrez$ENTREZID

    ## Run enrichment analyses
    results <- list()

    ## 1) Hallmark gene sets via enricher (MSigDB)
    tryCatch({
      hallmark_term2gene <- msigdb_hallmark %>% dplyr::select(gs_name, gene_symbol)
      enrich_hallmark <- enricher(
        gene          = gene_list,
        TERM2GENE     = hallmark_term2gene,
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 0.2,
        minGSSize     = 5,
        maxGSSize     = 500
      )
      if (!is.null(enrich_hallmark) && nrow(enrich_hallmark@result) > 0) {
        results$hallmark <- enrich_hallmark
        cat("[OK] Hallmark enrichment:", nrow(enrich_hallmark@result), "pathways\n")
      }
    }, error = function(e) {
      cat("[WARN] Hallmark enrichment failed:", conditionMessage(e), "\n")
    })

    ## 2) KEGG pathways
    tryCatch({
      enrich_kegg <- enrichKEGG(
        gene          = entrez_ids,
        organism      = "mmu",
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 0.2,
        minGSSize     = 5,
        maxGSSize     = 500
      )
      if (!is.null(enrich_kegg) && nrow(enrich_kegg@result) > 0) {
        results$kegg <- enrich_kegg
        cat("[OK] KEGG enrichment:", nrow(enrich_kegg@result), "pathways\n")
      }
    }, error = function(e) {
      cat("[WARN] KEGG enrichment failed:", conditionMessage(e), "\n")
    })

    ## 3) GO Biological Process
    tryCatch({
      enrich_go <- enrichGO(
        gene          = entrez_ids,
        OrgDb         = org.Mm.eg.db,
        ont           = "BP",
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 0.2,
        readable      = TRUE,
        minGSSize     = 5,
        maxGSSize     = 500
      )
      if (!is.null(enrich_go) && nrow(enrich_go@result) > 0) {
        ## Simplify GO terms to reduce redundancy
        enrich_go_simp <- simplify(enrich_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
        results$gobp <- enrich_go_simp
        cat("[OK] GO:BP enrichment:", nrow(enrich_go_simp@result), "terms (simplified)\n")
      }
    }, error = function(e) {
      cat("[WARN] GO:BP enrichment failed:", conditionMessage(e), "\n")
    })

    ## Save results and create plots
    if (length(results) > 0) {
      for (db_name in names(results)) {
        enrich_obj <- results[[db_name]]

        ## Filter by FDR
        sig_results <- enrich_obj@result %>%
          dplyr::filter(p.adjust < fdr_cutoff) %>%
          dplyr::arrange(p.adjust)

        if (nrow(sig_results) == 0) {
          cat("[INFO] No significant pathways in", db_name, "\n")
          next
        }

        ## Save table
        table_file <- paste0("Pathway_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".csv")
        write.csv(
          sig_results,
          file      = file.path(outdir, "tables", table_file),
          row.names = FALSE
        )
        cat("[OK] Saved pathway table:", table_file, "\n")

        ## Create publication-ready bar plot
        if (nrow(sig_results) > 0) {

          ## Take top 15 pathways by p.adjust
          plot_data <- sig_results %>%
            dplyr::slice_head(n = 15) %>%
            dplyr::mutate(
              Description = factor(Description, levels = rev(Description)),
              ## FIX: Safe division - handle potential zero denominator
              GeneRatio_numeric = sapply(strsplit(GeneRatio, "/"), function(x) {
                num <- as.numeric(x[1])
                denom <- as.numeric(x[2])
                if (is.na(denom) || denom == 0) return(0)
                return(num / denom)
              })
            )

          ## Choose color based on direction
          bar_color <- if (direction == "Up") "#D95F02" else "#1B9E77"

          p_pathway <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description, fill = p.adjust)) +
            geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
            scale_fill_viridis(
              option    = "plasma",
              direction = -1,
              name      = "Adjusted\nP-value",
              limits    = c(0, max(plot_data$p.adjust)),
              breaks    = pretty_breaks(n = 5)
            ) +
            labs(
              title = paste0(toupper(db_name), " Pathways (", direction, ")\n", contrast_name),
              x     = "Gene Ratio",
              y     = NULL
            ) +
            theme_classic(base_size = 12) +
            theme(
              plot.title        = element_text(face = "bold", hjust = 0.5, size = 14),
              axis.text.y       = element_text(size = 10, color = "black"),
              axis.text.x       = element_text(size = 10, color = "black", face = "bold"),
              axis.title.x      = element_text(size = 12, face = "bold"),
              legend.title      = element_text(size = 10, face = "bold"),
              legend.text       = element_text(size = 9),
              legend.position   = "right",
              panel.border      = element_rect(color = "black", fill = NA, linewidth = 1)
            )

          plot_file <- paste0("Pathway_barplot_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".png")
          ggsave(
            file.path(outdir, "plots", plot_file),
            plot   = p_pathway,
            width  = 10,
            height = max(6, nrow(plot_data) * 0.3),
            dpi    = 300,
            bg     = "white"
          )
          cat("[OK] Saved pathway bar plot:", plot_file, "\n")
        }
      }
    }

    return(results)
  }

  ## 4.14) Loop over contrasts: DE, volcano, heatmap --------

  for (cn in contrast_names) {

    cat("\n=== Contrast:", cn, "===\n")

    ## ---- DE fit + table -----------------------------------
    fit2 <- contrasts.fit(fit_tissue, cont_tissue[, cn, drop = FALSE])
    fit2 <- eBayes(fit2)

    tt <- topTable(
      fit2,
      coef    = 1,
      n       = Inf,
      sort.by = "P"
    )

    if (!"gene_name" %in% colnames(tt)) {
      tt$gene_name <- rownames(tt)
    }

    ## Store for later reference
    tt_list[[cn]] <- tt

    de_file <- paste0("DE_tissue_", cn, "_", run_tag, ".csv")
    write.csv(
      tt,
      file      = file.path(outdir, "tables", de_file),
      row.names = TRUE
    )

    n_sig <- sum(tt$adj.P.Val < fdr_cut_tissue, na.rm = TRUE)
    cat("DE table written:",
        file.path(outdir, "tables", de_file), "\n")
    cat("Number of genes with FDR < 0.05 (", cn, "):", n_sig, "\n")

    ## Count up/down regulated genes
    n_up <- sum(tt$adj.P.Val < fdr_cut_tissue & tt$logFC > logFC_cut_tissue, na.rm = TRUE)
    n_down <- sum(tt$adj.P.Val < fdr_cut_tissue & tt$logFC < -logFC_cut_tissue, na.rm = TRUE)

    deg_summary$N_DEG_FDR_lt_0_05[deg_summary$Contrast == cn] <- n_sig
    deg_summary$N_Up[deg_summary$Contrast == cn] <- n_up
    deg_summary$N_Down[deg_summary$Contrast == cn] <- n_down

    cat("Number of upregulated genes (", cn, "):", n_up, "\n")
    cat("Number of downregulated genes (", cn, "):", n_down, "\n")

    ## ---- Volcano plot (viridis + FDR + |FC| labels) -------

    tt_plot <- tt %>%
      dplyr::filter(
        is.finite(logFC),
        is.finite(adj.P.Val),
        adj.P.Val > 0
      ) %>%
      mutate(
        negLogFDR = -log10(adj.P.Val)
      )

    logFC_cut <- logFC_cut_tissue
    fdr_cut   <- fdr_cut_tissue

    tt_plot <- tt_plot %>%
      mutate(
        sig_cat = dplyr::case_when(
          adj.P.Val < fdr_cut & logFC >=  logFC_cut  ~ "Up",
          adj.P.Val < fdr_cut & logFC <= -logFC_cut  ~ "Down",
          TRUE                                       ~ "NS"
        ),
        sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
      )

    ## Filter to significant DEs for label selection
    tt_sig <- tt_plot %>%
      dplyr::filter(
        adj.P.Val < fdr_cut,
        abs(logFC) >= logFC_cut,
        sig_cat != "NS"
      )

    top_n_fdr <- 10
    top_n_fc  <- 10

    ## --- Upregulated ---
    up_sig <- tt_sig %>% dplyr::filter(sig_cat == "Up")

    up_top_fdr <- up_sig %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = top_n_fdr)

    up_top_fc <- up_sig %>%
      dplyr::arrange(dplyr::desc(logFC)) %>%
      dplyr::slice_head(n = top_n_fc)

    ## --- Downregulated ---
    down_sig <- tt_sig %>% dplyr::filter(sig_cat == "Down")

    down_top_fdr <- down_sig %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = top_n_fdr)

    down_top_fc <- down_sig %>%
      dplyr::arrange(logFC) %>%
      dplyr::slice_head(n = top_n_fc)

    ## Combine and dedupe
    label_genes <- dplyr::bind_rows(
      up_top_fdr,
      up_top_fc,
      down_top_fdr,
      down_top_fc
    ) %>%
      dplyr::distinct(gene_name, .keep_all = TRUE)

    vol_tissue <- ggplot(
      tt_plot,
      aes(
        x     = logFC,
        y     = negLogFDR,
        color = negLogFDR,
        shape = sig_cat
      )
    ) +
      geom_point(size = 2, alpha = 0.9) +
      scale_color_viridis(
        option    = "rocket",
        direction = 1,
        name      = "-log10(FDR)"
      ) +
      scale_shape_manual(
        values = c(
          "Up"   = 17,
          "Down" = 25,
          "NS"   = 16
        ),
        name = "Direction"
      ) +
      geom_point(
        data = label_genes,
        size = 2.8
      ) +
      geom_text_repel(
        data             = label_genes,
        aes(label        = gene_name),
        color            = "black",
        size             = 3.5,
        box.padding      = unit(0.25, "lines"),
        segment.color    = "gray25",
        segment.size  = 0.25,
        point.padding    = unit(0.2, "lines"),
        max.overlaps     = Inf
      ) +
      geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
      geom_vline(xintercept = c(logFC_cut, -logFC_cut), linetype = "dashed") +
      theme_bw(base_size = 14) +
      theme(
        panel.grid      = element_blank(),
        axis.text       = element_text(size = 10, face = "bold", colour = "black"),
        axis.title      = element_text(size = 16, face = "bold", colour = "black"),
        strip.text      = element_text(size = 10, face = "bold", colour = "black"),
        plot.title      = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom"
      ) +
      ggtitle(paste0("Volcano: ", cn, " (WAT tissue)")) +
      xlab(paste0("log2FC (", cn, ")")) +
      ylab("-log10 FDR")

    volcano_file <- paste0("Volcano_tissue_", cn, "_", run_tag, ".png")
    ggsave(
      file.path(outdir, "plots", volcano_file),
      plot   = vol_tissue,
      width  = 8,
      height = 6,
      dpi    = 300
    )

    cat("Volcano saved:",
        file.path(outdir, "plots", volcano_file), "\n")

    if (is.null(hero_volcano_file)) {
      hero_volcano_file <- volcano_file
    }

    ## ---- Pathway Analysis for Up/Down DEGs ----------------

    ## Extract upregulated genes
    up_genes <- tt %>%
      dplyr::filter(adj.P.Val < fdr_cut_tissue, logFC > logFC_cut_tissue) %>%
      dplyr::pull(gene_name)

    ## Extract downregulated genes
    down_genes <- tt %>%
      dplyr::filter(adj.P.Val < fdr_cut_tissue, logFC < -logFC_cut_tissue) %>%
      dplyr::pull(gene_name)

    ## Run pathway analysis for upregulated genes
    if (length(up_genes) >= 5) {
      pathway_up <- run_pathway_analysis(
        gene_list       = up_genes,
        direction       = "Up",
        contrast_name   = cn,
        run_tag         = run_tag,
        outdir          = outdir,
        fdr_cutoff      = 0.05,
        msigdb_hallmark = msigdb_hallmark,
        msigdb_kegg     = msigdb_kegg,
        msigdb_gobp     = msigdb_gobp
      )
      ## Store results
      pathway_results[[paste0(cn, "_Up")]] <- pathway_up
    } else {
      cat("[WARN] Too few upregulated genes for pathway analysis (", cn, ")\n")
    }

    ## Run pathway analysis for downregulated genes
    if (length(down_genes) >= 5) {
      pathway_down <- run_pathway_analysis(
        gene_list       = down_genes,
        direction       = "Down",
        contrast_name   = cn,
        run_tag         = run_tag,
        outdir          = outdir,
        fdr_cutoff      = 0.05,
        msigdb_hallmark = msigdb_hallmark,
        msigdb_kegg     = msigdb_kegg,
        msigdb_gobp     = msigdb_gobp
      )
      ## Store results
      pathway_results[[paste0(cn, "_Down")]] <- pathway_down
    } else {
      cat("[WARN] Too few downregulated genes for pathway analysis (", cn, ")\n")
    }

    ## ---- Heatmap of strongest DEGs ------------------------

    de_sig <- tt %>%
      dplyr::filter(adj.P.Val < fdr_cut_tissue,
                    abs(logFC) > logFC_cut_tissue)

    if (nrow(de_sig) >= 2) {

      de_sig <- de_sig[order(de_sig$adj.P.Val), , drop = FALSE]
      if (nrow(de_sig) > 300) {
        de_sig <- de_sig[1:300, , drop = FALSE]
      }
      gene_set <- rownames(de_sig)

      depot <- NA_character_
      sex   <- NA_character_
      if (grepl("^iWAT_F", cn)) {
        depot <- "iWAT"; sex <- "F"
      } else if (grepl("^iWAT_M", cn)) {
        depot <- "iWAT"; sex <- "M"
      } else if (grepl("^gWAT_F", cn)) {
        depot <- "gWAT"; sex <- "F"
      } else if (grepl("^gWAT_M", cn)) {
        depot <- "gWAT"; sex <- "M"
      }

      if (!is.na(depot) && !is.na(sex)) {

        idx_samples <- with(DGE_tissue_filt$samples,
                            Depot == depot & Sex == sex)
        samples_heat <- rownames(DGE_tissue_filt$samples)[idx_samples]

        mat <- logCPM_tissue_filt[gene_set, samples_heat, drop = FALSE]

        palette_hm <- colorRampPalette(c("blue", "white", "red"))(100)

        annot_col <- data.frame(
          Genotype = DGE_tissue_filt$samples$Genotype[idx_samples],
          Depot    = DGE_tissue_filt$samples$Depot[idx_samples],
          Sex      = DGE_tissue_filt$samples$Sex[idx_samples]
        )
        rownames(annot_col) <- samples_heat

        annot_colors <- list(
          Genotype = c(CTL = "#1B9E77", KAT8KD = "#D95F02"),
          Depot    = c(iWAT = "#377eb8", gWAT = "#e41a1c"),
          Sex      = c(F = "#984ea3", M = "#ff7f00")
        )

        heat_file <- paste0("Heatmap_tissue_", cn, "_", run_tag, ".png")
        png(file.path(outdir, "plots", heat_file),
            width = 3000, height = 3500, res = 300)
        pheatmap(
          mat,
          show_rownames     = TRUE,
          show_colnames     = TRUE,
          scale             = "row",
          cluster_cols      = FALSE,
          cluster_rows      = TRUE,
          color             = palette_hm,
          border_color      = NA,
          fontsize_row      = 6,
          fontsize_col      = 10,
          annotation_col    = annot_col,
          annotation_colors = annot_colors,
          annotation_legend = TRUE,
          main              = paste0("Top DEGs ", cn,
                                     " (FDR<0.05, |logFC|>", logFC_cut_tissue, ")")
        )
        dev.off()
        cat("Heatmap saved:",
            file.path(outdir, "plots", heat_file), "\n")

      } else {
        cat("Could not map contrast", cn, "to depot/sex for heatmap.\n")
      }
    } else {
      cat("Not enough DE genes for heatmap (", cn, ").\n")
    }

  } ## end contrast loop

  ## 4.15) DEG summary table -------------------------------

  deg_summary_file <- paste0("DEG_summary_tissue_", run_tag, ".csv")
  write.csv(
    deg_summary,
    file = file.path(outdir, "tables", deg_summary_file),
    row.names = FALSE
  )
  cat("\nDEG summary written:",
      file.path(outdir, "tables", deg_summary_file), "\n")

  ## ---- SANITY CHECK 7: DE results validation ----
  cat("\n--- SANITY CHECK 7: DE results validation ---\n")

  cat("[INFO] DEG Summary across all contrasts:\n")
  print(deg_summary)

  ## Check for contrasts with zero DEGs
  zero_deg_contrasts <- deg_summary$Contrast[deg_summary$N_DEG_FDR_lt_0_05 == 0]
  if (length(zero_deg_contrasts) > 0) {
    cat("[WARN] Contrasts with zero DEGs at FDR < 0.05:",
        paste(zero_deg_contrasts, collapse = ", "), "\n")
    cat("[INFO] This may indicate weak effects or need for adjusted thresholds\n")
  } else {
    cat("[PASS] All contrasts have at least one significant DEG\n")
  }

  ## Check for extremely high DEG counts (may indicate issues)
  high_deg_threshold <- nrow(DGE_tissue_filt$counts) * 0.3  ## 30% of genes
  high_deg_contrasts <- deg_summary$Contrast[deg_summary$N_DEG_FDR_lt_0_05 > high_deg_threshold]
  if (length(high_deg_contrasts) > 0) {
    cat("[WARN] Contrasts with unusually high DEG count (>30% of genes):",
        paste(high_deg_contrasts, collapse = ", "), "\n")
    cat("[INFO] Consider checking for batch effects or sample swaps\n")
  }

  ## Check for asymmetry (much more up than down or vice versa)
  for (i in seq_len(nrow(deg_summary))) {
    cn <- deg_summary$Contrast[i]
    n_up <- deg_summary$N_Up[i]
    n_down <- deg_summary$N_Down[i]
    total <- n_up + n_down

    if (total > 10) {  ## Only check if there are enough DEGs
      ratio <- max(n_up, n_down) / max(min(n_up, n_down), 1)
      if (ratio > 5) {
        cat("[INFO]", cn, ": Strong asymmetry in DEGs (Up:", n_up, ", Down:", n_down, ")\n")
      }
    }
  }

  cat("[PASS] DE results validation complete\n")

  ## 4.16) save_core bookkeeping ----------------------------

  if (!is.null(hero_volcano_file)) {
    hero_path <- file.path(outdir, "plots", hero_volcano_file)
    if (file.exists(hero_path)) {
      save_run_file(hero_path, tag = "hero_volcano")
    }
  }

  try(dedupe(outdir), silent = TRUE)

  ## 4.17) Save session information for reproducibility -----

  session_file <- paste0("sessionInfo_", run_tag, ".txt")
  sink(file.path(outdir, "logs", session_file))
  cat("=======================================================\n")
  cat("SESSION INFORMATION FOR REPRODUCIBILITY\n")
  cat("=======================================================\n")
  cat("\nScript:", run_info$script_name, "\n")
  cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Run tag:", run_tag, "\n")
  cat("\nThresholds used:\n")
  cat("  - logFC cutoff:", logFC_cut_tissue, "\n")
  cat("  - FDR cutoff:", fdr_cut_tissue, "\n")
  cat("\nOutlier samples removed:", paste(outlier_samples, collapse = ", "), "\n")
  cat("\n")
  print(sessionInfo())
  sink()

  cat("\n[OK] Session info saved:", file.path(outdir, "logs", session_file), "\n")

  ## 4.18) Final Sanity Check Summary -----------------------

  cat("\n=======================================================\n")
  cat("SANITY CHECK SUMMARY\n")
  cat("=======================================================\n")
  cat("1. Counts file structure:     Validated\n")
  cat("2. Sample annotation:         ", length(tissue_samples), "samples after outlier removal\n")
  cat("3. Library sizes:             ", round(min(lib_sizes)/1e6, 1), "-",
      round(max(lib_sizes)/1e6, 1), "M reads\n")
  cat("4. Gene filtering:            ", nrow(DGE_tissue_filt$counts), "genes retained (",
      round(retention_rate, 1), "%)\n")
  cat("5. Design matrix:             ", design_rank, "/", expected_rank, "rank (full rank)\n")
  cat("6. Voom transformation:       Expression range",
      round(voom_range[1], 1), "to", round(voom_range[2], 1), "\n")
  cat("7. DE analysis:               ", sum(deg_summary$N_DEG_FDR_lt_0_05), "total DEGs across",
      nrow(deg_summary), "contrasts\n")
  cat("=======================================================\n")

  ## Save sanity check summary to file
  sanity_file <- paste0("sanity_check_summary_", run_tag, ".txt")
  sanity_lines <- c(
    "=======================================================",
    "SANITY CHECK SUMMARY",
    "=======================================================",
    paste0("Run tag: ", run_tag),
    paste0("Generated: ", Sys.time()),
    "",
    "--- Input Data ---",
    paste0("Counts file: ", counts_file),
    paste0("Total genes in input: ", nrow(counts_raw)),
    paste0("Total samples in input: ", length(sample_cols)),
    "",
    "--- Sample Processing ---",
    paste0("Expected samples: ", length(expected_samples)),
    paste0("Outliers removed: ", paste(outlier_samples, collapse = ", ")),
    paste0("Final sample count: ", length(tissue_samples)),
    "",
    "--- Quality Metrics ---",
    paste0("Library size range: ", round(min(lib_sizes)/1e6, 2), " - ",
           round(max(lib_sizes)/1e6, 2), " million reads"),
    paste0("Library size CV: ", round(lib_cv, 1), "%"),
    paste0("Normalization factor range: ", round(min(norm_factors), 3), " - ",
           round(max(norm_factors), 3)),
    "",
    "--- Gene Filtering ---",
    paste0("Genes before filtering: ", length(keep_tissue)),
    paste0("Genes after filtering: ", nrow(DGE_tissue_filt$counts)),
    paste0("Retention rate: ", round(retention_rate, 1), "%"),
    paste0("Filter criteria: CPM > 1 in >= ", round(n_tissue/2), " samples"),
    "",
    "--- Design Matrix ---",
    paste0("Samples: ", nrow(design_tissue)),
    paste0("Groups: ", ncol(design_tissue)),
    paste0("Rank: ", design_rank, " of ", expected_rank),
    "",
    "--- Differential Expression ---",
    paste0("Contrasts tested: ", length(contrast_names)),
    paste0("logFC cutoff: ", logFC_cut_tissue),
    paste0("FDR cutoff: ", fdr_cut_tissue),
    "",
    "--- DEG Counts per Contrast ---"
  )

  for (i in seq_len(nrow(deg_summary))) {
    sanity_lines <- c(sanity_lines,
      paste0("  ", deg_summary$Contrast[i], ": ",
             deg_summary$N_DEG_FDR_lt_0_05[i], " total (",
             deg_summary$N_Up[i], " up, ",
             deg_summary$N_Down[i], " down)"))
  }

  sanity_lines <- c(sanity_lines,
    "",
    "--- Output Files Generated ---",
    paste0("Plots: ", length(list.files(file.path(outdir, "plots")))),
    paste0("Tables: ", length(list.files(file.path(outdir, "tables")))),
    paste0("Logs: ", length(list.files(file.path(outdir, "logs")))),
    "",
    "======================================================="
  )

  writeLines(sanity_lines, con = file.path(outdir, "logs", sanity_file))
  cat("[OK] Sanity check summary saved:", file.path(outdir, "logs", sanity_file), "\n")

  cat("\n=======================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("=======================================================\n")
  cat("Output directory:", normalizePath(outdir), "\n")
  cat("=======================================================\n\n")

}, error = function(e) {

  cat("\nERROR:\n")
  cat(conditionMessage(e), "\n")

  ## Write error log into run folder
  err_file <- paste0("ERROR_", run_tag, ".txt")
  writeLines(
    c("Script error:", conditionMessage(e)),
    con = file.path(outdir, "logs", err_file)
  )

  stop(e)

})
