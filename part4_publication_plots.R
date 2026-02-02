#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: PUBLICATION QUALITY PLOTS
## =========================================================
## 1. Dot plots for ORA pathway enrichment (GO:BP + KEGG)
## 2. Grouped heatmaps with two columns (CTL vs KAT8KD)
##
## Run AFTER parts 1-3 have completed
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "viridis")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "ComplexHeatmap",
  "circlize"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

## 2) Find most recent run directory -----------------------
cat("\n=== FINDING ANALYSIS RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run parts 1-3 first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No run directories found in savepoints/")
}

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")

## Use centralized parameters from parameters.R
ora_dotplot_params <- list(
  fdr_cutoff  = dotplot_params$fdr_cutoff,
  min_count   = dotplot_params$min_count,
  top_n       = list("GO:BP" = dotplot_params$top_n_gobp, "KEGG" = dotplot_params$top_n_kegg),
  default_top_n = dotplot_params$top_n_gobp,
  order_by    = "adj. p-value",
  tie_breaker = "gene ratio",
  gobp_simplify_cutoff = dotplot_params$gobp_simplify_cutoff
)

## heatmap_params is already defined in parameters.R


## ============================================================
## NOTE: DOT PLOTS AND BAR PLOTS ARE IN part4b_visualizations.R
## ============================================================
## For ORA bar plots, fGSEA bar plots, and dot plots,
## use part4b_visualizations.R instead of this script.
## This script now only handles HEATMAPS.
## ============================================================

cat("\n=== SKIPPING DOT PLOTS (use part4b_visualizations.R) ===\n")


## ============================================================
## SECTION 2: GROUPED HEATMAPS WITH LOG2FC COLORING
## ============================================================
## Heatmap style (journal-ready):
## - Two columns: CTL and KAT8KD
## - Rows grouped by custom gene categories
## - Values: DESeq2 log2FoldChange for coloring
## - Color scale: centered at 0, clipped at ±2
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## ============================================================
## HEATMAP SETTINGS
## ============================================================
heatmap_lfc_clip <- heatmap_params$lfc_clip       ## Clip log2FC at ±2 (use 1.5 for subtle effects)
## ============================================================

## ============================================================
## >>> PATHWAY-BASED GENE CATEGORIES FROM FGSEA RESULTS <<<
## ============================================================
## These genes are automatically extracted from FGSEA leading edge
## to show genes that actually drive pathway enrichment

## Function to extract leading edge genes from FGSEA CSV
extract_leading_edge <- function(fgsea_file, pathway_patterns, max_genes_per_pathway = 15) {
  if (!file.exists(fgsea_file)) {
    cat("[WARN] FGSEA file not found: ", fgsea_file, "\n")
    return(character(0))
  }

  fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)

  genes <- character(0)
  for (pattern in pathway_patterns) {
    ## Find pathway by partial match (case-insensitive)
    matched_rows <- grep(pattern, fgsea_data$pathway_name, ignore.case = TRUE)
    if (length(matched_rows) > 0) {
      ## Take the first (most significant) match
      row <- matched_rows[1]
      le_symbols <- fgsea_data$leadingEdge_symbols[row]
      if (!is.na(le_symbols) && le_symbols != "") {
        gene_vec <- unlist(strsplit(le_symbols, ","))
        ## Take top genes
        genes <- c(genes, head(gene_vec, max_genes_per_pathway))
      }
    }
  }
  unique(genes)
}

## Find FGSEA GO:BP files
fgsea_gobp_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)
## Also check project root (where user uploaded files)
fgsea_gobp_root <- list.files(getwd(), pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)
fgsea_gobp_files <- c(fgsea_gobp_files, fgsea_gobp_root)

## Initialize gene categories - will be populated from FGSEA or use defaults
gene_categories <- list()

## Try to extract genes from iWAT FGSEA (has strong mitochondrial signal)
iwat_fgsea <- fgsea_gobp_files[grep("iWAT", fgsea_gobp_files, ignore.case = TRUE)]
if (length(iwat_fgsea) > 0) {
  iwat_fgsea <- iwat_fgsea[1]
  cat("[INFO] Extracting genes from iWAT FGSEA: ", basename(iwat_fgsea), "\n")

  ## Mitochondrial/Energy genes (DOWN in iWAT)
  mito_genes <- extract_leading_edge(iwat_fgsea, c(
    "mitochondrial ATP synthesis",
    "mitochondrial translation",
    "respiratory chain complex",
    "fatty acid beta-oxidation"
  ), max_genes_per_pathway = 10)

  if (length(mito_genes) > 0) {
    gene_categories[["Mitochondrial /\nEnergy"]] <- head(mito_genes, 20)
    cat("[OK] Found ", length(mito_genes), " mitochondrial genes\n")
  }
}

## Try to extract genes from gWAT FGSEA (has strong immune signal)
gwat_fgsea <- fgsea_gobp_files[grep("gWAT", fgsea_gobp_files, ignore.case = TRUE)]
if (length(gwat_fgsea) > 0) {
  gwat_fgsea <- gwat_fgsea[1]
  cat("[INFO] Extracting genes from gWAT FGSEA: ", basename(gwat_fgsea), "\n")

  ## Immune/Inflammatory genes (UP in both tissues)
  immune_genes <- extract_leading_edge(gwat_fgsea, c(
    "leukocyte degranulation",
    "myeloid cell activation",
    "positive regulation of cytokine"
  ), max_genes_per_pathway = 10)

  if (length(immune_genes) > 0) {
    gene_categories[["Immune /\nInflammatory"]] <- head(immune_genes, 20)
    cat("[OK] Found ", length(immune_genes), " immune genes\n")
  }

  ## Lipid metabolism genes
  lipid_genes <- extract_leading_edge(gwat_fgsea, c(
    "lipid localization",
    "lipid catabolic",
    "fatty acid"
  ), max_genes_per_pathway = 10)

  if (length(lipid_genes) > 0) {
    gene_categories[["Lipid\nMetabolism"]] <- head(lipid_genes, 15)
    cat("[OK] Found ", length(lipid_genes), " lipid metabolism genes\n")
  }

  ## Tissue remodeling / ECM genes
  ecm_genes <- extract_leading_edge(gwat_fgsea, c(
    "tissue remodeling",
    "extracellular matrix",
    "cell adhesion"
  ), max_genes_per_pathway = 10)

  if (length(ecm_genes) > 0) {
    gene_categories[["ECM /\nRemodeling"]] <- head(ecm_genes, 15)
    cat("[OK] Found ", length(ecm_genes), " ECM/remodeling genes\n")
  }
}

## Fallback: if no FGSEA files found, use default categories
if (length(gene_categories) == 0) {
  cat("[INFO] No FGSEA files found, using default gene categories\n")
  gene_categories <- list(

    "Inflammatory\nResponse" = c(
      "Il1b", "Il6", "Tnf", "Il1a", "Ccl2", "Ccl7", "Cxcl12"
    ),

    "ECM /\nCollagens" = c(
      "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
      "Col5a1", "Col6a1", "Col6a2", "Col6a3"
    ),

    "Matrix\nRemodeling" = c(
      "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
      "Timp1", "Timp2", "Timp3"
    ),

    "Adipocyte\nMarkers" = c(
      "Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1"
    )
  )
}

cat("[INFO] Gene categories: ", paste(names(gene_categories), collapse = ", "), "\n")
cat("[INFO] Total genes: ", sum(sapply(gene_categories, length)), "\n")

## ============================================================
## >>> END OF CUSTOM GENE SECTION <<<
## ============================================================


## Function to create grouped heatmap with Z-scored VST expression
## UPDATED: Always uses Z-scored VST (shows actual expression patterns in BOTH groups)
## This avoids the "binary" appearance of CTL=0 reference approach
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    lfc_clip = heatmap_lfc_clip,
                                    output_dir = plots_dir,
                                    vst_data = NULL) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  if (!file.exists(de_file)) {
    cat("[WARN] DE file not found: ", de_file, "\n")
    return(NULL)
  }

  ## Check if VST data is available - REQUIRED for proper grouped heatmap
  if (is.null(vst_data) || is.null(vst_data$vst_mat) || is.null(vst_data$sample_info)) {
    cat("[WARN] VST data required for grouped heatmap (avoids binary CTL=0 display)\n")
    cat("[INFO] Re-run part1_main_analysis.R to generate VST data\n")
    return(NULL)
  }

  tissue <- sub("_.*", "", contrast_name)
  sample_info <- vst_data$sample_info

  ## Find tissue samples
  if ("Depot" %in% colnames(sample_info)) {
    tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
  } else if ("Tissue" %in% colnames(sample_info)) {
    tissue_samples <- rownames(sample_info)[sample_info$Tissue == tissue]
  } else {
    tissue_samples <- colnames(vst_data$vst_mat)[grepl(tissue, colnames(vst_data$vst_mat), ignore.case = TRUE)]
  }

  if (length(tissue_samples) == 0 || !"Genotype" %in% colnames(sample_info)) {
    cat("[WARN] Cannot find tissue samples or Genotype column\n")
    return(NULL)
  }

  vst_tissue <- vst_data$vst_mat[, colnames(vst_data$vst_mat) %in% tissue_samples, drop = FALSE]
  sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

  if (!all(c("CTL", "KAT8KD") %in% sample_info_tissue$Genotype)) {
    cat("[WARN] Both CTL and KAT8KD genotypes required\n")
    return(NULL)
  }

  ## Find genes in data
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(vst_tissue))

  cat("[INFO] Category genes found: ", length(genes_in_data), "/",
      length(all_category_genes), "\n", sep = "")

  if (length(genes_in_data) < 3) {
    cat("[WARN] Too few genes found for heatmap\n")
    return(NULL)
  }

  ## Map genes to categories
  gene_to_category <- character(length(genes_in_data))
  names(gene_to_category) <- genes_in_data

  for (cat_name in names(gene_categories)) {
    cat_genes <- intersect(gene_categories[[cat_name]], genes_in_data)
    for (g in cat_genes) {
      if (gene_to_category[g] == "") {
        gene_to_category[g] <- cat_name
      }
    }
  }

  gene_to_category <- gene_to_category[gene_to_category != ""]
  genes_to_plot <- names(gene_to_category)

  if (length(genes_to_plot) < 3) {
    cat("[WARN] Too few categorized genes for heatmap\n")
    return(NULL)
  }

  ## Calculate group means (VST expression)
  genotype_levels <- c("CTL", "KAT8KD")
  group_means <- sapply(genotype_levels, function(geno) {
    rowMeans(vst_tissue[, sample_info_tissue$Genotype == geno, drop = FALSE], na.rm = TRUE)
  })
  mat <- group_means[genes_to_plot, , drop = FALSE]

  ## Order by category
  gene_order <- c()
  for (cat_name in names(gene_categories)) {
    cat_genes <- names(gene_to_category[gene_to_category == cat_name])
    gene_order <- c(gene_order, cat_genes)
  }
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Z-score normalize ACROSS BOTH GROUPS (shows relative expression)
  mat <- t(scale(t(mat)))
  mat <- mat[!rowSums(is.na(mat)), , drop = FALSE]

  ## Clip at specified value
  z_clip <- lfc_clip
  mat[mat > z_clip] <- z_clip
  mat[mat < -z_clip] <- -z_clip

  ## TONED DOWN COLORS: Use muted diverging palette from parameters.R
  ## Blue-grey-red: publication-ready, less saturated
  col_fun <- colorRamp2(
    c(-z_clip, -z_clip/2, 0, z_clip/2, z_clip),
    heatmap_params$muted_colors
  )

  param_caption <- paste0(
    "Z-scored mean VST by genotype | Clipped at \u00b1", z_clip
  )
  legend_title <- "Z-score\n(VST mean)"
  legend_at <- c(-z_clip, 0, z_clip)
  legend_labels <- c(paste0("\u2264", -z_clip), "0", paste0("\u2265", z_clip))
  heatmap_name <- "Z-score"

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = heatmap_name,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_split = factor(gene_to_category, levels = names(gene_categories)),
    row_gap = unit(1.5, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title = contrast_name,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    width = unit(2.5, "cm"),
    heatmap_legend_param = list(
      title = legend_title,
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(2.5, "cm"),
      at = legend_at,
      labels = legend_labels
    )
  )

  ## Save
  plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = 600,
      height = max(400, n_genes * 18 + n_categories * 30 + 50),
      res = 150)
  draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
  ## Add parameter caption at bottom
  grid.text(
    param_caption,
    x = 0.5, y = unit(3, "mm"),
    just = "center",
    gp = gpar(fontsize = 7, col = "grey40")
  )
  dev.off()

  cat("[OK] Saved heatmap: ", plot_file, "\n", sep = "")

  return(ht)
}

## ============================================================
## Function to create INDIVIDUAL SAMPLE heatmap with z-scored VST
## This shows expression gradation across all samples
## ============================================================
create_individual_heatmap <- function(vst_data, de_file, contrast_name, gene_categories,
                                       output_dir = plots_dir) {

  cat("\n--- Creating individual sample heatmap: ", contrast_name, " ---\n", sep = "")

  if (is.null(vst_data)) {
    cat("[WARN] VST data not available\n")
    return(NULL)
  }

  vst_mat <- vst_data$vst_mat
  sample_info <- vst_data$sample_info

  ## Identify samples for this contrast (extract tissue from contrast name)
  tissue <- sub("_.*", "", contrast_name)  ## e.g., "iWAT" from "iWAT_KAT8KD_vs_CTL"

  ## Filter samples for this tissue
  if ("Depot" %in% colnames(sample_info)) {
    tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
  } else if ("Tissue" %in% colnames(sample_info)) {
    tissue_samples <- rownames(sample_info)[sample_info$Tissue == tissue]
  } else {
    ## Fallback: try to match sample names
    tissue_samples <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
  }

  if (length(tissue_samples) == 0) {
    cat("[WARN] No samples found for tissue: ", tissue, "\n")
    return(NULL)
  }

  ## Filter VST matrix to tissue samples
  vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples, drop = FALSE]
  sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

  ## Find genes in categories
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(vst_tissue))

  cat("[INFO] Category genes found in VST: ", length(genes_in_data), "/",
      length(all_category_genes), "\n", sep = "")

  if (length(genes_in_data) < 3) {
    cat("[WARN] Too few genes found for heatmap\n")
    return(NULL)
  }

  ## Map genes to categories
  gene_to_category <- character(length(genes_in_data))
  names(gene_to_category) <- genes_in_data

  for (cat_name in names(gene_categories)) {
    cat_genes <- intersect(gene_categories[[cat_name]], genes_in_data)
    for (g in cat_genes) {
      if (gene_to_category[g] == "") {
        gene_to_category[g] <- cat_name
      }
    }
  }

  gene_to_category <- gene_to_category[gene_to_category != ""]
  genes_to_plot <- names(gene_to_category)

  if (length(genes_to_plot) < 3) {
    cat("[WARN] Too few categorized genes for heatmap\n")
    return(NULL)
  }

  ## Extract VST values for selected genes
  mat <- vst_tissue[genes_to_plot, , drop = FALSE]

  ## Z-score normalize per gene (row)
  mat_scaled <- t(scale(t(mat)))

  ## Order genes by category
  gene_order <- c()
  for (cat_name in names(gene_categories)) {
    cat_genes <- names(gene_to_category[gene_to_category == cat_name])
    gene_order <- c(gene_order, cat_genes)
  }
  mat_scaled <- mat_scaled[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Order samples by genotype (CTL first, then KAT8KD)
  if ("Genotype" %in% colnames(sample_info_tissue)) {
    sample_order <- order(sample_info_tissue$Genotype)
    mat_scaled <- mat_scaled[, sample_order, drop = FALSE]
    sample_info_tissue <- sample_info_tissue[sample_order, , drop = FALSE]
  }

  ## Clip extreme z-scores for visualization
  z_clip <- 2.5
  mat_scaled[mat_scaled > z_clip] <- z_clip
  mat_scaled[mat_scaled < -z_clip] <- -z_clip

  ## TONED DOWN COLORS: Use muted diverging palette from parameters.R
  ## Blue-grey-red: publication-ready, less saturated
  col_fun <- colorRamp2(
    c(-z_clip, -z_clip/2, 0, z_clip/2, z_clip),
    heatmap_params$muted_colors
  )

  ## Column annotation for genotype
  if ("Genotype" %in% colnames(sample_info_tissue)) {
    col_anno <- HeatmapAnnotation(
      Genotype = sample_info_tissue$Genotype,
      col = list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02")),
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
      simple_anno_size = unit(4, "mm")
    )
  } else {
    col_anno <- NULL
  }

  ## Parameter caption
  param_caption <- paste0(
    "Z-scored VST expression | Individual samples | Clipped at \u00b1", z_clip
  )

  ## Create heatmap
  ht <- Heatmap(
    mat_scaled,
    col = col_fun,
    name = "Z-score",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    top_annotation = col_anno,
    row_split = factor(gene_to_category, levels = names(gene_categories)),
    row_gap = unit(1.5, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    column_title = paste0(contrast_name, " (Individual Samples)"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    heatmap_legend_param = list(
      title = "Z-score\n(VST)",
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(2.5, "cm"),
      at = c(-z_clip, 0, z_clip),
      labels = c(paste0("\u2264", -z_clip), "0", paste0("\u2265", z_clip))
    )
  )

  ## Save
  plot_file <- paste0("Heatmap_individual_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_samples <- ncol(mat_scaled)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = max(600, n_samples * 50 + 200),
      height = max(400, n_genes * 18 + n_categories * 30 + 100),
      res = 150)
  draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
  ## Add parameter caption at bottom
  grid.text(
    param_caption,
    x = 0.5, y = unit(3, "mm"),
    just = "center",
    gp = gpar(fontsize = 7, col = "grey40")
  )
  dev.off()

  cat("[OK] Saved individual heatmap: ", plot_file, "\n", sep = "")

  return(ht)
}


## ============================================================
## Create BOTH heatmap types for each contrast
## ============================================================

## Load VST data if available
vst_files <- list.files(tables_dir, pattern = "^VST_matrix_heatmap_.*\\.rds$", full.names = TRUE)
vst_data <- NULL
if (length(vst_files) > 0) {
  vst_file <- vst_files[order(file.info(vst_files)$mtime, decreasing = TRUE)][1]
  vst_data <- readRDS(vst_file)
  cat("[OK] Loaded VST matrix from: ", basename(vst_file), "\n", sep = "")
} else {
  cat("[INFO] VST matrix not found - individual sample heatmaps will be skipped\n")
  cat("[INFO] Re-run part1_main_analysis.R to generate VST data\n")
}

## Create heatmaps for each contrast
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

for (de_file in de_files) {
  filename <- basename(de_file)
  contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
  if (!is.na(contrast)) {
    ## Create grouped heatmap (log2FC: CTL=0, KAT8KD=log2FC)
    create_grouped_heatmap(de_file, contrast, gene_categories, vst_data = vst_data)

    ## Create individual sample heatmap (z-scored VST)
    if (!is.null(vst_data)) {
      create_individual_heatmap(vst_data, de_file, contrast, gene_categories)
    }
  }
}


cat("\n======================================\n")
cat("=== PART 4 COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - GO:BP enrichment dot plots\n")
cat("  - KEGG enrichment dot plots\n")
cat("  - Grouped gene heatmaps (CTL vs KAT8KD log2FC)\n")
if (!is.null(vst_data)) {
  cat("  - Individual sample heatmaps (z-scored VST)\n")
}
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
