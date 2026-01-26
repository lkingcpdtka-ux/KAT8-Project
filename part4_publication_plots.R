#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: PUBLICATION QUALITY PLOTS
## =========================================================
## This script creates publication-quality visualizations:
## 1. Dot plots for ORA pathway enrichment (GO:BP + KEGG)
## 2. Grouped heatmaps with custom gene categories
##
## Run AFTER parts 1-3 have completed
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr")
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
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

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


## ============================================================
## SECTION 1: DOT PLOTS FOR ORA PATHWAY ENRICHMENT
## ============================================================
## - Processes both GO:BP and KEGG results
## - Ordered by adjusted p-value (most significant at bottom)
## - Top 15 pathways
## - Small filled dots
## ============================================================

cat("\n=== CREATING ORA ENRICHMENT DOT PLOTS ===\n")

create_enrichment_dotplot <- function(ora_file, contrast_name, direction,
                                       database = "GO:BP",
                                       top_n = 15, output_dir = plots_dir) {

  if (!file.exists(ora_file)) {
    cat("[WARN] File not found: ", ora_file, "\n")
    return(NULL)
  }

  ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)

  if (nrow(ora_data) == 0) {
    cat("[WARN] No pathways in: ", basename(ora_file), "\n")
    return(NULL)
  }

  ## Parse GeneRatio to numeric
  if ("GeneRatio" %in% colnames(ora_data) && is.character(ora_data$GeneRatio)) {
    ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      num <- as.numeric(x[1])
      denom <- as.numeric(x[2])
      if (is.na(denom) || denom == 0) return(0)
      num / denom
    })
  }

  ## Get gene count
  if (!"Count" %in% colnames(ora_data)) {
    ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1])
    })
  }

  ## Select top pathways by p-value, order by p-value for display
  plot_data <- ora_data %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(desc(p.adjust)) %>%  ## Reverse so most significant at bottom
    dplyr::mutate(
      Description = str_wrap(Description, width = 42),
      Description = factor(Description, levels = Description)
    )

  if (nrow(plot_data) == 0) {
    cat("[WARN] No significant pathways for dotplot\n")
    return(NULL)
  }

  ## Direction-specific colors (orange=up, teal=down)
  if (direction == "Up") {
    color_low <- "#FDBE85"
    color_high <- "#D94701"
  } else {
    color_low <- "#9ECAE1"
    color_high <- "#08519C"
  }

  ## Create dot plot
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(
      low = color_high,
      high = color_low,
      name = "Adj.\nP-value"
    ) +
    scale_size_continuous(
      name = "Gene\nCount",
      range = c(1.5, 5),
      breaks = pretty(plot_data$Count, n = 3)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = paste0(database, ": ", contrast_name, " (", direction, ")"),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.title.x = element_text(size = 9, face = "bold"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.35, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    )

  ## Save
  db_short <- gsub(":", "", database)
  plot_file <- paste0("Dotplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 7,
    height = max(3.5, nrow(plot_data) * 0.22 + 1),
    dpi = 300,
    bg = "white"
  )
  cat("[OK] Saved dot plot: ", plot_file, "\n", sep = "")

  return(p)
}

## Process GO:BP results
cat("\n--- Processing GO:BP results ---\n")
gobp_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)

for (ora_file in gobp_files) {
  filename <- basename(ora_file)
  parts <- str_match(filename, "ORA_gobp_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
  if (!is.na(parts[1])) {
    contrast <- parts[2]
    direction <- parts[3]
    create_enrichment_dotplot(ora_file, contrast, direction, database = "GO:BP")
  }
}

## Process KEGG results
cat("\n--- Processing KEGG results ---\n")
kegg_files <- list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE)

for (ora_file in kegg_files) {
  filename <- basename(ora_file)
  parts <- str_match(filename, "ORA_kegg_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
  if (!is.na(parts[1])) {
    contrast <- parts[2]
    direction <- parts[3]
    create_enrichment_dotplot(ora_file, contrast, direction, database = "KEGG")
  }
}


## ============================================================
## SECTION 2: GROUPED HEATMAPS WITH CUSTOM GENE CATEGORIES
## ============================================================
## - If counts.txt exists: loads expression data, calculates z-scores
## - If not: shows log2 fold change (single column, correct representation)
## - Color scale: Teal (down) - White (center) - Orange (up)
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## ============================================================
## >>> ADD YOUR CUSTOM GENES HERE <<<
## ============================================================

gene_categories <- list(

  "Inflammatory\nResponse" = c(
    "Il1b", "Il6", "Tnf", "Il1a"
    ## ADD MORE HERE
  ),

  "Chemotaxis" = c(
    "Ccl2", "Ccl7", "Cxcl12"
    ## ADD MORE HERE
  ),

  "ECM /\nCollagens" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
    "Col5a1", "Col5a3", "Col6a1", "Col6a2", "Col6a3",
    "Col6a6", "Col15a1"
    ## ADD MORE HERE
  ),

  "Matrix\nRemodeling" = c(
    "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
    "Timp1", "Timp2", "Timp3", "Timp4"
    ## ADD MORE HERE
  ),

  "Adipocyte\nMarkers" = c(
    "Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1"
    ## ADD MORE HERE
  )

)

## ============================================================
## >>> END OF CUSTOM GENE SECTION <<<
## ============================================================


## Check if counts.txt exists for expression-based heatmaps
counts_file <- file.path(getwd(), "counts.txt")
has_counts <- file.exists(counts_file)

if (has_counts) {
  cat("[INFO] Found counts.txt - will create expression-based heatmaps\n")
} else {
  cat("[INFO] No counts.txt found - will show log2 fold change heatmaps\n")
}


## Function to create grouped heatmap
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    counts_file = NULL, output_dir = plots_dir) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  if (!file.exists(de_file)) {
    cat("[WARN] DE file not found: ", de_file, "\n")
    return(NULL)
  }

  de_data <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Find genes in data
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(de_data))

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

  ## Determine depot from contrast name
  depot <- NA
  if (grepl("iWAT", contrast_name)) depot <- "iWAT"
  if (grepl("gWAT", contrast_name)) depot <- "gWAT"

  ## Try to load expression data if counts file exists
  mat <- NULL
  col_labels <- NULL
  legend_title <- "log2FC"

  if (!is.null(counts_file) && file.exists(counts_file)) {
    cat("[INFO] Loading expression data from counts.txt\n")

    counts <- tryCatch({
      read.delim(counts_file, row.names = 1, check.names = FALSE)
    }, error = function(e) NULL)

    if (!is.null(counts) && all(genes_to_plot %in% rownames(counts))) {
      ## Log2 transform (add 1 pseudocount)
      log_counts <- log2(counts + 1)

      ## Get sample metadata from column names
      ## Expected format: JS_XX where samples can be identified
      ## For now, assume we need sample annotation - check if it's in the DE file path
      ## or we can infer from part1 parameters

      ## Simple approach: look for CTL vs KAT8KD patterns in sample names
      ## This is project-specific - adjust as needed
      sample_cols <- colnames(log_counts)

      ## Filter to genes of interest
      expr_mat <- log_counts[genes_to_plot, , drop = FALSE]

      ## Z-score per gene (across all samples)
      expr_z <- t(scale(t(expr_mat)))

      ## For a proper heatmap, we'd need sample group info
      ## For now, show per-sample z-scores
      mat <- expr_z
      legend_title <- "Z-score"
      cat("[INFO] Using z-scored expression data\n")
    }
  }

  ## Fallback: use log2FC from DE results
  if (is.null(mat)) {
    ## Get logFC column
    if ("log2FoldChange" %in% colnames(de_data)) {
      lfc_col <- "log2FoldChange"
    } else if ("logFC" %in% colnames(de_data)) {
      lfc_col <- "logFC"
    } else {
      cat("[WARN] Cannot find logFC column\n")
      return(NULL)
    }

    ## Single column showing fold change
    log2fc_values <- de_data[genes_to_plot, lfc_col]
    mat <- matrix(log2fc_values, ncol = 1)
    rownames(mat) <- genes_to_plot
    colnames(mat) <- "log2FC"
    legend_title <- "log2FC"
  }

  ## Order by category
  gene_order <- names(sort(gene_to_category))
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Color scale (symmetric, teal-white-orange)
  max_val <- max(abs(mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  max_val <- ceiling(max_val * 2) / 2  ## Round for cleaner legend

  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),
    c("#0072B2", "white", "#E69F00")
  )

  ## Adjust width based on number of columns
  heatmap_width <- unit(max(1.5, ncol(mat) * 0.4), "cm")

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = legend_title,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_split = factor(gene_to_category, levels = names(gene_categories)),
    row_gap = unit(1.5, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_title = contrast_name,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    rect_gp = gpar(col = "grey90", lwd = 0.3),
    width = heatmap_width,
    heatmap_legend_param = list(
      title = legend_title,
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(2.5, "cm"),
      at = c(-max_val, 0, max_val)
    )
  )

  ## Save
  plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = max(500, 300 + ncol(mat) * 80),
      height = max(400, n_genes * 18 + n_categories * 30),
      res = 150)
  draw(ht, padding = unit(c(2, 12, 2, 2), "mm"))
  dev.off()

  cat("[OK] Saved heatmap: ", plot_file, "\n", sep = "")

  return(ht)
}

## Create heatmaps for each contrast
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

for (de_file in de_files) {
  filename <- basename(de_file)
  contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
  if (!is.na(contrast)) {
    create_grouped_heatmap(
      de_file, contrast, gene_categories,
      counts_file = if (has_counts) counts_file else NULL
    )
  }
}


cat("\n======================================\n")
cat("=== PART 4 COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - GO:BP enrichment dot plots (top 15, ordered by p-value)\n")
cat("  - KEGG enrichment dot plots (top 15, ordered by p-value)\n")
cat("  - Grouped gene heatmaps\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
