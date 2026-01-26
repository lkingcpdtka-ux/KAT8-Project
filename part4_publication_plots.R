#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: PUBLICATION QUALITY PLOTS
## =========================================================
## This script creates publication-quality visualizations:
## 1. Dot plots for ORA pathway enrichment
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

cat("\n=== CREATING ORA ENRICHMENT DOT PLOTS ===\n")

create_enrichment_dotplot <- function(ora_file, contrast_name, direction,
                                       top_n = 20, output_dir = plots_dir) {

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

  ## Select top pathways
  plot_data <- ora_data %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_data) == 0) {
    cat("[WARN] No significant pathways for dotplot\n")
    return(NULL)
  }

  ## Clean descriptions
  plot_data <- plot_data %>%
    dplyr::mutate(
      Description = str_wrap(Description, width = 50),
      Description = factor(Description, levels = rev(Description))
    )

  ## Direction-specific colors
  if (direction == "Up") {
    color_low <- "#FEE0D2"
    color_high <- "#CB181D"
  } else {
    color_low <- "#DEEBF7"
    color_high <- "#08519C"
  }

  ## Create dot plot
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(
      low = color_high,
      high = color_low,
      name = "Adj. P-value",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(
      name = "Gene\nNumber",
      range = c(4, 12),
      breaks = pretty(plot_data$Count, n = 4)
    ) +
    labs(
      title = paste0("GO Enrichment: ", contrast_name, " (", direction, "-regulated)"),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )

  ## Save
  plot_file <- paste0("Dotplot_GO_", contrast_name, "_", direction, "_", run_tag, ".png")
  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 10,
    height = max(6, nrow(plot_data) * 0.35),
    dpi = 300,
    bg = "white"
  )
  cat("[OK] Saved dot plot: ", plot_file, "\n", sep = "")

  return(p)
}

## Process all ORA GO:BP results
ora_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)

for (ora_file in ora_files) {
  filename <- basename(ora_file)
  parts <- str_match(filename, "ORA_gobp_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
  if (!is.na(parts[1])) {
    contrast <- parts[2]
    direction <- parts[3]
    create_enrichment_dotplot(ora_file, contrast, direction)
  }
}


## ============================================================
## SECTION 2: GROUPED HEATMAPS WITH CUSTOM GENE CATEGORIES
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## ============================================================
## >>> ADD YOUR CUSTOM GENES HERE <<<
## ============================================================
##
## Instructions:
## 1. Each category is a named list element
## 2. The NAME becomes the row label in the heatmap
## 3. Use \n for line breaks in category names
## 4. Add gene symbols exactly as they appear in your DE results
##
## Example format:
##   "Category\nName" = c("Gene1", "Gene2", "Gene3", ...)
##
## ============================================================

gene_categories <- list(

  ## ------------------------------------------------------------
  ## INFLAMMATION GENES - add/remove genes as needed
  ## ------------------------------------------------------------
  "Inflammatory\nResponse" = c(
    "Il1b", "Il6", "Tnf", "Il1a"
    ## ADD MORE INFLAMMATION GENES HERE
  ),

 ## ------------------------------------------------------------
  ## CHEMOTAXIS / MIGRATION GENES
  ## ------------------------------------------------------------
  "Chemotaxis" = c(
    "Ccl2", "Ccl7", "Cxcl12"
    ## ADD MORE CHEMOTAXIS GENES HERE
  ),

  ## ------------------------------------------------------------
  ## ECM / COLLAGEN GENES
  ## ------------------------------------------------------------
  "ECM /\nCollagens" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
    "Col5a1", "Col5a3", "Col6a1", "Col6a2", "Col6a3",
    "Col6a6", "Col15a1"
    ## ADD MORE ECM GENES HERE
  ),

  ## ------------------------------------------------------------
  ## MATRIX REMODELING (MMPs, TIMPs)
  ## ------------------------------------------------------------
  "Matrix\nRemodeling" = c(
    "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
    "Timp1", "Timp2", "Timp3", "Timp4"
    ## ADD MORE MMP/TIMP GENES HERE
  ),

  ## ------------------------------------------------------------
  ## ADIPOCYTE MARKERS
  ## ------------------------------------------------------------
  "Adipocyte\nMarkers" = c(
    "Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1"
    ## ADD MORE ADIPOCYTE GENES HERE
  )

  ## ------------------------------------------------------------
  ## ADD MORE CATEGORIES BELOW
  ## ------------------------------------------------------------
  ##
  ## "Your Category\nName" = c(
  ##   "Gene1", "Gene2", "Gene3"
  ## )

)

## ============================================================
## >>> END OF CUSTOM GENE SECTION <<<
## ============================================================


## Function to create grouped heatmap
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    output_dir = plots_dir) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  if (!file.exists(de_file)) {
    cat("[WARN] DE file not found: ", de_file, "\n")
    return(NULL)
  }

  de_data <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Find genes in data
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(de_data))

  cat("[INFO] Category genes found in data: ", length(genes_in_data), "/",
      length(all_category_genes), "\n", sep = "")

  ## Show which genes are missing
  missing_genes <- setdiff(all_category_genes, rownames(de_data))
  if (length(missing_genes) > 0 && length(missing_genes) <= 20) {
    cat("[INFO] Missing genes: ", paste(missing_genes, collapse = ", "), "\n", sep = "")
  } else if (length(missing_genes) > 20) {
    cat("[INFO] Missing ", length(missing_genes), " genes (too many to list)\n", sep = "")
  }

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

  ## Get logFC column
  if ("log2FoldChange" %in% colnames(de_data)) {
    lfc_col <- "log2FoldChange"
  } else if ("logFC" %in% colnames(de_data)) {
    lfc_col <- "logFC"
  } else {
    cat("[WARN] Cannot find logFC column\n")
    return(NULL)
  }

  ## Get FDR column
  if ("padj" %in% colnames(de_data)) {
    fdr_col <- "padj"
  } else if ("adj.P.Val" %in% colnames(de_data)) {
    fdr_col <- "adj.P.Val"
  } else {
    fdr_col <- NULL
  }

  ## Create matrix
  mat <- as.matrix(de_data[genes_to_plot, lfc_col, drop = FALSE])
  colnames(mat) <- "log2FC"

  ## Order by category
  gene_order <- names(sort(gene_to_category))
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Significance markers
  sig_markers <- rep("", length(gene_order))
  if (!is.null(fdr_col)) {
    fdr_vals <- de_data[gene_order, fdr_col]
    sig_markers <- ifelse(is.na(fdr_vals), "",
                          ifelse(fdr_vals < 0.001, "***",
                                 ifelse(fdr_vals < 0.01, "**",
                                        ifelse(fdr_vals < 0.05, "*", ""))))
  }

  ## Color scale
  max_val <- max(abs(mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),
    c("#0072B2", "white", "#E69F00")
  )

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = "log2FC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_split = factor(gene_to_category, levels = names(gene_categories)),
    row_gap = unit(3, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 9),
    column_title = paste0(contrast_name),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    border = TRUE,
    width = unit(3, "cm"),
    heatmap_legend_param = list(
      title = "log2FC",
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(4, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_markers[i] != "") {
        grid.text(sig_markers[i], x, y, gp = gpar(fontsize = 8, fontface = "bold"))
      }
    }
  )

  ## Save
  plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = 800,
      height = max(600, n_genes * 22 + n_categories * 40),
      res = 150)
  draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
  dev.off()

  cat("[OK] Saved grouped heatmap: ", plot_file, "\n", sep = "")
  cat("[INFO] Genes plotted: ", n_genes, " in ", n_categories, " categories\n", sep = "")

  return(ht)
}

## Create heatmaps for each contrast
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

for (de_file in de_files) {
  filename <- basename(de_file)
  contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
  if (!is.na(contrast)) {
    create_grouped_heatmap(de_file, contrast, gene_categories)
  }
}


cat("\n======================================\n")
cat("=== PART 4 (Publication Plots) COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - ORA enrichment dot plots\n")
cat("  - Grouped gene heatmaps\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
cat("To customize gene categories, edit the gene_categories list\n")
cat("starting at line ~165 in this script.\n\n")
