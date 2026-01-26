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
## Nature / Cell Metabolism style:
## - Size = Gene Ratio
## - Color = -log10(FDR)
## - Most significant at TOP
## - Separate plots for Up/Down
## - GO:BP: top 15, KEGG: top 10
## ============================================================

cat("\n=== CREATING ORA ENRICHMENT DOT PLOTS ===\n")

create_enrichment_dotplot <- function(ora_file, contrast_name, direction,
                                       database = "GO:BP",
                                       output_dir = plots_dir) {

  if (!file.exists(ora_file)) {
    cat("[WARN] File not found: ", ora_file, "\n")
    return(NULL)
  }

  ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)

  if (nrow(ora_data) == 0) {
    cat("[WARN] No pathways in: ", basename(ora_file), "\n")
    return(NULL)
  }

  ## Set top_n based on database
  top_n <- if (database == "KEGG") 10 else 15

  ## Parse GeneRatio to numeric
  if ("GeneRatio" %in% colnames(ora_data) && is.character(ora_data$GeneRatio)) {
    ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      num <- as.numeric(x[1])
      denom <- as.numeric(x[2])
      if (is.na(denom) || denom == 0) return(0)
      num / denom
    })
  }

  ## Get gene count if not present
  if (!"Count" %in% colnames(ora_data)) {
    ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1])
    })
  }

  ## Calculate -log10(FDR) for color mapping
  ora_data$neg_log10_fdr <- -log10(ora_data$p.adjust)

  ## Filter and select top pathways
  ## - FDR < 0.05
  ## - Count >= 3
  ## - Rank by p.adjust (ascending), break ties by GeneRatio (descending)
  plot_data <- ora_data %>%
    dplyr::filter(p.adjust < 0.05, Count >= 3) %>%
    dplyr::arrange(p.adjust, desc(GeneRatio_numeric)) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_data) == 0) {
    cat("[WARN] No significant pathways meeting criteria (FDR<0.05, Count>=3)\n")
    return(NULL)
  }

  ## Clip extreme -log10(FDR) values to prevent one pathway dominating color scale
  ## Cap at 99th percentile or max of 10, whichever is smaller
  fdr_cap <- min(quantile(plot_data$neg_log10_fdr, 0.99, na.rm = TRUE), 10)
  plot_data$neg_log10_fdr_clipped <- pmin(plot_data$neg_log10_fdr, fdr_cap)

  ## Order y-axis: most significant at TOP
  ## Factor levels in reverse order of p-value (smallest p = top)
  plot_data <- plot_data %>%
    dplyr::arrange(desc(p.adjust)) %>%
    dplyr::mutate(
      Description = str_wrap(Description, width = 45),
      Description = factor(Description, levels = Description)
    )

  ## Direction-specific color gradient (cool to warm, light to dark)
  if (direction == "Up") {
    ## Orange gradient for upregulated
    color_low <- "#FEE8C8"   ## Light
    color_high <- "#B30000"  ## Dark red-orange
  } else {
    ## Blue gradient for downregulated
    color_low <- "#DEEBF7"   ## Light
    color_high <- "#08306B"  ## Dark blue
  }

  ## Build title in journal format
  ## Example: "GO:BP Enrichment (KAT8KD vs CTL, iWAT, Upregulated)"
  direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
  plot_title <- paste0(database, " Enrichment\n(", contrast_name, ", ", direction_label, ")")

  ## Create dot plot
  ## Size = GeneRatio, Color = -log10(FDR)
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = GeneRatio_numeric, color = neg_log10_fdr_clipped)) +
    scale_color_gradient(
      low = color_low,
      high = color_high,
      name = expression(-log[10](FDR)),
      limits = c(min(plot_data$neg_log10_fdr_clipped), fdr_cap)
    ) +
    scale_size_continuous(
      name = "Gene Ratio",
      range = c(2, 6),
      breaks = pretty(plot_data$GeneRatio_numeric, n = 3)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.1)),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    labs(
      title = plot_title,
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      ## Title
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11, lineheight = 1.1),
      ## Axes
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 5)),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      ## Grid - minimal
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      ## Border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      ## Legend
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      ## Margins
      plot.margin = margin(10, 10, 10, 10)
    ) +
    guides(
      color = guide_colorbar(order = 1, barwidth = 0.8, barheight = 3),
      size = guide_legend(order = 2)
    )

  ## Save
  db_short <- gsub(":", "", database)
  plot_file <- paste0("Dotplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 7.5,
    height = max(4, nrow(plot_data) * 0.25 + 1.5),
    dpi = 300,
    bg = "white"
  )
  cat("[OK] Saved: ", plot_file, " (", nrow(plot_data), " pathways)\n", sep = "")

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
## SECTION 2: GROUPED HEATMAPS WITH TWO COLUMNS
## ============================================================
## Shows mean expression per group (CTL vs KAT8KD)
## Z-score normalized across both groups
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## ============================================================
## >>> ADD YOUR CUSTOM GENES HERE <<<
## ============================================================

gene_categories <- list(

  "Inflammatory\nResponse" = c(
    "Il1b", "Il6", "Tnf", "Il1a"
  ),

  "Chemotaxis" = c(
    "Ccl2", "Ccl7", "Cxcl12"
  ),

  "ECM /\nCollagens" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
    "Col5a1", "Col5a3", "Col6a1", "Col6a2", "Col6a3",
    "Col6a6", "Col15a1"
  ),

  "Matrix\nRemodeling" = c(
    "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
    "Timp1", "Timp2", "Timp3", "Timp4"
  ),

  "Adipocyte\nMarkers" = c(
    "Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1"
  )

)

## ============================================================
## >>> END OF CUSTOM GENE SECTION <<<
## ============================================================


## Sample annotation (matching part1_main_analysis.R)
sample_ids <- paste0("JS_", sprintf("%02d", 1:40))
sample_annot <- data.frame(
  Sample   = sample_ids,
  Depot    = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
  Sex      = c(
    rep("F", 5), rep("F", 5),
    rep("M", 5), rep("M", 5),
    rep("F", 5), rep("F", 5),
    rep("M", 5), rep("M", 5)
  ),
  Genotype = c(
    rep("CTL", 5), rep("KAT8KD", 5),
    rep("CTL", 5), rep("KAT8KD", 5),
    rep("CTL", 5), rep("KAT8KD", 5),
    rep("CTL", 5), rep("KAT8KD", 5)
  ),
  stringsAsFactors = FALSE
)
rownames(sample_annot) <- sample_annot$Sample

## Outliers to exclude
outlier_samples <- c("JS_08", "JS_28")


## Check if counts.txt exists
counts_file <- file.path(getwd(), "counts.txt")
has_counts <- file.exists(counts_file)

if (has_counts) {
  cat("[INFO] Found counts.txt - loading expression data\n")
  counts_raw <- tryCatch({
    read.delim(counts_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    cat("[WARN] Failed to load counts.txt: ", conditionMessage(e), "\n")
    NULL
  })
} else {
  cat("[INFO] No counts.txt found - will use log2FC values\n")
  counts_raw <- NULL
}


## Function to create grouped heatmap with two columns
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    counts_raw = NULL, sample_annot = NULL,
                                    outliers = c(), output_dir = plots_dir) {

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

  mat <- NULL
  legend_title <- "Z-score"

  ## Try to use expression data if available
  if (!is.null(counts_raw) && !is.null(sample_annot) && !is.na(depot)) {

    ## Check if gene_id column exists
    if ("gene_id" %in% colnames(counts_raw)) {
      rownames(counts_raw) <- counts_raw$gene_id
    }

    ## Check if genes are in counts
    genes_in_counts <- intersect(genes_to_plot, rownames(counts_raw))

    if (length(genes_in_counts) >= 3) {
      cat("[INFO] Using expression data for ", length(genes_in_counts), " genes\n", sep = "")

      ## Get samples for this depot (excluding outliers)
      depot_samples <- sample_annot %>%
        dplyr::filter(Depot == depot, !Sample %in% outliers)

      ctl_samples <- depot_samples %>% dplyr::filter(Genotype == "CTL") %>% dplyr::pull(Sample)
      kd_samples <- depot_samples %>% dplyr::filter(Genotype == "KAT8KD") %>% dplyr::pull(Sample)

      ## Check samples exist in counts
      ctl_samples <- intersect(ctl_samples, colnames(counts_raw))
      kd_samples <- intersect(kd_samples, colnames(counts_raw))

      if (length(ctl_samples) >= 2 && length(kd_samples) >= 2) {

        ## Extract expression for genes of interest
        expr_raw <- counts_raw[genes_in_counts, c(ctl_samples, kd_samples), drop = FALSE]
        expr_raw <- as.matrix(expr_raw)

        ## Log2 transform (pseudocount of 1)
        expr_log <- log2(expr_raw + 1)

        ## Calculate group means
        ctl_mean <- rowMeans(expr_log[, ctl_samples, drop = FALSE], na.rm = TRUE)
        kd_mean <- rowMeans(expr_log[, kd_samples, drop = FALSE], na.rm = TRUE)

        ## Z-score across both group means
        combined <- cbind(ctl_mean, kd_mean)
        row_means <- rowMeans(combined, na.rm = TRUE)
        row_sds <- apply(combined, 1, sd, na.rm = TRUE)
        row_sds[row_sds == 0] <- 1  ## Avoid division by zero

        mat <- (combined - row_means) / row_sds
        colnames(mat) <- c("CTL", "KAT8KD")

        ## Update genes_to_plot to only include those in counts
        genes_to_plot <- genes_in_counts
        gene_to_category <- gene_to_category[genes_to_plot]

        cat("[INFO] Created z-scored group means heatmap\n")
      }
    }
  }

  ## Fallback: use log2FC with pseudo two-column approach
  if (is.null(mat)) {
    cat("[INFO] Using log2FC-based heatmap\n")

    if ("log2FoldChange" %in% colnames(de_data)) {
      lfc_col <- "log2FoldChange"
    } else if ("logFC" %in% colnames(de_data)) {
      lfc_col <- "logFC"
    } else {
      cat("[WARN] Cannot find logFC column\n")
      return(NULL)
    }

    ## Create two columns: CTL = -FC/2, KD = +FC/2
    ## This centers the data so both columns show color
    log2fc_values <- de_data[genes_to_plot, lfc_col]
    mat <- cbind(
      CTL = -log2fc_values / 2,
      KAT8KD = log2fc_values / 2
    )
    rownames(mat) <- genes_to_plot
    legend_title <- "Relative\nExpression"
  }

  ## Order by category
  gene_order <- names(sort(gene_to_category))
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Color scale (symmetric, teal-white-orange)
  max_val <- max(abs(mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  max_val <- ceiling(max_val * 2) / 2

  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),
    c("#0072B2", "white", "#E69F00")
  )

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = legend_title,
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
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title = contrast_name,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    rect_gp = gpar(col = "white", lwd = 0.5),
    width = unit(2.5, "cm"),
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
      width = 600,
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
      counts_raw = counts_raw,
      sample_annot = sample_annot,
      outliers = outlier_samples
    )
  }
}


cat("\n======================================\n")
cat("=== PART 4 COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - GO:BP enrichment dot plots\n")
cat("  - KEGG enrichment dot plots\n")
cat("  - Grouped gene heatmaps (CTL vs KAT8KD)\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
