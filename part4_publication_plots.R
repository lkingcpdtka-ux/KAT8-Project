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

ora_dotplot_params <- list(
  fdr_cutoff  = 0.05,
  min_count   = 3,
  top_n       = list("GO:BP" = 15, "KEGG" = 15),
  default_top_n = 15,
  order_by    = "adj. p-value",
  tie_breaker = "gene ratio",
  gobp_simplify_cutoff = 0.7  ## GO:BP simplification cutoff (from part2_ora.R)
)

heatmap_params <- list(
  lfc_clip    = 2.0,
  lfc_label   = "log2FC (DESeq2)",
  reference   = "CTL=0 reference"
)


## ============================================================
## SECTION 1: DOT PLOTS FOR ORA PATHWAY ENRICHMENT
## ============================================================
## Style:
## - Size = Gene Count (number of genes)
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
  top_n <- ora_dotplot_params$top_n[[database]]
  if (is.null(top_n)) {
    top_n <- ora_dotplot_params$default_top_n
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

  ## Get gene count if not present
  if (!"Count" %in% colnames(ora_data)) {
    ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1])
    })
  }

  ## Calculate -log10(FDR) for color mapping
  ora_data$neg_log10_fdr <- -log10(ora_data$p.adjust)

  ## Determine thresholds from file metadata if available
  plot_fdr_cutoff <- ora_dotplot_params$fdr_cutoff
  if ("fdr_cutoff" %in% colnames(ora_data)) {
    fdr_vals <- unique(na.omit(ora_data$fdr_cutoff))
    if (length(fdr_vals) > 0) {
      plot_fdr_cutoff <- fdr_vals[1]
    }
  }

  plot_logfc_cutoff <- NULL
  if ("logFC_cutoff" %in% colnames(ora_data)) {
    logfc_vals <- unique(na.omit(ora_data$logFC_cutoff))
    if (length(logfc_vals) > 0) {
      plot_logfc_cutoff <- logfc_vals[1]
    }
  }

  ## Filter and select top pathways
  ## - FDR < cutoff
  ## - Count >= 3 (minimum genes for reliable enrichment)
  ## - Rank by p.adjust (ascending), break ties by GeneRatio (descending)
  plot_data <- ora_data %>%
    dplyr::filter(p.adjust < plot_fdr_cutoff, Count >= ora_dotplot_params$min_count) %>%
    dplyr::arrange(p.adjust, desc(GeneRatio_numeric)) %>%
    dplyr::slice_head(n = top_n)

  if (nrow(plot_data) == 0) {
    cat(
      "[WARN] No significant pathways meeting criteria (FDR<",
      plot_fdr_cutoff,
      ", Count>=",
      ora_dotplot_params$min_count,
      ")\n",
      sep = ""
    )
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
    color_low <- "#FDD49E"   ## Light orange
    color_high <- "#E69F00"  ## Orange
  } else {
    ## Blue gradient for downregulated
    color_low <- "#9ECAE1"   ## Light blue
    color_high <- "#0072B2"  ## Blue
  }

  ## Build title in journal format
  direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
  plot_title <- paste0(database, " Enrichment\n(", contrast_name, ", ", direction_label, ")")

  ## Parameter caption for bottom of plot
  ## Add simplification info for GO:BP
  simplify_label <- if (database == "GO:BP") {
    paste0(" | Simplified (cutoff=", ora_dotplot_params$gobp_simplify_cutoff, ")")
  } else {
    ""
  }
  if (!is.null(plot_logfc_cutoff)) {
    param_caption <- paste(
      paste0(
        "FDR < ", plot_fdr_cutoff,
        " | |log2FC| > ", plot_logfc_cutoff,
        " | Top ", top_n, " pathways", simplify_label
      ),
      paste0("Ordered by ", ora_dotplot_params$order_by),
      sep = "\n"
    )
  } else {
    param_caption <- paste(
      paste0(
        "FDR < ", plot_fdr_cutoff,
        " | Top ", top_n, " pathways", simplify_label
      ),
      paste0("Ordered by ", ora_dotplot_params$order_by),
      sep = "\n"
    )
  }

  ## Create dot plot
  ## Size = Count, Color = -log10(FDR)
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = neg_log10_fdr_clipped)) +
    scale_color_gradient(
      low = color_low,
      high = color_high,
      name = expression(-log[10](FDR)),
      limits = c(min(plot_data$neg_log10_fdr_clipped), fdr_cap)
    ) +
    scale_size_continuous(
      name = "Gene Count",
      range = c(2, 6),
      breaks = pretty(plot_data$Count, n = 4)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      title = plot_title,
      x = "Gene Ratio",
      y = NULL,
      caption = param_caption
    ) +
    theme_bw(base_size = 10) +
    theme(
      ## Title
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11, lineheight = 1.1),
      ## Axes
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.title.x = element_text(size = 8, margin = margin(t = 5)),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      ## Grid - light lines for both axes
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      ## Border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      ## Legend
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      ## Caption (parameters)
      plot.caption = element_text(size = 7, color = "grey40", hjust = 0.5, margin = margin(t = 8)),
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
    width = 6.5,
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


## Function to create grouped heatmap with log2FC coloring
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    lfc_clip = heatmap_lfc_clip,
                                    output_dir = plots_dir) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  if (!file.exists(de_file)) {
    cat("[WARN] DE file not found: ", de_file, "\n")
    return(NULL)
  }

  de_data <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Get log2FC column
  if ("log2FoldChange" %in% colnames(de_data)) {
    lfc_col <- "log2FoldChange"
  } else if ("logFC" %in% colnames(de_data)) {
    lfc_col <- "logFC"
  } else {
    cat("[WARN] Cannot find log2FC column\n")
    return(NULL)
  }

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

  ## Create matrix using log2FC values
  ## CTL = 0 (reference), KAT8KD = log2FC
  log2fc_values <- de_data[genes_to_plot, lfc_col]
  mat <- cbind(
    CTL = rep(0, length(genes_to_plot)),
    KAT8KD = log2fc_values
  )
  rownames(mat) <- genes_to_plot

  ## Order by category
  gene_order <- c()
  for (cat_name in names(gene_categories)) {
    cat_genes <- names(gene_to_category[gene_to_category == cat_name])
    gene_order <- c(gene_order, cat_genes)
  }
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Clip values at ±lfc_clip
  mat[mat > lfc_clip] <- lfc_clip
  mat[mat < -lfc_clip] <- -lfc_clip

  ## Color scale: blue-white-orange, centered at 0
  col_fun <- colorRamp2(
    c(-lfc_clip, 0, lfc_clip),
    c("#2166AC", "white", "#E69F00")  ## Blue-white-orange
  )

  ## Parameter caption
  param_caption <- paste0(
    heatmap_params$lfc_label,
    " | ",
    heatmap_params$reference,
    " | Clipped at \u00b1",
    lfc_clip
  )

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = "log2FC",
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
      title = "log2 fold change\n(DESeq2)",
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(2.5, "cm"),
      at = c(-lfc_clip, 0, lfc_clip),
      labels = c(paste0("\u2264", -lfc_clip), "0", paste0("\u2265", lfc_clip))
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
cat("=== PART 4 COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - GO:BP enrichment dot plots\n")
cat("  - KEGG enrichment dot plots\n")
cat("  - Grouped gene heatmaps (CTL vs KAT8KD)\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
