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
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "ComplexHeatmap",
  "circlize",
  "SummarizedExperiment"  ## For loading VST data
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
  library(viridis)  ## For viridis color scale
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

  ## Standardize column names - handle both 'Description' and 'pathway_name'
  if (!"Description" %in% colnames(ora_data) && "pathway_name" %in% colnames(ora_data)) {
    ora_data$Description <- ora_data$pathway_name
    cat("[INFO] Using 'pathway_name' column as Description\n")
  }

  if (!"Description" %in% colnames(ora_data)) {
    cat("[WARN] No Description or pathway_name column found in: ", basename(ora_file), "\n")
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
  ## Cap at 99th percentile or max of 15, whichever is smaller
  fdr_cap <- min(quantile(plot_data$neg_log10_fdr, 0.99, na.rm = TRUE), 15)
  fdr_cap <- max(fdr_cap, 3)  ## Ensure minimum range for color scale
  plot_data$neg_log10_fdr_clipped <- pmin(plot_data$neg_log10_fdr, fdr_cap)

  ## Order y-axis: most significant at TOP
  ## Factor levels in reverse order of p-value (smallest p = top)
  plot_data <- plot_data %>%
    dplyr::arrange(desc(p.adjust)) %>%
    dplyr::mutate(
      Description = str_wrap(Description, width = 45),
      Description = factor(Description, levels = Description)
    )

  ## Viridis color scale (higher -log10 p-value = more significant = darker/purple)
  ## This is cleaner than direction-specific colors and universally recognized

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
  ## Size = Count, Color = -log10(FDR) using mako scale
  ## Higher -log10(FDR) = more significant = darker in mako

  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = neg_log10_fdr_clipped)) +
    scale_color_viridis_c(
      option = "mako",
      direction = -1,  ## Higher values = darker (more significant)
      name = expression(-log[10](FDR)),
      limits = c(0, ceiling(fdr_cap)),  ## Start at 0 for proper scale
      breaks = scales::pretty_breaks(n = 4)
    ) +
    scale_size_continuous(
      name = "Gene Count",
      range = c(2.5, 6),  ## Smaller dots to prevent clipping
      breaks = pretty(plot_data$Count, n = 4)
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.08, 0.08))  ## More padding for dots at edges
    ) +
    labs(
      title = plot_title,
      x = "Gene Ratio",
      y = NULL,
      caption = param_caption
    ) +
    theme_bw(base_size = 12) +
    theme(
      ## Title
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13, lineheight = 1.1),
      ## Axes - larger fonts
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, margin = margin(t = 5)),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      ## Grid - light lines for both axes
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      ## Border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      ## Legend - larger fonts
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      ## Caption (parameters)
      plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8)),
      ## Margins - tighter for less whitespace
      plot.margin = margin(8, 5, 8, 5)
    ) +
    guides(
      color = guide_colorbar(order = 1, barwidth = 1, barheight = 3.5),
      size = guide_legend(order = 2)
    )

  ## Save - balanced width (not too wide, not clipping)
  db_short <- gsub(":", "", database)
  plot_file <- paste0("Dotplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 5.5,  ## Balanced width to avoid clipping
    height = max(4, nrow(plot_data) * 0.3 + 1.5),
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
## SECTION 2: GROUPED HEATMAPS WITH VST Z-SCORES
## ============================================================
## Heatmap style (journal-ready):
## - Two columns: mean CTL and mean KAT8KD expression
## - Rows grouped by custom gene categories
## - Values: Z-scored VST expression (like Part 1 heatmaps)
## - Color scale: centered at 0, diverging scale
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## ============================================================
## LOAD VST DATA FROM PART 1
## ============================================================
cat("[INFO] Loading VST data from Part 1...\n")

## Load cached VST matrix
cache_dir <- file.path(outdir, "cache")
vst_cache_file <- file.path(cache_dir, "vst_tissue_heatmap.rds")

if (!file.exists(vst_cache_file)) {
  cat("[WARN] VST cache not found: ", vst_cache_file, "\n")
  cat("[WARN] Falling back to log2FC-based heatmaps\n")
  use_vst <- FALSE
} else {
  vst_tissue <- readRDS(vst_cache_file)
  vst_mat <- SummarizedExperiment::assay(vst_tissue)
  vst_coldata <- as.data.frame(SummarizedExperiment::colData(vst_tissue))
  cat("[OK] Loaded VST matrix: ", nrow(vst_mat), " genes x ", ncol(vst_mat), " samples\n", sep = "")
  use_vst <- TRUE
}

## ============================================================
## HEATMAP SETTINGS
## ============================================================
heatmap_zscore_clip <- 2.5  ## Clip z-scores at ±2.5 for visualization
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


## Function to create grouped heatmap with VST z-scores (journal style)
## Shows mean expression per treatment group
create_grouped_heatmap <- function(contrast_name, gene_categories,
                                    zscore_clip = heatmap_zscore_clip,
                                    output_dir = plots_dir) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  if (!use_vst) {
    cat("[WARN] VST data not available, skipping heatmap\n")
    return(NULL)
  }

  ## Determine depot from contrast name
  depot <- NA_character_
  if (grepl("^iWAT", contrast_name)) depot <- "iWAT"
  else if (grepl("^gWAT", contrast_name)) depot <- "gWAT"

  if (is.na(depot)) {
    cat("[WARN] Could not determine depot from contrast: ", contrast_name, "\n")
    return(NULL)
  }

  ## Get samples for this depot
  depot_samples <- rownames(vst_coldata)[vst_coldata$Depot == depot]
  depot_meta <- vst_coldata[depot_samples, , drop = FALSE]

  cat("[INFO] Found ", length(depot_samples), " samples for ", depot, "\n", sep = "")

  ## Find genes in VST data
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(vst_mat))

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

  ## Get VST expression for these genes
  vst_subset <- vst_mat[genes_to_plot, depot_samples, drop = FALSE]

  ## Calculate mean expression per treatment group
  ctl_samples <- depot_samples[depot_meta$Genotype == "CTL"]
  kd_samples <- depot_samples[depot_meta$Genotype == "KAT8KD"]

  mean_ctl <- rowMeans(vst_subset[, ctl_samples, drop = FALSE])
  mean_kd <- rowMeans(vst_subset[, kd_samples, drop = FALSE])

  ## Create matrix with mean values
  mat_means <- cbind(CTL = mean_ctl, KAT8KD = mean_kd)

  ## Z-score across rows (each gene scaled to mean=0, sd=1)
  ## This shows relative expression difference between conditions
  mat_zscore <- t(scale(t(mat_means)))

  ## Order by category
  gene_order <- c()
  for (cat_name in names(gene_categories)) {
    cat_genes <- names(gene_to_category[gene_to_category == cat_name])
    gene_order <- c(gene_order, cat_genes)
  }
  mat_zscore <- mat_zscore[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Clip z-scores for visualization
  mat_zscore[mat_zscore > zscore_clip] <- zscore_clip
  mat_zscore[mat_zscore < -zscore_clip] <- -zscore_clip

  ## Color scale: mako-inspired diverging scale centered at 0
  col_fun <- colorRamp2(
    c(-zscore_clip, 0, zscore_clip),
    c("#0B0405", "white", "#DEF5E5")  ## Mako dark -> white -> mako light
  )

  ## Parameter caption
  param_caption <- paste0(
    "VST expression (z-scored) | ",
    depot, " samples | ",
    "n=", length(ctl_samples), " CTL, n=", length(kd_samples), " KAT8KD"
  )

  ## Create heatmap
  ht <- Heatmap(
    mat_zscore,
    col = col_fun,
    name = "Z-score",
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
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = paste0(contrast_name, "\n(Mean VST Expression)"),
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.8),
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    width = unit(2.8, "cm"),
    heatmap_legend_param = list(
      title = "Relative\nExpression\n(Z-score)",
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(2.5, "cm"),
      at = c(-zscore_clip, 0, zscore_clip),
      labels = c(paste0("\u2264", -zscore_clip), "0", paste0("\u2265", zscore_clip))
    )
  )

  ## Save
  plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = 650,
      height = max(400, n_genes * 18 + n_categories * 30 + 80),
      res = 150)
  draw(ht, padding = unit(c(5, 12, 12, 5), "mm"))
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
    create_grouped_heatmap(contrast, gene_categories)
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
