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

  ## Parameter caption for bottom of plot
  param_caption <- paste0(
    "Cutoffs: FDR < 0.05 | Count \u2265 3 | Top ", top_n, " pathways | ",
    "Ordered by adj. p-value"
  )

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
      y = NULL,
      caption = param_caption
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
## SECTION 2: PATHWAY-BASED HEATMAPS WITH LOG2FC
## ============================================================
## Heatmap style (journal-ready):
## - Two columns: CTL (reference=0) and KAT8KD (log2FC)
## - Rows grouped by pathway blocks from ORA results
## - Values: DESeq2 log2FoldChange (NOT z-scores or expression)
## - Color scale: centered at 0, clipped at ±2
## - Pathway order: by FDR ascending (most significant first)
## - Gene order within pathway: padj ascending, then |log2FC| descending
## ============================================================

cat("\n=== CREATING PATHWAY-BASED HEATMAPS ===\n")

## ============================================================
## HEATMAP SETTINGS
## ============================================================
heatmap_lfc_clip <- 2.0       ## Clip log2FC at ±2 (use 1.5 for subtle effects)
heatmap_top_pathways <- 8     ## Number of pathways to show
heatmap_max_genes_per_pathway <- 15  ## Max genes per pathway block
## ============================================================

## Function to create pathway-based heatmap using ORA results
create_pathway_heatmap <- function(contrast_name, direction,
                                   database = "GO:BP",
                                   top_pathways = heatmap_top_pathways,
                                   max_genes = heatmap_max_genes_per_pathway,
                                   lfc_clip = heatmap_lfc_clip,
                                   output_dir = plots_dir) {

  cat("\n--- Creating pathway heatmap: ", contrast_name, " ", direction, " (", database, ") ---\n", sep = "")

  ## Find ORA file
  db_short <- tolower(gsub(":", "", database))
  ora_pattern <- paste0("^ORA_", db_short, "_", contrast_name, "_", direction, "_.*\\.csv$")
  ora_files <- list.files(tables_dir, pattern = ora_pattern, full.names = TRUE)

  if (length(ora_files) == 0) {
    cat("[WARN] No ORA file found for pattern: ", ora_pattern, "\n")
    return(NULL)
  }
  ora_file <- ora_files[1]

  ## Find DE file
  de_pattern <- paste0("^DE_tissue_", contrast_name, "_.*\\.csv$")
  de_files <- list.files(tables_dir, pattern = de_pattern, full.names = TRUE)

  if (length(de_files) == 0) {
    cat("[WARN] No DE file found for pattern: ", de_pattern, "\n")
    return(NULL)
  }
  de_file <- de_files[1]

  ## Load data
  ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
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

  ## Get padj column
  if ("padj" %in% colnames(de_data)) {
    padj_col <- "padj"
  } else if ("adj.P.Val" %in% colnames(de_data)) {
    padj_col <- "adj.P.Val"
  } else {
    padj_col <- NULL
  }

  ## Filter significant pathways and select top by FDR
  sig_pathways <- ora_data %>%
    dplyr::filter(p.adjust < 0.05, Count >= 3) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_pathways)

  if (nrow(sig_pathways) == 0) {
    cat("[WARN] No significant pathways (FDR<0.05, Count>=3)\n")
    return(NULL)
  }

  cat("[INFO] Using top ", nrow(sig_pathways), " pathways\n", sep = "")

  ## Build gene-to-pathway mapping
  gene_pathway_list <- list()
  pathway_order <- c()

  for (i in seq_len(nrow(sig_pathways))) {
    pathway_name <- sig_pathways$Description[i]
    pathway_fdr <- sig_pathways$p.adjust[i]
    gene_str <- sig_pathways$geneID[i]

    if (is.na(gene_str) || gene_str == "") next

    ## Parse genes (typically "/" separated)
    genes <- unlist(strsplit(gene_str, "/"))
    genes <- trimws(genes)
    genes <- genes[genes %in% rownames(de_data)]

    if (length(genes) == 0) next

    ## Create gene info with ordering metrics
    gene_info <- data.frame(
      gene = genes,
      log2FC = de_data[genes, lfc_col],
      stringsAsFactors = FALSE
    )

    if (!is.null(padj_col)) {
      gene_info$padj <- de_data[genes, padj_col]
    } else {
      gene_info$padj <- 1
    }

    gene_info$abs_log2FC <- abs(gene_info$log2FC)

    ## Order genes: padj ascending, then |log2FC| descending
    gene_info <- gene_info %>%
      dplyr::arrange(padj, desc(abs_log2FC)) %>%
      dplyr::slice_head(n = max_genes)

    ## Wrap long pathway names
    pathway_label <- str_wrap(pathway_name, width = 35)

    gene_pathway_list[[pathway_label]] <- gene_info
    pathway_order <- c(pathway_order, pathway_label)
  }

  if (length(gene_pathway_list) == 0) {
    cat("[WARN] No genes found in pathways\n")
    return(NULL)
  }

  ## Build matrix: genes as rows, CTL and KAT8KD as columns
  ## CTL = 0 (reference), KAT8KD = log2FC
  all_genes <- c()
  gene_to_pathway <- c()
  lfc_values <- c()

  for (pathway in pathway_order) {
    gene_info <- gene_pathway_list[[pathway]]
    all_genes <- c(all_genes, gene_info$gene)
    gene_to_pathway <- c(gene_to_pathway, rep(pathway, nrow(gene_info)))
    lfc_values <- c(lfc_values, gene_info$log2FC)
  }

  ## Create matrix
  mat <- cbind(
    CTL = rep(0, length(all_genes)),
    KAT8KD = lfc_values
  )
  rownames(mat) <- all_genes

  ## Handle duplicate gene names (can occur across pathways)
  if (any(duplicated(rownames(mat)))) {
    rownames(mat) <- make.unique(rownames(mat), sep = "_")
  }

  ## Clip values at ±lfc_clip
  mat[mat > lfc_clip] <- lfc_clip
  mat[mat < -lfc_clip] <- -lfc_clip

  cat("[INFO] Heatmap: ", length(all_genes), " genes across ",
      length(unique(gene_to_pathway)), " pathways\n", sep = "")

  ## Color scale: blue-white-red, centered at 0
  col_fun <- colorRamp2(
    c(-lfc_clip, 0, lfc_clip),
    c("#2166AC", "white", "#B2182B")  ## Blue-white-red
  )

  ## Create pathway factor for row splitting
  pathway_factor <- factor(gene_to_pathway, levels = pathway_order)

  ## Direction label for title
  direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
  plot_title <- paste0(contrast_name, "\n", direction_label, " Pathways (", database, ")")

  ## Parameter caption for bottom annotation
  param_caption <- paste0(
    "Cutoffs: Pathway FDR < 0.05 | Count \u2265 3 | Top ", top_pathways, " pathways | ",
    "Max ", max_genes, " genes/pathway | log2FC clipped at \u00b1", lfc_clip
  )

  ## Create heatmap
  ht <- Heatmap(
    mat,
    col = col_fun,
    name = "log2FC\n(DESeq2)",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    row_split = pathway_factor,
    row_gap = unit(2, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 7, fontface = "bold", lineheight = 0.9),
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title = plot_title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold", lineheight = 1.1),
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 0.5),
    rect_gp = gpar(col = "grey80", lwd = 0.3),
    width = unit(2.2, "cm"),
    heatmap_legend_param = list(
      title = "log2 fold change\n(DESeq2)",
      title_position = "topcenter",
      title_gp = gpar(fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontsize = 6),
      legend_height = unit(2.5, "cm"),
      at = c(-lfc_clip, 0, lfc_clip),
      labels = c(paste0("≤", -lfc_clip), "0", paste0("≥", lfc_clip))
    )
  )

  ## Save
  db_tag <- gsub(":", "", database)
  plot_file <- paste0("Heatmap_pathways_", db_tag, "_", contrast_name, "_", direction, "_", run_tag, ".png")
  n_genes <- length(all_genes)
  n_pathways <- length(unique(gene_to_pathway))

  png(file.path(output_dir, plot_file),
      width = 650,
      height = max(450, n_genes * 14 + n_pathways * 40 + 120),
      res = 150)
  draw(ht, padding = unit(c(5, 15, 12, 5), "mm"))
  ## Add parameter caption at bottom
  grid.text(
    param_caption,
    x = 0.5, y = unit(3, "mm"),
    just = "center",
    gp = gpar(fontsize = 6, col = "grey40")
  )
  dev.off()

  cat("[OK] Saved: ", plot_file, "\n", sep = "")

  return(ht)
}

## Create pathway heatmaps for all contrasts and directions
cat("\n--- Generating pathway heatmaps ---\n")

## Get unique contrasts from DE files
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
contrasts <- unique(sapply(de_files, function(f) {
  str_match(basename(f), "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
}))
contrasts <- contrasts[!is.na(contrasts)]

## Generate heatmaps for GO:BP (both Up and Down)
for (contrast in contrasts) {
  for (direction in c("Up", "Down")) {
    create_pathway_heatmap(contrast, direction, database = "GO:BP")
  }
}

## Generate heatmaps for KEGG (both Up and Down)
for (contrast in contrasts) {
  for (direction in c("Up", "Down")) {
    create_pathway_heatmap(contrast, direction, database = "KEGG")
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
