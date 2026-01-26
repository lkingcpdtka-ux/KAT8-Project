#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: PUBLICATION QUALITY PLOTS
## =========================================================
## This script creates publication-quality visualizations:
## 1. Dot plots for pathway enrichment (like PDIA3 paper Figure 3E)
## 2. Grouped heatmaps with gene categories (like Figure 3A/3D)
## 3. GSEA enrichment plots
##
## Run AFTER parts 1-3 have completed
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "ComplexHeatmap",
  "circlize",
  "org.Mm.eg.db",
  "AnnotationDbi"
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
## SECTION 1: DOT PLOTS FOR PATHWAY ENRICHMENT
## ============================================================
## Creates dot plots similar to PDIA3 paper Figure 3E
## - Y-axis: Pathway/term names
## - X-axis: Gene ratio (proportion of DEGs in pathway)
## - Dot size: Number of genes (GeneNumber)
## - Dot color: P-value (gradient)
## ============================================================

cat("\n=== CREATING ENRICHMENT DOT PLOTS ===\n")

create_enrichment_dotplot <- function(ora_file, contrast_name, direction,
                                       top_n = 20, output_dir = plots_dir) {

  ## Read ORA results
  if (!file.exists(ora_file)) {
    cat("[WARN] File not found: ", ora_file, "\n")
    return(NULL)
  }

  ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)

  if (nrow(ora_data) == 0) {
    cat("[WARN] No pathways in: ", basename(ora_file), "\n")
    return(NULL)
  }

  ## Parse GeneRatio to numeric if not already
  if ("GeneRatio" %in% colnames(ora_data) && is.character(ora_data$GeneRatio)) {
    ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      num <- as.numeric(x[1])
      denom <- as.numeric(x[2])
      if (is.na(denom) || denom == 0) return(0)
      num / denom
    })
  } else if ("GeneRatio_numeric" %in% colnames(ora_data)) {
    ## Already numeric
  } else {
    cat("[WARN] Cannot find GeneRatio column\n")
    return(NULL)
  }

  ## Calculate gene count from GeneRatio if Count not present
  if (!"Count" %in% colnames(ora_data)) {
    if ("GeneRatio" %in% colnames(ora_data)) {
      ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
        as.numeric(x[1])
      })
    } else {
      ora_data$Count <- 10  ## Default
    }
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

  ## Clean up pathway descriptions
  plot_data <- plot_data %>%
    dplyr::mutate(
      Description = str_wrap(Description, width = 50),
      Description = factor(Description, levels = rev(Description))
    )

  ## Color palette based on direction
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

## Find and process all ORA GO:BP results
ora_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)

for (ora_file in ora_files) {
  filename <- basename(ora_file)
  ## Extract contrast and direction from filename
  ## Pattern: ORA_gobp_[contrast]_[direction]_[timestamp].csv
  parts <- str_match(filename, "ORA_gobp_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
  if (!is.na(parts[1])) {
    contrast <- parts[2]
    direction <- parts[3]
    create_enrichment_dotplot(ora_file, contrast, direction)
  }
}

## ============================================================
## SECTION 2: GROUPED HEATMAPS WITH GENE CATEGORIES
## ============================================================
## Creates heatmaps similar to PDIA3 paper Figure 3A/3D
## - Rows: Genes grouped by biological category
## - Columns: Samples split by condition
## - Row splits show category labels (Inflammatory, ECM, etc.)
## ============================================================

cat("\n=== CREATING GROUPED GENE HEATMAPS ===\n")

## Define gene categories for KAT8 study
## Customize these based on your pathway analysis results!
gene_categories <- list(

  ## Inflammatory response genes (from GO:BP inflammatory response)
  "Inflammatory\nResponse" = c(
    "Il1b", "Il6", "Tnf", "Il1a", "Il18", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
    "Cxcl1", "Cxcl2", "Cxcl10", "Nlrp3", "Nfkb1", "Ptgs2", "Nos2",
    "Cd14", "Tlr2", "Tlr4", "Myd88"
  ),

  ## Chemotaxis / Migration genes
  "Chemotaxis /\nMigration" = c(
    "Ccl2", "Ccl7", "Cxcl12", "Cxcl14", "Ccr2", "Ccr5", "Cxcr4",
    "Cx3cl1", "Cx3cr1", "Itgam", "Itgb2", "Mmp9", "Mmp2"
  ),

  ## ECM / Fibrosis genes
  "ECM /\nFibrosis" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col5a1", "Col5a3",
    "Col6a1", "Col6a2", "Col6a3", "Col6a6", "Col15a1",
    "Fn1", "Lama4", "Lamb1", "Lamc1",
    "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14",
    "Timp1", "Timp2", "Timp3", "Timp4",
    "Tgfb1", "Tgfb2", "Tgfb3", "Acta2", "Ctgf"
  ),

  ## Adipocyte / Metabolic genes
  "Adipocyte /\nMetabolic" = c(
    "Pparg", "Ppara", "Ppargc1a", "Srebf1", "Fabp4", "Adipoq", "Lep",
    "Lpl", "Fasn", "Acaca", "Scd1", "Dgat1", "Dgat2",
    "Glut4", "Slc2a4", "Irs1", "Irs2", "Plin1", "Cidea", "Ucp1"
  ),

  ## Macrophage markers
  "Macrophage\nMarkers" = c(
    "Cd68", "Adgre1", "Itgam", "Mrc1", "Cd163", "Cd86", "Cd80",
    "H2-Aa", "H2-Ab1", "Csf1r", "Arg1", "Nos2", "Cd11c", "Itgax"
  )
)

## Alternative: Simpler categories for KAT8 (ECM-focused)
gene_categories_ecm_focused <- list(

  "Collagens" = c(
    "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
    "Col5a1", "Col5a3", "Col6a1", "Col6a2", "Col6a3",
    "Col6a6", "Col15a1", "Col14a1", "Col8a1"
  ),

  "ECM\nRemodeling" = c(
    "Fn1", "Mmp2", "Mmp3", "Mmp9", "Mmp12", "Mmp14", "Mmp16",
    "Timp1", "Timp2", "Timp3", "Timp4",
    "Adamts1", "Adamts4", "Lox", "Loxl1", "Loxl2"
  ),

  "Laminins /\nBasement\nMembrane" = c(
    "Lama1", "Lama2", "Lama4", "Lama5",
    "Lamb1", "Lamb2", "Lamc1", "Lamc2",
    "Nid1", "Nid2", "Hspg2"
  ),

  "Adipocyte\nMarkers" = c(
    "Lep", "Adipoq", "Pparg", "Ppargc1a", "Ppara", "Srebf1",
    "Fabp4", "Plin1", "Cidea", "Ucp1", "Sorbs1"
  ),

  "Metabolism" = c(
    "Lpl", "Fasn", "Acaca", "Scd1", "Dgat1", "Dgat2",
    "Pck1", "G6pc", "Acly", "Hmgcr"
  )
)

## Function to create grouped heatmap
create_grouped_heatmap <- function(de_file, contrast_name, gene_categories,
                                    vst_file = NULL, output_dir = plots_dir) {

  cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

  ## Load DE results
  if (!file.exists(de_file)) {
    cat("[WARN] DE file not found: ", de_file, "\n")
    return(NULL)
  }

  de_data <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Flatten gene categories and find which genes are in data
  all_category_genes <- unique(unlist(gene_categories))
  genes_in_data <- intersect(all_category_genes, rownames(de_data))

  cat("[INFO] Category genes found in data: ", length(genes_in_data), "/",
      length(all_category_genes), "\n", sep = "")

  if (length(genes_in_data) < 5) {
    cat("[WARN] Too few genes found for heatmap\n")
    return(NULL)
  }

  ## Create gene-to-category mapping (use first category if gene in multiple)
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

  ## Remove genes not assigned to any category
  gene_to_category <- gene_to_category[gene_to_category != ""]
  genes_to_plot <- names(gene_to_category)

  if (length(genes_to_plot) < 5) {
    cat("[WARN] Too few categorized genes for heatmap\n")
    return(NULL)
  }

  ## Get log2FC and significance for heatmap
  ## Use log2FoldChange or logFC column
  if ("log2FoldChange" %in% colnames(de_data)) {
    lfc_col <- "log2FoldChange"
  } else if ("logFC" %in% colnames(de_data)) {
    lfc_col <- "logFC"
  } else {
    cat("[WARN] Cannot find logFC column\n")
    return(NULL)
  }

  if ("padj" %in% colnames(de_data)) {
    fdr_col <- "padj"
  } else if ("adj.P.Val" %in% colnames(de_data)) {
    fdr_col <- "adj.P.Val"
  } else {
    fdr_col <- NULL
  }

  ## Create matrix for heatmap (log2FC values)
  mat <- as.matrix(de_data[genes_to_plot, lfc_col, drop = FALSE])
  colnames(mat) <- "log2FC"

  ## Order genes by category
  gene_order <- names(sort(gene_to_category))
  mat <- mat[gene_order, , drop = FALSE]
  gene_to_category <- gene_to_category[gene_order]

  ## Create significance markers
  sig_markers <- NULL
  if (!is.null(fdr_col)) {
    fdr_vals <- de_data[gene_order, fdr_col]
    sig_markers <- ifelse(fdr_vals < 0.001, "***",
                          ifelse(fdr_vals < 0.01, "**",
                                 ifelse(fdr_vals < 0.05, "*", "")))
  }

  ## Color scale (teal-white-orange like existing heatmaps)
  max_val <- max(abs(mat), na.rm = TRUE)
  if (max_val == 0) max_val <- 1
  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),
    c("#0072B2", "white", "#E69F00")
  )

  ## Row annotation for significance
  if (!is.null(sig_markers)) {
    row_ha <- rowAnnotation(
      Sig = anno_text(sig_markers, gp = gpar(fontsize = 8))
    )
  } else {
    row_ha <- NULL
  }

  ## Create heatmap with row splits by category
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
    column_title = paste0("Grouped Gene Expression: ", contrast_name),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    border = TRUE,
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title = "log2FC",
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      legend_height = unit(4, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ## Add text for significant genes
      if (!is.null(sig_markers)) {
        grid.text(sig_markers[i], x, y, gp = gpar(fontsize = 7))
      }
    }
  )

  ## Save
  plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
  n_genes <- length(genes_to_plot)
  n_categories <- length(unique(gene_to_category))

  png(file.path(output_dir, plot_file),
      width = 1200,
      height = max(800, n_genes * 25 + n_categories * 50),
      res = 150)
  draw(ht)
  dev.off()

  cat("[OK] Saved grouped heatmap: ", plot_file, "\n", sep = "")

  return(ht)
}

## Find DE tables and create grouped heatmaps
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

for (de_file in de_files) {
  filename <- basename(de_file)
  contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
  if (!is.na(contrast)) {
    ## Create heatmap with ECM-focused categories (more relevant for KAT8)
    create_grouped_heatmap(de_file, contrast, gene_categories_ecm_focused)
  }
}

## ============================================================
## SECTION 3: COMBINED DOT PLOT (Multiple pathways in one)
## ============================================================
## Creates a combined dot plot showing top pathways from
## multiple databases or contrasts for comparison
## ============================================================

cat("\n=== CREATING COMBINED ENRICHMENT DOT PLOT ===\n")

create_combined_dotplot <- function(tables_dir, contrasts = NULL, databases = c("gobp"),
                                     direction = "Down", top_n_per = 10, output_dir = plots_dir) {

  all_data <- data.frame()

  for (db in databases) {
    pattern <- paste0("^ORA_", db, "_.*_", direction, "_.*\\.csv$")
    files <- list.files(tables_dir, pattern = pattern, full.names = TRUE)

    for (f in files) {
      filename <- basename(f)
      parts <- str_match(filename, paste0("ORA_", db, "_(.+)_", direction, "_[0-9]+_[0-9]+\\.csv"))

      if (!is.na(parts[1])) {
        contrast <- parts[2]

        if (!is.null(contrasts) && !(contrast %in% contrasts)) next

        data <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
        if (is.null(data) || nrow(data) == 0) next

        ## Parse GeneRatio
        if ("GeneRatio" %in% colnames(data) && is.character(data$GeneRatio)) {
          data$GeneRatio_numeric <- sapply(strsplit(data$GeneRatio, "/"), function(x) {
            as.numeric(x[1]) / as.numeric(x[2])
          })
        }

        if (!"Count" %in% colnames(data)) {
          data$Count <- sapply(strsplit(data$GeneRatio, "/"), function(x) as.numeric(x[1]))
        }

        ## Add metadata
        data$Contrast <- contrast
        data$Database <- toupper(db)

        ## Select top pathways
        top_data <- data %>%
          dplyr::filter(p.adjust < 0.05) %>%
          dplyr::arrange(p.adjust) %>%
          dplyr::slice_head(n = top_n_per)

        all_data <- rbind(all_data, top_data)
      }
    }
  }

  if (nrow(all_data) == 0) {
    cat("[WARN] No data found for combined dot plot\n")
    return(NULL)
  }

  ## Clean up descriptions and create unique labels
  all_data <- all_data %>%
    dplyr::mutate(
      Description_short = str_trunc(Description, 60),
      Label = paste0(Description_short, " (", Contrast, ")")
    ) %>%
    dplyr::arrange(Contrast, p.adjust) %>%
    dplyr::mutate(Label = factor(Label, levels = rev(unique(Label))))

  ## Color scale
  if (direction == "Down") {
    color_low <- "#DEEBF7"
    color_high <- "#08519C"
  } else {
    color_low <- "#FEE0D2"
    color_high <- "#CB181D"
  }

  ## Create plot
  p <- ggplot(all_data, aes(x = GeneRatio_numeric, y = Label)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(
      low = color_high,
      high = color_low,
      name = "Adj.\nP-value",
      guide = guide_colorbar(reverse = TRUE)
    ) +
    scale_size_continuous(
      name = "Gene\nNumber",
      range = c(4, 12)
    ) +
    facet_grid(Contrast ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = paste0("GO Enrichment Comparison (", direction, "-regulated)"),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      strip.text.y = element_text(angle = 0, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )

  ## Save
  plot_file <- paste0("Dotplot_combined_", direction, "_", run_tag, ".png")
  n_terms <- length(unique(all_data$Label))

  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 12,
    height = max(8, n_terms * 0.3),
    dpi = 300,
    bg = "white"
  )

  cat("[OK] Saved combined dot plot: ", plot_file, "\n", sep = "")

  return(p)
}

## Create combined plots for up and down-regulated genes
create_combined_dotplot(tables_dir, direction = "Down")
create_combined_dotplot(tables_dir, direction = "Up")

## ============================================================
## SECTION 4: FGSEA-STYLE ENRICHMENT DOT PLOT
## ============================================================
## Creates a dot plot from FGSEA results showing NES
## ============================================================

cat("\n=== CREATING FGSEA DOT PLOTS ===\n")

create_fgsea_dotplot <- function(tables_dir, database = "hallmark",
                                  contrast = NULL, top_n = 15, output_dir = plots_dir) {

  pattern <- paste0("^fgsea_", database, "_.*\\.csv$")
  files <- list.files(tables_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    cat("[WARN] No FGSEA files found for: ", database, "\n")
    return(NULL)
  }

  for (f in files) {
    filename <- basename(f)
    parts <- str_match(filename, paste0("fgsea_", database, "_(.+)_[0-9]+_[0-9]+\\.csv"))

    if (is.na(parts[1])) next

    contrast_name <- parts[2]
    if (!is.null(contrast) && contrast_name != contrast) next

    data <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(data) || nrow(data) == 0) {
      cat("[WARN] No data in: ", basename(f), "\n")
      next
    }

    ## Get top pathways by absolute NES
    top_up <- data %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = top_n/2)
    top_down <- data %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = top_n/2)

    plot_data <- rbind(top_down, top_up) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        Description_clean = ifelse(is.na(Description) | Description == "", pathway, Description),
        Description_clean = str_trunc(Description_clean, 50),
        Description_clean = factor(Description_clean, levels = Description_clean),
        Direction = ifelse(NES > 0, "Enriched (Up)", "Enriched (Down)")
      )

    if (nrow(plot_data) == 0) next

    ## Create dot plot with NES
    p <- ggplot(plot_data, aes(x = NES, y = Description_clean)) +
      geom_point(aes(size = size, color = padj)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_gradient(
        low = "#08519C",
        high = "#DEEBF7",
        name = "Adj.\nP-value",
        guide = guide_colorbar(reverse = TRUE)
      ) +
      scale_size_continuous(
        name = "Gene Set\nSize",
        range = c(4, 12)
      ) +
      labs(
        title = paste0("GSEA: ", toupper(database), " (", contrast_name, ")"),
        subtitle = "Normalized Enrichment Score (NES)",
        x = "NES",
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major.y = element_line(color = "grey90", linetype = "dotted")
      )

    ## Save
    plot_file <- paste0("Dotplot_fgsea_", database, "_", contrast_name, "_", run_tag, ".png")
    ggsave(
      file.path(output_dir, plot_file),
      plot = p,
      width = 10,
      height = max(6, nrow(plot_data) * 0.35),
      dpi = 300,
      bg = "white"
    )

    cat("[OK] Saved FGSEA dot plot: ", plot_file, "\n", sep = "")
  }
}

## Create FGSEA dot plots for each database
for (db in c("hallmark", "kegg", "gobp", "wikipathways")) {
  create_fgsea_dotplot(tables_dir, database = db)
}

cat("\n======================================\n")
cat("=== PART 4 (Publication Plots) COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - Enrichment dot plots (ORA results)\n")
cat("  - Grouped gene heatmaps (by category)\n")
cat("  - Combined comparison dot plots\n")
cat("  - FGSEA dot plots\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
