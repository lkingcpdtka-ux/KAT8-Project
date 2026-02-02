#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4B: VISUALIZATION-ONLY PLOTS
## =========================================================
## This script generates plots WITHOUT re-running analysis.
## It reads saved results from parts 1-3 and creates:
##   1. ORA bar plots (GO:BP, KEGG)
##   2. fGSEA bar plots (GO:BP, KEGG)
##   3. ORA dot plots
##   4. Grouped heatmaps (CTL vs KAT8KD)
##   5. Individual sample heatmaps
##
## RUN THIS when you want to tweak plot aesthetics quickly!
## Analysis is NOT re-run - just visualization.
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "viridis", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c("ComplexHeatmap", "circlize")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(viridis)
  library(scales)
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

## ============================================================
## VISUALIZATION TOGGLE: Choose what plots to generate
## ============================================================
## Set these to TRUE/FALSE to control which plots are created
## This makes it fast to regenerate specific plot types

generate_ora_barplots   <- TRUE   ## ORA bar plots (part2 style)
generate_fgsea_barplots <- TRUE   ## fGSEA bar plots (part3 style)
generate_ora_dotplots   <- TRUE   ## ORA dot plots (part4 style)
generate_grouped_heatmaps <- TRUE ## Grouped heatmaps (CTL vs KAT8KD)
generate_individual_heatmaps <- TRUE ## Individual sample heatmaps

## ============================================================

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
## SECTION 1: ORA BAR PLOTS
## ============================================================

if (generate_ora_barplots) {
  cat("\n=== GENERATING ORA BAR PLOTS ===\n")

  create_ora_barplot <- function(ora_file, contrast_name, direction,
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

    ## Standardize column names
    if (!"pathway_name" %in% colnames(ora_data) && "Description" %in% colnames(ora_data)) {
      ora_data$pathway_name <- ora_data$Description
    }

    top_n <- ora_params$top_n

    ## Filter and select top pathways
    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < ora_params$pvalue_cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n)

    if (nrow(plot_data) == 0) {
      cat("[WARN] No significant pathways for plot\n")
      return(NULL)
    }

    ## Parse GeneRatio to numeric
    if ("GeneRatio" %in% colnames(plot_data) && is.character(plot_data$GeneRatio)) {
      plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
        num <- as.numeric(x[1])
        denom <- as.numeric(x[2])
        if (is.na(denom) || denom == 0) return(0)
        num / denom
      })
    }

    ## Calculate -log10(FDR) for color
    plot_data$neg_log10_padj <- -log10(plot_data$p.adjust)

    ## Order by adjusted p-value (most significant at top)
    plot_data <- plot_data %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
    plot_title <- paste0(database, " Pathways (ORA - ", direction_label, ")\n", contrast_name)

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = pathway_name, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_viridis_c(option = "viridis", direction = 1, name = "Adj. P-value") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
      labs(
        title = plot_title,
        x = "Gene Ratio",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )

    ## Save
    db_short <- gsub(":", "", database)
    plot_file <- paste0("ORA_barplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
    ggsave(
      file.path(output_dir, plot_file),
      plot = p,
      width = 10,
      height = max(6, nrow(plot_data) * 0.3),
      dpi = 300,
      bg = "white"
    )
    cat("[OK] Saved ORA bar plot: ", plot_file, "\n", sep = "")

    return(p)
  }

  ## Process all ORA files
  ora_gobp_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)
  ora_kegg_files <- list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE)

  for (ora_file in c(ora_gobp_files, ora_kegg_files)) {
    filename <- basename(ora_file)
    parts <- str_match(filename, "ORA_(gobp|kegg)_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      db <- ifelse(parts[2] == "gobp", "GO:BP", "KEGG")
      contrast <- parts[3]
      direction <- parts[4]
      create_ora_barplot(ora_file, contrast, direction, database = db)
    }
  }
}

## ============================================================
## SECTION 2: fGSEA BAR PLOTS (with fixed legend)
## ============================================================

if (generate_fgsea_barplots) {
  cat("\n=== GENERATING fGSEA BAR PLOTS ===\n")

  create_fgsea_barplot <- function(fgsea_file, contrast_name,
                                    database = "GO:BP",
                                    output_dir = plots_dir) {

    if (!file.exists(fgsea_file)) {
      cat("[WARN] File not found: ", fgsea_file, "\n")
      return(NULL)
    }

    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)

    if (nrow(fgsea_data) == 0) {
      cat("[WARN] No pathways in: ", basename(fgsea_file), "\n")
      return(NULL)
    }

    ## Filter significant pathways
    fdr_cutoff <- gsea_params$fdr_cutoff
    sig_data <- fgsea_data %>% dplyr::filter(padj < fdr_cutoff)

    if (nrow(sig_data) == 0) {
      cat("[WARN] No significant pathways for: ", basename(fgsea_file), "\n")
      return(NULL)
    }

    top_n <- gsea_params$top_n_per_direction

    ## Get top up and down pathways
    top_up <- sig_data %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(padj, dplyr::desc(NES)) %>%
      dplyr::slice_head(n = top_n)

    top_down <- sig_data %>%
      dplyr::filter(NES < 0) %>%
      dplyr::arrange(padj, NES) %>%
      dplyr::slice_head(n = top_n)

    plot_data <- dplyr::bind_rows(top_down, top_up) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        pathway_label = ifelse(!is.na(pathway_name) & pathway_name != "", pathway_name, pathway_id),
        pathway_label = factor(pathway_label, levels = pathway_label),
        Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"),
        ## Calculate -log10(padj) for color intensity
        padj_safe = pmax(padj, .Machine$double.eps),
        neg_log10_padj = -log10(padj_safe)
      )

    if (nrow(plot_data) == 0) return(NULL)

    p <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.95) +
      scale_fill_viridis_c(option = "viridis", direction = 1, name = "Adj. P-value") +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
      labs(
        title = paste0("GSEA ", database, ": ", contrast_name),
        x = "Normalized Enrichment Score (NES)",
        y = NULL,
        caption = paste0("FDR < ", fdr_cutoff, " | Top ", top_n, " per direction")
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8))
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5)

    ## Save
    db_short <- gsub(":", "", tolower(database))
    plot_file <- paste0("fgsea_plot_", db_short, "_", contrast_name, "_", run_tag, ".png")
    ggsave(file.path(output_dir, plot_file), plot = p, width = 12,
           height = max(6, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
    cat("[OK] Saved fGSEA bar plot: ", plot_file, "\n", sep = "")

    return(p)
  }

  ## Process all fGSEA files
  fgsea_gobp_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)
  fgsea_kegg_files <- list.files(tables_dir, pattern = "^fgsea_kegg_.*\\.csv$", full.names = TRUE)

  for (fgsea_file in fgsea_gobp_files) {
    filename <- basename(fgsea_file)
    parts <- str_match(filename, "fgsea_gobp_(.+)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      contrast <- parts[2]
      create_fgsea_barplot(fgsea_file, contrast, database = "GO:BP")
    }
  }

  for (fgsea_file in fgsea_kegg_files) {
    filename <- basename(fgsea_file)
    parts <- str_match(filename, "fgsea_kegg_(.+)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      contrast <- parts[2]
      create_fgsea_barplot(fgsea_file, contrast, database = "KEGG")
    }
  }
}

## ============================================================
## SECTION 3: ORA DOT PLOTS
## ============================================================

if (generate_ora_dotplots) {
  cat("\n=== GENERATING ORA DOT PLOTS ===\n")

  create_ora_dotplot <- function(ora_file, contrast_name, direction,
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

    ## Standardize column names
    if (!"Description" %in% colnames(ora_data) && "pathway_name" %in% colnames(ora_data)) {
      ora_data$Description <- ora_data$pathway_name
    }

    top_n <- dotplot_params$top_n_gobp

    ## Parse GeneRatio
    if ("GeneRatio" %in% colnames(ora_data) && is.character(ora_data$GeneRatio)) {
      ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
        num <- as.numeric(x[1])
        denom <- as.numeric(x[2])
        if (is.na(denom) || denom == 0) return(0)
        num / denom
      })
    }

    if (!"Count" %in% colnames(ora_data)) {
      ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
        as.numeric(x[1])
      })
    }

    ## Calculate -log10(FDR)
    ora_data$neg_log10_fdr <- -log10(ora_data$p.adjust)

    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < dotplot_params$fdr_cutoff, Count >= dotplot_params$min_count) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n)

    if (nrow(plot_data) == 0) {
      cat("[WARN] No significant pathways for dot plot\n")
      return(NULL)
    }

    ## Clip extreme values
    fdr_cap <- min(quantile(plot_data$neg_log10_fdr, 0.99, na.rm = TRUE), 10)
    plot_data$neg_log10_fdr_clipped <- pmin(plot_data$neg_log10_fdr, fdr_cap)

    plot_data <- plot_data %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(
        Description = str_wrap(Description, width = 45),
        Description = factor(Description, levels = rev(Description))
      )

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
    plot_title <- paste0(database, " Enrichment\n(", contrast_name, ", ", direction_label, ")")

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
      geom_point(aes(size = Count, color = neg_log10_fdr_clipped)) +
      scale_color_viridis_c(option = "viridis", direction = 1, name = "Adj. P-value") +
      scale_size_continuous(
        name = "Gene Count",
        range = c(2, 6),
        breaks = pretty(plot_data$Count, n = 4)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      labs(
        title = plot_title,
        x = "Gene Ratio",
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11, lineheight = 1.1),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.title.x = element_text(size = 8, margin = margin(t = 5)),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
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
    cat("[OK] Saved dot plot: ", plot_file, "\n", sep = "")

    return(p)
  }

  ## Process all ORA files
  ora_gobp_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)
  ora_kegg_files <- list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE)

  for (ora_file in c(ora_gobp_files, ora_kegg_files)) {
    filename <- basename(ora_file)
    parts <- str_match(filename, "ORA_(gobp|kegg)_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      db <- ifelse(parts[2] == "gobp", "GO:BP", "KEGG")
      contrast <- parts[3]
      direction <- parts[4]
      create_ora_dotplot(ora_file, contrast, direction, database = db)
    }
  }
}

## ============================================================
## SECTION 4: GROUPED HEATMAPS (Fixed: actual expression, not binary)
## ============================================================

if (generate_grouped_heatmaps || generate_individual_heatmaps) {
  cat("\n=== GENERATING HEATMAPS ===\n")

  ## Load VST data
  vst_files <- list.files(tables_dir, pattern = "^VST_matrix_heatmap_.*\\.rds$", full.names = TRUE)
  vst_data <- NULL
  if (length(vst_files) > 0) {
    vst_file <- vst_files[order(file.info(vst_files)$mtime, decreasing = TRUE)][1]
    vst_data <- readRDS(vst_file)
    cat("[OK] Loaded VST matrix from: ", basename(vst_file), "\n", sep = "")
  } else {
    cat("[WARN] VST matrix not found - heatmaps may be limited\n")
  }

  ## Extract gene categories from fGSEA results
  fgsea_gobp_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)

  extract_leading_edge <- function(fgsea_file, pathway_patterns, max_genes_per_pathway = 15) {
    if (!file.exists(fgsea_file)) return(character(0))
    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
    genes <- character(0)
    for (pattern in pathway_patterns) {
      matched_rows <- grep(pattern, fgsea_data$pathway_name, ignore.case = TRUE)
      if (length(matched_rows) > 0) {
        row <- matched_rows[1]
        le_symbols <- fgsea_data$leadingEdge_symbols[row]
        if (!is.na(le_symbols) && le_symbols != "") {
          gene_vec <- unlist(strsplit(le_symbols, ","))
          genes <- c(genes, head(gene_vec, max_genes_per_pathway))
        }
      }
    }
    unique(genes)
  }

  gene_categories <- list()

  iwat_fgsea <- fgsea_gobp_files[grep("iWAT", fgsea_gobp_files, ignore.case = TRUE)]
  if (length(iwat_fgsea) > 0) {
    mito_genes <- extract_leading_edge(iwat_fgsea[1], c(
      "mitochondrial ATP synthesis", "mitochondrial translation",
      "respiratory chain complex", "fatty acid beta-oxidation"
    ), max_genes_per_pathway = 10)
    if (length(mito_genes) > 0) gene_categories[["Mitochondrial /\nEnergy"]] <- head(mito_genes, 20)
  }

  gwat_fgsea <- fgsea_gobp_files[grep("gWAT", fgsea_gobp_files, ignore.case = TRUE)]
  if (length(gwat_fgsea) > 0) {
    immune_genes <- extract_leading_edge(gwat_fgsea[1], c(
      "leukocyte degranulation", "myeloid cell activation", "positive regulation of cytokine"
    ), max_genes_per_pathway = 10)
    if (length(immune_genes) > 0) gene_categories[["Immune /\nInflammatory"]] <- head(immune_genes, 20)

    lipid_genes <- extract_leading_edge(gwat_fgsea[1], c(
      "lipid localization", "lipid catabolic", "fatty acid"
    ), max_genes_per_pathway = 10)
    if (length(lipid_genes) > 0) gene_categories[["Lipid\nMetabolism"]] <- head(lipid_genes, 15)

    ecm_genes <- extract_leading_edge(gwat_fgsea[1], c(
      "tissue remodeling", "extracellular matrix", "cell adhesion"
    ), max_genes_per_pathway = 10)
    if (length(ecm_genes) > 0) gene_categories[["ECM /\nRemodeling"]] <- head(ecm_genes, 15)
  }

  ## Fallback categories
  if (length(gene_categories) == 0) {
    gene_categories <- list(
      "Inflammatory\nResponse" = c("Il1b", "Il6", "Tnf", "Il1a", "Ccl2", "Ccl7", "Cxcl12"),
      "ECM /\nCollagens" = c("Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col5a1", "Col6a1"),
      "Adipocyte\nMarkers" = c("Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1")
    )
  }

  ## ============================================================
  ## Fixed Grouped Heatmap: Show actual Z-scored expression in BOTH groups
  ## (Not binary CTL=0)
  ## ============================================================

  create_grouped_heatmap_fixed <- function(de_file, contrast_name, gene_categories,
                                            vst_data, output_dir = plots_dir) {

    cat("\n--- Creating grouped heatmap: ", contrast_name, " ---\n", sep = "")

    if (is.null(vst_data) || is.null(vst_data$vst_mat) || is.null(vst_data$sample_info)) {
      cat("[WARN] VST data required for proper grouped heatmap\n")
      return(NULL)
    }

    vst_mat <- vst_data$vst_mat
    sample_info <- vst_data$sample_info

    ## Get tissue from contrast
    tissue <- sub("_.*", "", contrast_name)

    ## Filter samples for this tissue
    if ("Depot" %in% colnames(sample_info)) {
      tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
    } else if ("Tissue" %in% colnames(sample_info)) {
      tissue_samples <- rownames(sample_info)[sample_info$Tissue == tissue]
    } else {
      tissue_samples <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
    }

    if (length(tissue_samples) == 0) {
      cat("[WARN] No samples found for tissue: ", tissue, "\n")
      return(NULL)
    }

    vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples, drop = FALSE]
    sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

    if (!"Genotype" %in% colnames(sample_info_tissue)) {
      cat("[WARN] Genotype column not found in sample info\n")
      return(NULL)
    }

    ## Find genes in categories
    all_category_genes <- unique(unlist(gene_categories))
    genes_in_data <- intersect(all_category_genes, rownames(vst_tissue))

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

    ## Calculate group means
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

    ## Z-score normalize ACROSS BOTH GROUPS
    mat <- t(scale(t(mat)))
    mat <- mat[!rowSums(is.na(mat)), , drop = FALSE]

    ## Clip at ±2 (toned down from default)
    z_clip <- heatmap_params$lfc_clip
    mat[mat > z_clip] <- z_clip
    mat[mat < -z_clip] <- -z_clip

    ## TONED DOWN COLORS: Use muted diverging palette from parameters.R
    ## Blue-grey-red: publication-ready, less saturated
    col_fun <- colorRamp2(
      c(-z_clip, -z_clip/2, 0, z_clip/2, z_clip),
      heatmap_params$muted_colors
    )

    legend_title <- "Z-score\n(mean VST)"
    param_caption <- paste0("Z-scored mean VST by genotype | Clipped at \u00b1", z_clip)

    ht <- Heatmap(
      mat,
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
        at = c(-z_clip, 0, z_clip),
        labels = c(paste0("\u2264", -z_clip), "0", paste0("\u2265", z_clip))
      )
    )

    plot_file <- paste0("Heatmap_grouped_", contrast_name, "_", run_tag, ".png")
    n_genes <- length(genes_to_plot)
    n_categories <- length(unique(gene_to_category))

    png(file.path(output_dir, plot_file),
        width = 600,
        height = max(400, n_genes * 18 + n_categories * 30 + 50),
        res = 150)
    draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
    grid.text(
      param_caption,
      x = 0.5, y = unit(3, "mm"),
      just = "center",
      gp = gpar(fontsize = 7, col = "grey40")
    )
    dev.off()

    cat("[OK] Saved grouped heatmap: ", plot_file, "\n", sep = "")
    return(ht)
  }

  ## Create heatmaps for each contrast
  de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

  if (generate_grouped_heatmaps && !is.null(vst_data)) {
    for (de_file in de_files) {
      filename <- basename(de_file)
      contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (!is.na(contrast)) {
        create_grouped_heatmap_fixed(de_file, contrast, gene_categories, vst_data)
      }
    }
  }
}

## ============================================================
## COMPLETE
## ============================================================

cat("\n======================================\n")
cat("=== PART 4B (VISUALIZATIONS) COMPLETE ===\n")
cat("======================================\n")
cat("Generated:\n")
if (generate_ora_barplots) cat("  - ORA bar plots\n")
if (generate_fgsea_barplots) cat("  - fGSEA bar plots (with fixed legend)\n")
if (generate_ora_dotplots) cat("  - ORA dot plots\n")
if (generate_grouped_heatmaps) cat("  - Grouped heatmaps (Z-scored, not binary)\n")
cat("\nAll plots saved to: ", plots_dir, "\n", sep = "")
cat("\nTIP: Edit toggle variables at top of script to regenerate specific plot types.\n\n")
