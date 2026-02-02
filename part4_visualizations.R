#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: ALL VISUALIZATIONS
## =========================================================
## This script generates ALL plots WITHOUT re-running analysis.
## It reads saved results from parts 1-3 and creates:
##
##   SECTION A: QC PLOTS
##     - PCA plot
##     - MDS plot
##     - Density plots (before/after filtering)
##
##   SECTION B: ENRICHMENT PLOTS
##     - ORA bar plots (GO:BP, KEGG)
##     - fGSEA bar plots (GO:BP, KEGG)
##     - ORA dot plots
##
##   SECTION C: EXPRESSION HEATMAPS
##     - Individual sample heatmaps (z-scored)
##     - Top DEG heatmaps
##     - Genes of Interest heatmaps
##
## RUN THIS when you want to tweak plot aesthetics quickly!
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "viridis", "scales", "ggrepel")
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
  library(ggrepel)
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
## VISUALIZATION TOGGLES
## ============================================================

## Section A: QC Plots
generate_pca_plot       <- TRUE
generate_mds_plot       <- TRUE
generate_density_plot   <- TRUE

## Section B: Enrichment Plots
generate_ora_barplots   <- TRUE
generate_fgsea_barplots <- TRUE
generate_ora_dotplots   <- TRUE

## Section C: Heatmaps
generate_individual_heatmaps <- TRUE
generate_deg_heatmaps        <- TRUE
generate_goi_heatmaps        <- TRUE

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
## SECTION A: QC PLOTS (PCA, MDS, Density)
## ============================================================

if (generate_pca_plot || generate_mds_plot || generate_density_plot) {
  cat("\n=== SECTION A: QC PLOTS ===\n")

  ## Load QC data
  qc_files <- list.files(tables_dir, pattern = "^QC_data_.*\\.rds$", full.names = TRUE)
  if (length(qc_files) == 0) {
    cat("[WARN] QC data not found. Skipping QC plots.\n")
  } else {
    qc_file <- qc_files[order(file.info(qc_files)$mtime, decreasing = TRUE)][1]
    qc_data <- readRDS(qc_file)
    cat("[OK] Loaded QC data from: ", basename(qc_file), "\n", sep = "")

    vst_mat <- qc_data$vst_mat_qc
    sample_info <- qc_data$sample_info
    depot_sex_fill <- qc_data$depot_sex_fill

    ## A.1) PCA Plot
    if (generate_pca_plot && !is.null(vst_mat)) {
      cat("\n--- Creating PCA plot ---\n")

      pca_result <- prcomp(t(vst_mat), scale. = TRUE)
      pca_var <- summary(pca_result)$importance[2, 1:2] * 100

      pca_df <- data.frame(
        Sample   = colnames(vst_mat),
        Genotype = sample_info$Genotype,
        DepotSex = sample_info$DepotSex,
        PC1      = pca_result$x[, 1],
        PC2      = pca_result$x[, 2],
        stringsAsFactors = FALSE
      )

      p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = DepotSex, fill = DepotSex,
                                   shape = Genotype, label = Sample)) +
        geom_point(size = 3, stroke = 1) +
        geom_text_repel(size = 3, max.overlaps = 60, color = "black") +
        scale_color_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_fill_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
        labs(
          title = "PCA",
          x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
          y = paste0("PC2 (", round(pca_var[2], 1), "%)")
        ) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right") +
        guides(
          color = "none",
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "white", color = "black"))
        )

      pca_file <- paste0("PCA_tissue_VST_", run_tag, ".png")
      ggsave(file.path(plots_dir, pca_file), plot = p_pca, width = 8, height = 6, dpi = 300)
      cat("[OK] Saved: ", pca_file, "\n", sep = "")
    }

    ## A.2) MDS Plot
    if (generate_mds_plot && !is.null(vst_mat)) {
      cat("\n--- Creating MDS plot ---\n")

      dist_mat <- dist(t(vst_mat))
      mds_coords <- cmdscale(dist_mat, k = 2)

      mds_df <- data.frame(
        Sample   = colnames(vst_mat),
        Genotype = sample_info$Genotype,
        DepotSex = sample_info$DepotSex,
        MDS1     = mds_coords[, 1],
        MDS2     = mds_coords[, 2],
        stringsAsFactors = FALSE
      )

      p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = DepotSex, fill = DepotSex,
                                   shape = Genotype, label = Sample)) +
        geom_point(size = 3, stroke = 1) +
        geom_text_repel(size = 3, max.overlaps = 60, color = "black") +
        scale_color_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_fill_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
        labs(title = "MDS", x = "MDS1", y = "MDS2") +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right") +
        guides(
          color = "none",
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "white", color = "black"))
        )

      mds_file <- paste0("MDS_tissue_VST_", run_tag, ".png")
      ggsave(file.path(plots_dir, mds_file), plot = p_mds, width = 8, height = 6, dpi = 300)
      cat("[OK] Saved: ", mds_file, "\n", sep = "")
    }

    ## A.3) Density Plot
    if (generate_density_plot && !is.null(qc_data$vst_mat_before)) {
      cat("\n--- Creating density plot ---\n")

      vst_before <- qc_data$vst_mat_before
      vst_after <- vst_mat

      density_file <- paste0("VST_density_before_after_", run_tag, ".png")
      png(file.path(plots_dir, density_file), width = 1600, height = 900, res = 150)

      par(mfrow = c(1, 2), mar = c(5, 5, 5, 2), oma = c(0, 0, 2, 0))

      dens_before <- apply(vst_before, 2, density)
      dens_after  <- apply(vst_after, 2, density)

      xlim_global <- range(c(
        range(sapply(dens_before, function(d) range(d$x))),
        range(sapply(dens_after, function(d) range(d$x)))
      ))

      max_y_before <- max(sapply(dens_before, function(d) max(d$y)))
      max_y_after  <- max(sapply(dens_after, function(d) max(d$y)))

      sample_colors <- depot_sex_fill[sample_info$DepotSex]

      plot(dens_before[[1]],
           main = "Before Filtering",
           sub  = paste0("n=", nrow(vst_before), " genes"),
           xlab = "VST", lwd = 2, col = sample_colors[1],
           xlim = xlim_global, ylim = c(0, max_y_before * 1.15))
      for (i in 2:ncol(vst_before)) {
        lines(dens_before[[i]], lwd = 2, col = sample_colors[i])
      }

      plot(dens_after[[1]],
           main = "After Filtering",
           sub  = paste0("n=", nrow(vst_after), " genes"),
           xlab = "VST", lwd = 2, col = sample_colors[1],
           xlim = xlim_global, ylim = c(0, max_y_after * 1.15))
      for (i in 2:ncol(vst_after)) {
        lines(dens_after[[i]], lwd = 2, col = sample_colors[i])
      }

      dev.off()
      cat("[OK] Saved: ", density_file, "\n", sep = "")
    }
  }
}

## ============================================================
## SECTION B: ENRICHMENT PLOTS
## ============================================================

## B.1) ORA BAR PLOTS
if (generate_ora_barplots) {
  cat("\n=== GENERATING ORA BAR PLOTS ===\n")

  create_ora_barplot <- function(ora_file, contrast_name, direction,
                                  database = "GO:BP", output_dir = plots_dir) {

    if (!file.exists(ora_file)) return(NULL)

    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) return(NULL)

    if (!"pathway_name" %in% colnames(ora_data) && "Description" %in% colnames(ora_data)) {
      ora_data$pathway_name <- ora_data$Description
    }

    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < ora_params$pvalue_cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = ora_params$top_n)

    if (nrow(plot_data) == 0) return(NULL)

    if ("GeneRatio" %in% colnames(plot_data) && is.character(plot_data$GeneRatio)) {
      plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      })
    }

    plot_data$neg_log10_padj <- -log10(plot_data$p.adjust)

    plot_data <- plot_data %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = pathway_name, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_gradientn(
        colors = enrichment_colors$gradient,
        name = expression(-log[10]*"(FDR)")
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
      labs(
        title = paste0(database, " - ", direction_label, "\n", contrast_name),
        x = "Gene Ratio",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      )

    db_short <- gsub(":", "", database)
    plot_file <- paste0("ORA_barplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
    ggsave(file.path(output_dir, plot_file), plot = p,
           width = 10, height = max(5, nrow(plot_data) * 0.3), dpi = 300, bg = "white")
    cat("[OK] Saved: ", plot_file, "\n", sep = "")
    return(p)
  }

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

## B.2) fGSEA BAR PLOTS
if (generate_fgsea_barplots) {
  cat("\n=== GENERATING fGSEA BAR PLOTS ===\n")

  create_fgsea_barplot <- function(fgsea_file, contrast_name,
                                    database = "GO:BP", output_dir = plots_dir) {

    if (!file.exists(fgsea_file)) return(NULL)

    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
    if (nrow(fgsea_data) == 0) return(NULL)

    fdr_cutoff <- gsea_params$fdr_cutoff
    sig_data <- fgsea_data %>% dplyr::filter(padj < fdr_cutoff)
    if (nrow(sig_data) == 0) return(NULL)

    top_n <- gsea_params$top_n_per_direction

    top_up <- sig_data %>%
      dplyr::filter(NES > 0) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = top_n)

    top_down <- sig_data %>%
      dplyr::filter(NES < 0) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = top_n)

    plot_data <- dplyr::bind_rows(top_down, top_up) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        pathway_label = ifelse(!is.na(pathway_name) & pathway_name != "", pathway_name, pathway_id),
        pathway_label = factor(pathway_label, levels = pathway_label),
        neg_log10_padj = -log10(pmax(padj, .Machine$double.eps))
      )

    if (nrow(plot_data) == 0) return(NULL)

    p <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_gradientn(
        colors = enrichment_colors$gradient,
        name = expression(-log[10]*"(FDR)")
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
      labs(
        title = paste0("GSEA ", database, "\n", contrast_name),
        x = "NES",
        y = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 11),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

    db_short <- gsub(":", "", tolower(database))
    plot_file <- paste0("fgsea_plot_", db_short, "_", contrast_name, "_", run_tag, ".png")
    ggsave(file.path(output_dir, plot_file), plot = p,
           width = 11, height = max(5, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
    cat("[OK] Saved: ", plot_file, "\n", sep = "")
    return(p)
  }

  fgsea_gobp_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)
  fgsea_kegg_files <- list.files(tables_dir, pattern = "^fgsea_kegg_.*\\.csv$", full.names = TRUE)

  for (fgsea_file in fgsea_gobp_files) {
    filename <- basename(fgsea_file)
    parts <- str_match(filename, "fgsea_gobp_(.+)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      create_fgsea_barplot(fgsea_file, parts[2], database = "GO:BP")
    }
  }

  for (fgsea_file in fgsea_kegg_files) {
    filename <- basename(fgsea_file)
    parts <- str_match(filename, "fgsea_kegg_(.+)_[0-9]+_[0-9]+\\.csv")
    if (!is.na(parts[1])) {
      create_fgsea_barplot(fgsea_file, parts[2], database = "KEGG")
    }
  }
}

## B.3) ORA DOT PLOTS
if (generate_ora_dotplots) {
  cat("\n=== GENERATING ORA DOT PLOTS ===\n")

  create_ora_dotplot <- function(ora_file, contrast_name, direction,
                                  database = "GO:BP", output_dir = plots_dir) {

    if (!file.exists(ora_file)) return(NULL)

    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) return(NULL)

    if (!"Description" %in% colnames(ora_data) && "pathway_name" %in% colnames(ora_data)) {
      ora_data$Description <- ora_data$pathway_name
    }

    if ("GeneRatio" %in% colnames(ora_data) && is.character(ora_data$GeneRatio)) {
      ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      })
    }

    if (!"Count" %in% colnames(ora_data)) {
      ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
        as.numeric(x[1])
      })
    }

    ora_data$neg_log10_fdr <- -log10(ora_data$p.adjust)

    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < dotplot_params$fdr_cutoff, Count >= dotplot_params$min_count) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = dotplot_params$top_n_gobp)

    if (nrow(plot_data) == 0) return(NULL)

    fdr_cap <- min(quantile(plot_data$neg_log10_fdr, 0.99, na.rm = TRUE), 10)
    plot_data$neg_log10_fdr_clipped <- pmin(plot_data$neg_log10_fdr, fdr_cap)

    plot_data <- plot_data %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(
        Description = str_wrap(Description, width = 45),
        Description = factor(Description, levels = rev(Description))
      )

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
      geom_point(aes(size = Count, color = neg_log10_fdr_clipped)) +
      scale_color_gradientn(
        colors = enrichment_colors$gradient,
        name = expression(-log[10]*"(FDR)")
      ) +
      scale_size_continuous(name = "Count", range = c(2, 6)) +
      scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      labs(
        title = paste0(database, " - ", direction_label, "\n", contrast_name),
        x = "Gene Ratio",
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      ) +
      guides(
        color = guide_colorbar(order = 1, barwidth = 0.8, barheight = 3),
        size = guide_legend(order = 2)
      )

    db_short <- gsub(":", "", database)
    plot_file <- paste0("Dotplot_", db_short, "_", contrast_name, "_", direction, "_", run_tag, ".png")
    ggsave(file.path(output_dir, plot_file), plot = p,
           width = 6.5, height = max(4, nrow(plot_data) * 0.25 + 1.5), dpi = 300, bg = "white")
    cat("[OK] Saved: ", plot_file, "\n", sep = "")
    return(p)
  }

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
## SECTION C: EXPRESSION HEATMAPS (Z-scored, individual samples)
## ============================================================

if (generate_individual_heatmaps) {
  cat("\n=== SECTION C: HEATMAPS ===\n")

  ## Load VST data
  vst_files <- list.files(tables_dir, pattern = "^VST_matrix_heatmap_.*\\.rds$", full.names = TRUE)
  vst_data <- NULL
  if (length(vst_files) > 0) {
    vst_file <- vst_files[order(file.info(vst_files)$mtime, decreasing = TRUE)][1]
    vst_data <- readRDS(vst_file)
    cat("[OK] Loaded VST matrix from: ", basename(vst_file), "\n", sep = "")
  } else {
    cat("[WARN] VST matrix not found - skipping heatmaps\n")
  }

  if (!is.null(vst_data)) {
    ## Extract gene categories from fGSEA results
    fgsea_gobp_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)

    extract_leading_edge <- function(fgsea_file, pathway_patterns, max_genes = 15) {
      if (!file.exists(fgsea_file)) return(character(0))
      fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
      genes <- character(0)
      for (pattern in pathway_patterns) {
        matched <- grep(pattern, fgsea_data$pathway_name, ignore.case = TRUE)
        if (length(matched) > 0) {
          le <- fgsea_data$leadingEdge_symbols[matched[1]]
          if (!is.na(le) && le != "") {
            genes <- c(genes, head(unlist(strsplit(le, ",")), max_genes))
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
        "respiratory chain", "fatty acid beta-oxidation"
      ), 10)
      if (length(mito_genes) > 0) gene_categories[["Mitochondrial"]] <- head(mito_genes, 20)
    }

    gwat_fgsea <- fgsea_gobp_files[grep("gWAT", fgsea_gobp_files, ignore.case = TRUE)]
    if (length(gwat_fgsea) > 0) {
      immune_genes <- extract_leading_edge(gwat_fgsea[1], c(
        "leukocyte degranulation", "myeloid cell activation", "cytokine"
      ), 10)
      if (length(immune_genes) > 0) gene_categories[["Immune"]] <- head(immune_genes, 20)

      lipid_genes <- extract_leading_edge(gwat_fgsea[1], c(
        "lipid localization", "lipid catabolic", "fatty acid"
      ), 10)
      if (length(lipid_genes) > 0) gene_categories[["Lipid"]] <- head(lipid_genes, 15)
    }

    if (length(gene_categories) == 0) {
      gene_categories <- list(
        "Inflammatory" = c("Il1b", "Il6", "Tnf", "Ccl2", "Cxcl12"),
        "Collagen" = c("Col1a1", "Col1a2", "Col3a1", "Col4a1"),
        "Adipocyte" = c("Lep", "Adipoq", "Pparg", "Fabp4", "Plin1")
      )
    }

    ## Create heatmap function
    create_individual_heatmap <- function(vst_data, contrast_name, gene_categories,
                                           output_dir = plots_dir) {

      cat("\n  Processing: ", contrast_name, "\n", sep = "")

      vst_mat <- vst_data$vst_mat
      sample_info <- vst_data$sample_info

      tissue <- sub("_.*", "", contrast_name)

      if ("Depot" %in% colnames(sample_info)) {
        tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
      } else {
        tissue_samples <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
      }

      if (length(tissue_samples) == 0) return(NULL)

      vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples, drop = FALSE]
      sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

      all_genes <- unique(unlist(gene_categories))
      genes_in_data <- intersect(all_genes, rownames(vst_tissue))
      if (length(genes_in_data) < 3) return(NULL)

      gene_to_cat <- character(length(genes_in_data))
      names(gene_to_cat) <- genes_in_data
      for (cat_name in names(gene_categories)) {
        for (g in intersect(gene_categories[[cat_name]], genes_in_data)) {
          if (gene_to_cat[g] == "") gene_to_cat[g] <- cat_name
        }
      }
      gene_to_cat <- gene_to_cat[gene_to_cat != ""]
      genes_to_plot <- names(gene_to_cat)
      if (length(genes_to_plot) < 3) return(NULL)

      mat <- vst_tissue[genes_to_plot, , drop = FALSE]

      gene_order <- c()
      for (cat_name in names(gene_categories)) {
        gene_order <- c(gene_order, names(gene_to_cat[gene_to_cat == cat_name]))
      }
      mat <- mat[gene_order, , drop = FALSE]
      gene_to_cat <- gene_to_cat[gene_order]

      if ("Genotype" %in% colnames(sample_info_tissue)) {
        sample_order <- order(sample_info_tissue$Genotype)
        mat <- mat[, sample_order, drop = FALSE]
        sample_info_tissue <- sample_info_tissue[sample_order, , drop = FALSE]
      }

      ## Z-score normalization (row-wise)
      mat <- t(scale(t(mat)))
      z_clip <- heatmap_params$lfc_clip
      mat[mat > z_clip] <- z_clip
      mat[mat < -z_clip] <- -z_clip

      ## Color scale using custom gradient
      col_fun <- colorRamp2(
        seq(-z_clip, z_clip, length.out = 5),
        heatmap_params$custom_gradient
      )

      col_anno <- NULL
      if ("Genotype" %in% colnames(sample_info_tissue)) {
        col_anno <- HeatmapAnnotation(
          Genotype = sample_info_tissue$Genotype,
          col = list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02")),
          annotation_name_side = "left",
          simple_anno_size = unit(4, "mm")
        )
      }

      ht <- Heatmap(
        mat,
        col = col_fun,
        name = "Z-score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_rot = 45,
        top_annotation = col_anno,
        row_split = factor(gene_to_cat, levels = names(gene_categories)),
        row_gap = unit(1.5, "mm"),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 9, fontface = "bold"),
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        column_names_gp = gpar(fontsize = 8),
        column_title = contrast_name,
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        border = TRUE,
        heatmap_legend_param = list(
          title = "Z-score",
          title_gp = gpar(fontsize = 9),
          labels_gp = gpar(fontsize = 8),
          legend_height = unit(2.5, "cm")
        )
      )

      plot_file <- paste0("Heatmap_", contrast_name, "_", run_tag, ".png")
      n_genes <- length(genes_to_plot)
      n_samples <- ncol(mat)

      png(file.path(output_dir, plot_file),
          width = max(500, n_samples * 40 + 200),
          height = max(400, n_genes * 16 + 100),
          res = 150)
      draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
      dev.off()

      cat("[OK] Saved: ", plot_file, "\n", sep = "")
      return(ht)
    }

    ## Process each contrast
    de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

    for (de_file in de_files) {
      filename <- basename(de_file)
      contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (!is.na(contrast)) {
        create_individual_heatmap(vst_data, contrast, gene_categories)
      }
    }
  }
}

## ============================================================
## COMPLETE
## ============================================================

cat("\n======================================\n")
cat("=== PART 4 (VISUALIZATIONS) COMPLETE ===\n")
cat("======================================\n")
cat("Generated:\n")
if (generate_pca_plot) cat("  - PCA plot\n")
if (generate_mds_plot) cat("  - MDS plot\n")
if (generate_density_plot) cat("  - Density plot\n")
if (generate_ora_barplots) cat("  - ORA bar plots\n")
if (generate_fgsea_barplots) cat("  - fGSEA bar plots\n")
if (generate_ora_dotplots) cat("  - ORA dot plots\n")
if (generate_individual_heatmaps) cat("  - Individual sample heatmaps (z-scored)\n")
cat("\nAll plots saved to: ", plots_dir, "\n\n", sep = "")
