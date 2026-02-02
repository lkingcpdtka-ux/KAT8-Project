#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 4: ALL VISUALIZATIONS
## =========================================================
## This script generates ALL plots WITHOUT re-running analysis.
## NO PLOTS should be made in parts 1-3!
##
##   SECTION A: QC PLOTS
##     - PCA, MDS, Density
##
##   SECTION B: VOLCANO PLOTS
##
##   SECTION C: ENRICHMENT PLOTS
##     - ORA bar/dot plots, fGSEA bar plots
##
##   SECTION D: HEATMAPS
##     - Top DEG heatmaps (simple, top genes by significance)
##     - Individual sample heatmaps (pathway-grouped)
##
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
  stop("parameters.R not found.")
}

## ============================================================
## VISUALIZATION TOGGLES
## ============================================================

## Section A: QC Plots
generate_pca_plot       <- TRUE
generate_mds_plot       <- TRUE
generate_density_plot   <- TRUE

## Section B: Volcano Plots
generate_volcano_plots  <- TRUE

## Section C: Enrichment Plots
generate_ora_barplots   <- TRUE
generate_fgsea_barplots <- TRUE
generate_ora_dotplots   <- TRUE

## Section D: Heatmaps
generate_top_deg_heatmaps     <- TRUE   ## Simple top DEGs (no pathway sections)
generate_individual_heatmaps  <- TRUE   ## Individual samples with pathway groups

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

  qc_files <- list.files(tables_dir, pattern = "^QC_data_.*\\.rds$", full.names = TRUE)
  if (length(qc_files) == 0) {
    cat("[WARN] QC data not found. Skipping QC plots.\n")
  } else {
    qc_file <- qc_files[order(file.info(qc_files)$mtime, decreasing = TRUE)][1]
    qc_data <- readRDS(qc_file)
    cat("[OK] Loaded QC data\n")

    vst_mat <- qc_data$vst_mat_qc
    sample_info <- qc_data$sample_info
    depot_sex_levels <- unique(sample_info$DepotSex)
    if (!is.null(qc_plot_params$depot_sex_manual) &&
        qc_plot_params$depot_sex_palette == "manual") {
      depot_sex_fill <- qc_plot_params$depot_sex_manual[depot_sex_levels]
    } else if (!is.null(qc_plot_params$depot_sex_palette)) {
      depot_sex_fill <- setNames(
        viridis(length(depot_sex_levels), option = qc_plot_params$depot_sex_palette),
        depot_sex_levels
      )
    } else {
      depot_sex_fill <- qc_data$depot_sex_fill
    }

    ## A.1) PCA Plot
    if (generate_pca_plot && !is.null(vst_mat)) {
      cat("--- Creating PCA plot ---\n")
      pca_result <- prcomp(t(vst_mat), scale. = TRUE)
      pca_var <- summary(pca_result)$importance[2, 1:2] * 100

      pca_df <- data.frame(
        Sample = colnames(vst_mat), Genotype = sample_info$Genotype,
        DepotSex = sample_info$DepotSex,
        PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2]
      )

      p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = DepotSex, fill = DepotSex,
                                   shape = Genotype, label = Sample)) +
        geom_point(size = 3, stroke = 1) +
        geom_text_repel(size = 3, max.overlaps = 60, color = "black") +
        scale_color_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_fill_manual(values = depot_sex_fill, name = "Depot/Sex") +
        scale_shape_manual(values = c("CTL" = 21, "KAT8KD" = 24), name = "Genotype") +
        labs(title = "PCA", x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
             y = paste0("PC2 (", round(pca_var[2], 1), "%)")) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        guides(color = "none", fill = guide_legend(override.aes = list(shape = 21)),
               shape = guide_legend(override.aes = list(fill = "white", color = "black")))

      ggsave(file.path(plots_dir, paste0("PCA_", run_tag, ".png")), p_pca, width = 8, height = 6, dpi = 300)
      cat("[OK] Saved PCA plot\n")
    }

    ## A.2) MDS Plot
    if (generate_mds_plot && !is.null(vst_mat)) {
      cat("--- Creating MDS plot ---\n")
      dist_mat <- dist(t(vst_mat))
      mds_coords <- cmdscale(dist_mat, k = 2)

      mds_df <- data.frame(
        Sample = colnames(vst_mat), Genotype = sample_info$Genotype,
        DepotSex = sample_info$DepotSex,
        MDS1 = mds_coords[, 1], MDS2 = mds_coords[, 2]
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
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        guides(color = "none", fill = guide_legend(override.aes = list(shape = 21)),
               shape = guide_legend(override.aes = list(fill = "white", color = "black")))

      ggsave(file.path(plots_dir, paste0("MDS_", run_tag, ".png")), p_mds, width = 8, height = 6, dpi = 300)
      cat("[OK] Saved MDS plot\n")
    }

    ## A.3) Density Plot
    if (generate_density_plot && !is.null(qc_data$vst_mat_before)) {
      cat("--- Creating density plot ---\n")
      vst_before <- qc_data$vst_mat_before
      vst_after <- vst_mat

      png(file.path(plots_dir, paste0("Density_", run_tag, ".png")), width = 1600, height = 900, res = 150)
      par(mfrow = c(1, 2), mar = c(5, 5, 5, 2))

      dens_before <- apply(vst_before, 2, density)
      dens_after  <- apply(vst_after, 2, density)
      xlim_global <- range(c(range(sapply(dens_before, function(d) range(d$x))),
                              range(sapply(dens_after, function(d) range(d$x)))))
      sample_colors <- depot_sex_fill[sample_info$DepotSex]

      plot(dens_before[[1]], main = "Before Filtering", xlab = "VST", lwd = 2,
           col = sample_colors[1], xlim = xlim_global)
      for (i in 2:ncol(vst_before)) lines(dens_before[[i]], lwd = 2, col = sample_colors[i])

      plot(dens_after[[1]], main = "After Filtering", xlab = "VST", lwd = 2,
           col = sample_colors[1], xlim = xlim_global)
      for (i in 2:ncol(vst_after)) lines(dens_after[[i]], lwd = 2, col = sample_colors[i])

      dev.off()
      cat("[OK] Saved density plot\n")
    }
  }
}

## ============================================================
## SECTION B: VOLCANO PLOTS
## ============================================================

if (generate_volcano_plots) {
  cat("\n=== SECTION B: VOLCANO PLOTS ===\n")
  de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

  for (de_file in de_files) {
    filename <- basename(de_file)
    contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
    if (is.na(contrast)) next

    cat("--- Volcano: ", contrast, " ---\n", sep = "")

    tt <- read.csv(de_file, stringsAsFactors = FALSE)

    ## Get thresholds
    thresholds <- get_tissue_thresholds(contrast)
    cn_logFC_cut <- thresholds$logFC_cut
    cn_fdr_cut <- thresholds$fdr_cut

    tt_plot <- tt %>%
      dplyr::filter(is.finite(logFC), is.finite(padj), padj > 0) %>%
      dplyr::mutate(
        negLogFDR = -log10(padj),
        sig_cat = dplyr::case_when(
          padj < cn_fdr_cut & logFC >=  cn_logFC_cut ~ "Up",
          padj < cn_fdr_cut & logFC <= -cn_logFC_cut ~ "Down",
          TRUE ~ "NS"
        ),
        sig_cat = factor(sig_cat, levels = c("Up", "Down", "NS"))
      )

    ## Get top genes to label
    tt_sig <- tt_plot %>% dplyr::filter(sig_cat != "NS")

    top_up <- tt_sig %>% dplyr::filter(sig_cat == "Up") %>%
      dplyr::arrange(padj) %>% dplyr::slice_head(n = 10)
    top_down <- tt_sig %>% dplyr::filter(sig_cat == "Down") %>%
      dplyr::arrange(padj) %>% dplyr::slice_head(n = 10)
    label_genes <- dplyr::bind_rows(top_up, top_down)

    p <- ggplot(tt_plot, aes(x = logFC, y = negLogFDR, color = sig_cat)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(
        values = c("Up" = volcano_colors$up, "Down" = volcano_colors$down, "NS" = volcano_colors$ns),
        name = "Significance"
      ) +
      geom_point(data = label_genes, size = 3) +
      geom_text_repel(data = label_genes, aes(label = gene_name), color = "black",
                      size = 3.5, box.padding = 0.25, max.overlaps = Inf) +
      geom_hline(yintercept = -log10(cn_fdr_cut), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(cn_logFC_cut, -cn_logFC_cut), linetype = "dashed", color = "grey40") +
      theme_bw(base_size = 14) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5)) +
      labs(title = paste0("Volcano: ", contrast),
           x = expression(log[2]~Fold~Change), y = expression(-log[10]~FDR))

    ggsave(file.path(plots_dir, paste0("Volcano_", contrast, "_", run_tag, ".png")),
           p, width = 8, height = 6, dpi = 300)
    cat("[OK] Saved volcano plot\n")
  }
}

## ============================================================
## SECTION C: ENRICHMENT PLOTS
## ============================================================

## C.1) ORA BAR PLOTS
if (generate_ora_barplots) {
  cat("\n=== GENERATING ORA BAR PLOTS ===\n")

  ora_files <- c(
    list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE),
    list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE)
  )

  for (ora_file in ora_files) {
    filename <- basename(ora_file)
    parts <- str_match(filename, "ORA_(gobp|kegg)_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
    if (is.na(parts[1])) next

    db <- ifelse(parts[2] == "gobp", "GO:BP", "KEGG")
    contrast <- parts[3]
    direction <- parts[4]

    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) next

    if (!"pathway_name" %in% colnames(ora_data)) ora_data$pathway_name <- ora_data$Description

    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < ora_params$pvalue_cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = ora_params$top_n)

    if (nrow(plot_data) == 0) next

    plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    plot_data$neg_log10_padj <- -log10(plot_data$p.adjust)
    plot_data <- plot_data %>%
      dplyr::mutate(pathway_name = factor(pathway_name, levels = rev(pathway_name)))

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = pathway_name, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_gradientn(colors = enrichment_colors$gradient, name = expression(-log[10]*"(FDR)")) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
      labs(title = paste0(db, " - ", direction_label, "\n", contrast), x = "Gene Ratio", y = NULL) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.margin = margin(6, 10, 6, 10))

    db_short <- gsub(":", "", db)
    ggsave(file.path(plots_dir, paste0("ORA_", db_short, "_", contrast, "_", direction, "_", run_tag, ".png")),
           p, width = barplot_params$plot_width, height = barplot_params$plot_height,
           dpi = 300, bg = "white")
    cat("[OK] Saved ORA bar plot: ", db, " ", contrast, " ", direction, "\n", sep = "")
  }
}

## C.2) fGSEA BAR PLOTS
if (generate_fgsea_barplots) {
  cat("\n=== GENERATING fGSEA BAR PLOTS ===\n")

  fgsea_files <- c(
    list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE),
    list.files(tables_dir, pattern = "^fgsea_kegg_.*\\.csv$", full.names = TRUE)
  )

  for (fgsea_file in fgsea_files) {
    filename <- basename(fgsea_file)
    parts <- str_match(filename, "fgsea_(gobp|kegg)_(.+)_[0-9]+_[0-9]+\\.csv")
    if (is.na(parts[1])) next

    db <- ifelse(parts[2] == "gobp", "GO:BP", "KEGG")
    contrast <- parts[3]

    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
    sig_data <- fgsea_data %>% dplyr::filter(padj < gsea_params$fdr_cutoff)
    if (nrow(sig_data) == 0) next

    top_n <- gsea_params$top_n_per_direction
    top_up <- sig_data %>% dplyr::filter(NES > 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n)
    top_down <- sig_data %>% dplyr::filter(NES < 0) %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n)

    plot_data <- dplyr::bind_rows(top_down, top_up) %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(
        pathway_label = ifelse(!is.na(pathway_name) & pathway_name != "", pathway_name, pathway_id),
        pathway_label = factor(pathway_label, levels = pathway_label),
        neg_log10_padj = -log10(pmax(padj, .Machine$double.eps))
      )

    if (nrow(plot_data) == 0) next

    p <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.75) +
      scale_fill_gradientn(colors = enrichment_colors$gradient, name = expression(-log[10]*"(FDR)")) +
      scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
      labs(title = paste0("GSEA ", db, "\n", contrast), x = "NES", y = NULL) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.margin = margin(6, 10, 6, 10)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

    db_short <- tolower(gsub(":", "", db))
    ggsave(file.path(plots_dir, paste0("fGSEA_", db_short, "_", contrast, "_", run_tag, ".png")),
           p,
           width = gsea_plot_params$plot_width,
           height = gsea_plot_params$plot_height,
           dpi = 300, bg = "white")
    cat("[OK] Saved fGSEA plot: ", db, " ", contrast, "\n", sep = "")
  }
}

## C.3) ORA DOT PLOTS
if (generate_ora_dotplots) {
  cat("\n=== GENERATING ORA DOT PLOTS ===\n")

  ora_files <- c(
    list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE),
    list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE)
  )

  for (ora_file in ora_files) {
    filename <- basename(ora_file)
    parts <- str_match(filename, "ORA_(gobp|kegg)_(.+)_(Up|Down)_[0-9]+_[0-9]+\\.csv")
    if (is.na(parts[1])) next

    db <- ifelse(parts[2] == "gobp", "GO:BP", "KEGG")
    contrast <- parts[3]
    direction <- parts[4]

    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) next

    if (!"Description" %in% colnames(ora_data)) ora_data$Description <- ora_data$pathway_name

    ora_data$GeneRatio_numeric <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    ora_data$Count <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) as.numeric(x[1]))
    ora_data$neg_log10_fdr <- -log10(ora_data$p.adjust)

    plot_data <- ora_data %>%
      dplyr::filter(p.adjust < dotplot_params$fdr_cutoff, Count >= dotplot_params$min_count) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = dotplot_params$top_n_gobp)

    if (nrow(plot_data) == 0) next

    plot_data <- plot_data %>%
      dplyr::mutate(Description = str_wrap(Description, width = 32),
                    Description = factor(Description, levels = rev(Description)))

    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"

    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
      geom_point(aes(size = Count, color = neg_log10_fdr)) +
      scale_color_gradientn(colors = enrichment_colors$gradient, name = expression(-log[10]*"(FDR)")) +
      scale_size_continuous(name = "Count", range = c(2, 6)) +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
      labs(title = paste0(db, " - ", direction_label, "\n", contrast), x = "Gene Ratio", y = NULL) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
            panel.grid.minor = element_blank(),
            plot.margin = margin(6, 10, 6, 10))

    db_short <- gsub(":", "", db)
    ggsave(file.path(plots_dir, paste0("Dotplot_", db_short, "_", contrast, "_", direction, "_", run_tag, ".png")),
           p, width = dotplot_params$plot_width, height = dotplot_params$plot_height,
           dpi = 300, bg = "white")
    cat("[OK] Saved dot plot: ", db, " ", contrast, " ", direction, "\n", sep = "")
  }
}

## ============================================================
## SECTION D: HEATMAPS
## ============================================================

cat("\n=== SECTION D: HEATMAPS ===\n")

## Load VST data
vst_files <- list.files(tables_dir, pattern = "^VST_matrix_heatmap_.*\\.rds$", full.names = TRUE)
vst_data <- NULL
if (length(vst_files) > 0) {
  vst_file <- vst_files[order(file.info(vst_files)$mtime, decreasing = TRUE)][1]
  vst_data <- readRDS(vst_file)
  cat("[OK] Loaded VST matrix\n")
}

if (!is.null(vst_data)) {

  ## D.1) TOP DEG HEATMAPS (Simple - no pathway sections)
  if (generate_top_deg_heatmaps) {
    cat("\n--- Creating Top DEG Heatmaps ---\n")

    de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

    for (de_file in de_files) {
      filename <- basename(de_file)
      contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (is.na(contrast)) next

      cat("  Top DEG heatmap: ", contrast, "\n", sep = "")

      tt <- read.csv(de_file, stringsAsFactors = FALSE)
      tissue <- sub("_.*", "", contrast)

      ## Get samples for this tissue
      sample_info <- vst_data$sample_info
      vst_mat <- vst_data$vst_mat

      if ("Depot" %in% colnames(sample_info)) {
        tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
      } else {
        tissue_samples <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
      }

      if (length(tissue_samples) == 0) next

      vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples, drop = FALSE]
      sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

      ## Get top DEGs by significance
      thresholds <- get_tissue_thresholds(contrast)
      top_degs <- tt %>%
        dplyr::filter(padj < thresholds$fdr_cut, abs(logFC) >= thresholds$logFC_cut) %>%
        dplyr::arrange(padj) %>%
        dplyr::slice_head(n = 100)

      if (nrow(top_degs) < 5) next

      genes_to_plot <- intersect(top_degs$gene_name, rownames(vst_tissue))
      if (length(genes_to_plot) < 5) next

      mat <- vst_tissue[genes_to_plot, , drop = FALSE]

      ## Order samples by genotype
      if ("Genotype" %in% colnames(sample_info_tissue)) {
        sample_order <- order(sample_info_tissue$Genotype)
        mat <- mat[, sample_order, drop = FALSE]
        sample_info_tissue <- sample_info_tissue[sample_order, , drop = FALSE]
      }

      ## Z-score normalize
      mat <- t(scale(t(mat)))
      z_clip <- heatmap_params$lfc_clip
      mat[mat > z_clip] <- z_clip
      mat[mat < -z_clip] <- -z_clip

      col_fun <- colorRamp2(seq(-z_clip, z_clip, length.out = 5), heatmap_params$custom_gradient)

      col_anno <- NULL
      if ("Genotype" %in% colnames(sample_info_tissue)) {
        col_anno <- HeatmapAnnotation(
          Genotype = sample_info_tissue$Genotype,
          col = list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02")),
          simple_anno_size = unit(4, "mm")
        )
      }

      ht <- Heatmap(
        mat, col = col_fun, name = "Z-score",
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE,
        column_names_rot = 45,
        top_annotation = col_anno,
        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
        column_names_gp = gpar(fontsize = 8),
        column_title = paste0("Top DEGs: ", contrast),
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        border = TRUE,
        heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 9))
      )

      plot_file <- paste0("Heatmap_TopDEG_", contrast, "_", run_tag, ".png")
      png(file.path(plots_dir, plot_file),
          width = max(500, ncol(mat) * 40 + 200),
          height = max(400, nrow(mat) * 12 + 100), res = 150)
      draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
      dev.off()
      cat("[OK] Saved: ", plot_file, "\n", sep = "")
    }
  }

  ## D.2) INDIVIDUAL SAMPLE HEATMAPS (with pathway groups)
  if (generate_individual_heatmaps) {
    cat("\n--- Creating Individual Sample Heatmaps (pathway-grouped) ---\n")

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
        "respiratory chain complex", "fatty acid beta-oxidation"
      ), 10)
      if (length(mito_genes) > 0) {
        gene_categories[["Mitochondrial /\nEnergy"]] <- head(mito_genes, 20)
      }
    }

    gwat_fgsea <- fgsea_gobp_files[grep("gWAT", fgsea_gobp_files, ignore.case = TRUE)]
    if (length(gwat_fgsea) > 0) {
      immune_genes <- extract_leading_edge(gwat_fgsea[1], c(
        "leukocyte degranulation", "myeloid cell activation", "positive regulation of cytokine"
      ), 10)
      if (length(immune_genes) > 0) {
        gene_categories[["Immune /\nInflammatory"]] <- head(immune_genes, 20)
      }

      lipid_genes <- extract_leading_edge(gwat_fgsea[1], c(
        "lipid localization", "lipid catabolic", "fatty acid"
      ), 10)
      if (length(lipid_genes) > 0) {
        gene_categories[["Lipid\nMetabolism"]] <- head(lipid_genes, 15)
      }

      ecm_genes <- extract_leading_edge(gwat_fgsea[1], c(
        "tissue remodeling", "extracellular matrix", "cell adhesion"
      ), 10)
      if (length(ecm_genes) > 0) {
        gene_categories[["ECM /\nRemodeling"]] <- head(ecm_genes, 15)
      }
    }

    if (length(gene_categories) == 0) {
      gene_categories <- list(
        "Inflammatory\nResponse" = c("Il1b", "Il6", "Tnf", "Il1a", "Ccl2", "Ccl7", "Cxcl12"),
        "ECM /\nCollagens" = c("Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col5a1", "Col6a1"),
        "Adipocyte\nMarkers" = c("Lep", "Adipoq", "Pparg", "Ppargc1a", "Fabp4", "Plin1")
      )
    }

    de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)

    for (de_file in de_files) {
      filename <- basename(de_file)
      contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (is.na(contrast)) next

      cat("  Individual heatmap: ", contrast, "\n", sep = "")

      tissue <- sub("_.*", "", contrast)
      sample_info <- vst_data$sample_info
      vst_mat <- vst_data$vst_mat

      if ("Depot" %in% colnames(sample_info)) {
        tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
      } else {
        tissue_samples <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
      }

      if (length(tissue_samples) == 0) next

      vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples, drop = FALSE]
      sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]

      all_genes <- unique(unlist(gene_categories))
      genes_in_data <- intersect(all_genes, rownames(vst_tissue))
      if (length(genes_in_data) < 3) next

      gene_to_cat <- character(length(genes_in_data))
      names(gene_to_cat) <- genes_in_data
      for (cat_name in names(gene_categories)) {
        for (g in intersect(gene_categories[[cat_name]], genes_in_data)) {
          if (gene_to_cat[g] == "") gene_to_cat[g] <- cat_name
        }
      }
      gene_to_cat <- gene_to_cat[gene_to_cat != ""]
      genes_to_plot <- names(gene_to_cat)
      if (length(genes_to_plot) < 3) next

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

      mat <- t(scale(t(mat)))
      z_clip <- heatmap_params$lfc_clip
      mat[mat > z_clip] <- z_clip
      mat[mat < -z_clip] <- -z_clip

      col_fun <- colorRamp2(seq(-z_clip, z_clip, length.out = 5), heatmap_params$custom_gradient)

      col_anno <- NULL
      if ("Genotype" %in% colnames(sample_info_tissue)) {
        col_anno <- HeatmapAnnotation(
          Genotype = sample_info_tissue$Genotype,
          col = list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02")),
          simple_anno_size = unit(4, "mm")
        )
      }

      ht <- Heatmap(
        mat, col = col_fun, name = "Z-score",
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE,
        column_names_rot = 45,
        top_annotation = col_anno,
        row_split = factor(gene_to_cat, levels = names(gene_categories)),
        row_gap = unit(1.5, "mm"),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 9, fontface = "bold"),
        row_names_gp = gpar(fontsize = 8, fontface = "italic"),
        column_names_gp = gpar(fontsize = 8),
        column_title = contrast,
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        border = TRUE,
        heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 9))
      )

      plot_file <- paste0("Heatmap_individual_", contrast, "_", run_tag, ".png")
      png(file.path(plots_dir, plot_file),
          width = max(500, ncol(mat) * 40 + 200),
          height = max(400, nrow(mat) * 16 + 100), res = 150)
      draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
      dev.off()
      cat("[OK] Saved: ", plot_file, "\n", sep = "")
    }
  }
}

## ============================================================
## COMPLETE
## ============================================================

cat("\n======================================\n")
cat("=== PART 4 (ALL VISUALIZATIONS) COMPLETE ===\n")
cat("======================================\n")
cat("All plots saved to: ", plots_dir, "\n\n", sep = "")
