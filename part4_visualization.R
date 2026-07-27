#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   Sections E/F: box plots, forest panels, combined + cross-depot
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## KAT8 bulk RNA-seq - PART 4: ALL VISUALIZATIONS
## =========================================================
## This script generates ALL plots WITHOUT re-running analysis.
## NO PLOTS should be made in parts 1-3!
##
## Part 1 uses a unified DESeq2 model: ~ 0 + GroupDepot
## (all 38 tissue samples, per-depot KAT8 effects via contrasts)
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
generate_avg4col_heatmaps     <- TRUE   ## Averaged per genotype x sex (4-column summary)

## Per-gene expression box plots for GENES_OF_INTEREST (Section E)
generate_gene_boxplots        <- TRUE

## Marker-panel forest plots: log2FC + CI per gene (Section F)
generate_marker_forest_plots  <- TRUE

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

## Resolve the run folder Part 1 created.
## BUGFIX: previously selected by file mtime. This project lives in OneDrive and
## file sync updates mtimes, which can silently promote a STALE run folder --
## causing this script to read the wrong DE tables and write results into the
## wrong run. Prefer the explicit pointer written by Part 1; fall back to a
## FOLDER-NAME sort (RUN_YYYYMMDD_HHMMSS sorts chronologically and is immune
## to mtime churn). Never sort by mtime.
outdir <- NULL
latest_ptr <- file.path(savepoint_dir, "LATEST_RUN.txt")
if (file.exists(latest_ptr)) {
  cand <- trimws(readLines(latest_ptr, warn = FALSE))
  cand <- cand[nzchar(cand)]
  if (length(cand) >= 1 && dir.exists(cand[1])) {
    outdir <- cand[1]
    cat("[OK] Run folder from LATEST_RUN.txt: ", outdir, "\n", sep = "")
  } else {
    cat("[WARN] LATEST_RUN.txt present but path missing; falling back to name sort\n")
  }
}
if (is.null(outdir)) {
  run_dirs_sorted <- run_dirs[order(basename(run_dirs), decreasing = TRUE)]
  outdir <- run_dirs_sorted[1]
  cat("[INFO] Using newest RUN_ folder by name: ", outdir, "\n", sep = "")
}
## Derive a clean run tag: the run folder may be nested (RUN_<tag>/tissue) and
## its name may carry a label (RUN_TISSUE_PANELS_20260515_163857). Both broke
## downstream filename parsing -- run_tag became "tissue", and contrast names
## picked up "_TISSUE_PANELS", which stopped the authoritative DEG lists from
## matching. Take the trailing YYYYMMDD_HHMMSS when present.
.clean_run_tag <- function(dir) {
  b <- basename(dir)
  if (b %in% c("tissue", "cells")) b <- basename(dirname(dir))
  b <- sub("^RUN_", "", b)
  m <- regmatches(b, regexpr("[0-9]{8}_[0-9]{6}$", b))
  if (length(m) == 1 && nzchar(m)) m else b
}
run_tag <- .clean_run_tag(outdir)

cat("[INFO] Using results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")

## ---- ONE CONTRAST PER DEPOT ---------------------------------------------
## Everything here globs "DE_tissue_*.csv". A run folder that also holds DE
## tables written by an older version of Part 1 under a different naming
## scheme therefore yielded TWO contrasts per depot -- iWAT_KD_vs_CTL and
## iWAT_KD_vs_CTL_TISSUE_PANELS -- so every volcano, heatmap, box plot and
## forest plot was produced TWICE, half of them from stale input.
##
## Part 1 now moves superseded tables aside, but that only helps on a run
## where Part 1 actually executes. This is the backstop: keep one file per
## depot, preferring the plainest contrast name (a decorated name is the
## leftover), and say out loud what was ignored.
.unique_de_files <- function(files) {
  if (length(files) < 2) return(files)
  cn  <- sub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(files))
  dep <- sub("_.*$", "", cn)
  keep <- logical(length(files))
  for (d in unique(dep)) {
    idx <- which(dep == d)
    best <- idx[order(nchar(cn[idx]),
                      -as.numeric(file.info(files[idx])$mtime))][1]
    keep[best] <- TRUE
    if (length(idx) > 1)
      cat("[WARN] ", length(idx), " DE tables for ", d, "; using '", cn[best],
          "', ignoring: ", paste(cn[setdiff(idx, best)], collapse = ", "), "\n", sep = "")
  }
  files[keep]
}

## ---- SIDE-BY-SIDE FIGURE LAYOUT -----------------------------------------
## Two finished ggplots as one image, each KEEPING ITS OWN axis and colour
## scale -- iWAT effects reach |log2FC| ~ 2 while gWAT sits under 1, and a
## shared axis would squash the smaller panel into a stripe.
## patchwork is preferred; gridExtra is the fallback; if neither is installed
## the single-panel figures are still written and only this extra is skipped.
.save_side_by_side <- function(plots, path, title = NULL, subtitle = NULL,
                               width = 15, height = 8) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- Reduce(`|`, plots)
    if (!is.null(title))
      comb <- comb + patchwork::plot_annotation(
        title = title, subtitle = subtitle,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 15),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35")))
    ggsave(path, comb, width = width, height = height, dpi = 300,
           bg = "white", limitsize = FALSE)
    return(TRUE)
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    png(path, width = width, height = height, units = "in", res = 300, bg = "white")
    do.call(gridExtra::grid.arrange,
            c(plots, list(ncol = length(plots), top = title)))
    dev.off()
    return(TRUE)
  }
  cat("[SKIP] side-by-side figure needs the 'patchwork' package:",
      " install.packages(\"patchwork\")\n", sep = "")
  FALSE
}

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
    ## Align the annotation to the MATRIX COLUMNS explicitly. The QC plots read
    ## Genotype/DepotSex positionally, so if the row order ever differed from
    ## colnames(vst_mat) every point would be mislabelled without any error.
    if (!is.null(vst_mat) && "Sample" %in% colnames(sample_info)) {
      .ord <- match(colnames(vst_mat), sample_info$Sample)
      if (anyNA(.ord)) {
        stop("QC annotation does not cover every sample in the VST matrix: ",
             paste(setdiff(colnames(vst_mat), sample_info$Sample), collapse = ", "),
             call. = FALSE)
      }
      sample_info <- sample_info[.ord, , drop = FALSE]
      cat("[OK] QC annotation aligned to ", ncol(vst_mat), " matrix columns\n", sep = "")
    }
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
      ## PCA gene selection + scaling are set in parameters.R
      ## (qc_plot_params$pca_ntop / $pca_scale). This pipeline uses
      ## ntop = Inf, scale = TRUE -- every filtered gene, unit-scaled. See the
      ## note in parameters.R for why that convention is appropriate here.
      ## NOTE: the SIGN of a principal component is arbitrary, so a plot can
      ## appear mirrored between runs without anything being wrong.
      .ntop  <- if (!is.null(qc_plot_params$pca_ntop))  qc_plot_params$pca_ntop  else Inf
      .scale <- if (!is.null(qc_plot_params$pca_scale)) qc_plot_params$pca_scale else TRUE
      rv <- apply(vst_mat, 1, stats::var)
      ## scale. = TRUE errors on a constant column, so zero-variance genes must
      ## go first. They carry no information under either convention.
      pca_mat <- vst_mat
      if (isTRUE(.scale)) {
        keep_var <- is.finite(rv) & rv > 0
        if (any(!keep_var)) {
          cat("[INFO] dropped ", sum(!keep_var), " zero-variance genes before scaling\n", sep = "")
          pca_mat <- pca_mat[keep_var, , drop = FALSE]
          rv <- rv[keep_var]
        }
      }
      keep_n <- if (is.infinite(.ntop)) length(rv) else min(.ntop, length(rv))
      top_var_idx <- order(rv, decreasing = TRUE)[seq_len(keep_n)]
      cat("[INFO] PCA on ", keep_n, " genes, scale = ", .scale, "\n", sep = "")
      pca_result <- prcomp(t(pca_mat[top_var_idx, , drop = FALSE]), scale. = .scale)
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
## Gene labeling strategy: ORA-informed selection
## Genes are labeled if they are DEGs AND appear in significant ORA
## pathways. Ranked by composite score:
##   score = z(pathway_count) + z(|logFC|) + z(-log10(FDR))
## When no ORA pathways exist for a direction, falls back to DE ranking.
## See parameters.R  volcano_gene_selection for all tunables.
## ============================================================

## --- Helper: build per-gene pathway-count table from ORA KEGG files ----
##
## Returns a data.frame with columns: gene_name, pathway_count, direction
## where pathway_count = number of significant ORA KEGG pathways
## that contain the gene, counted separately per direction (Up / Down).

## Matched on the CONTRAST, not on a tissue prefix.
##
## This previously took a `tissue` string built by sub("_KD_vs_CTL$", "", contrast).
## Contrasts are named "iWAT_KAT8KD_vs_CTL", where the character before "KD" is
## "8" rather than "_", so that sub() never matched, `tissue` came through as the
## full contrast, and the filename pattern (which still demanded a "_<something>_"
## between tissue and direction) matched nothing. Every volcano then silently
## fell back to DE-only labelling and reported "0 ORA-supported".
##
## Part 2 writes `contrast` and `direction` as COLUMNS in each ORA table, so
## those are used directly and no filename has to be parsed. Filename parsing
## is kept only as a fallback for tables written by older versions.
build_ora_gene_counts <- function(tables_dir, contrast, ora_padj_cutoff) {
  ora_gene_list <- list()

  ## Both KEGG and GO:BP are used. KEGG alone yields very few significant
  ## pathways in adipose, so most DEGs had no pathway support even when the
  ## matching worked; GO:BP is where the adipogenic/metabolic terms are.
  ora_files <- c(
    list.files(tables_dir, pattern = "^ORA_kegg_.*\\.csv$", full.names = TRUE),
    list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)
  )
  ora_files <- ora_files[!grepl("ORA_RELAXED", basename(ora_files))]

  for (f in ora_files) {
    ora <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(ora) == 0) next

    ## Keep only rows belonging to THIS contrast.
    if ("contrast" %in% colnames(ora)) {
      ora <- ora[ora$contrast == contrast, , drop = FALSE]
    } else {
      fn_contrast <- str_match(basename(f),
                               "^ORA_(?:kegg|gobp)_(.+)_(?:Up|Down)_[0-9]{8}_[0-9]{6}\\.csv$")[2]
      if (is.na(fn_contrast) || fn_contrast != contrast) next
    }
    if (nrow(ora) == 0) next

    ora <- ora[is.finite(ora$p.adjust) & ora$p.adjust < ora_padj_cutoff, , drop = FALSE]
    if (nrow(ora) == 0) next

    ## Direction: column first, filename second.
    if ("direction" %in% colnames(ora)) {
      dirs <- as.character(ora$direction)
    } else {
      dirs <- rep(str_match(basename(f), "_(Up|Down)_[0-9]{8}_[0-9]{6}\\.csv$")[2], nrow(ora))
    }

    ## KEGG stores Entrez IDs in geneID; symbols are in geneID_symbols.
    ## enrichGO was run with readable = TRUE, so GO:BP geneID is already symbols.
    sym_col <- if ("geneID_symbols" %in% colnames(ora)) "geneID_symbols" else "geneID"
    for (i in seq_len(nrow(ora))) {
      genes <- trimws(unlist(strsplit(ora[[sym_col]][i], "[,/]")))
      genes <- genes[nzchar(genes)]
      if (!length(genes)) next
      ora_gene_list[[length(ora_gene_list) + 1]] <- data.frame(
        gene_name = genes, direction = dirs[i], stringsAsFactors = FALSE
      )
    }
  }

  if (length(ora_gene_list) == 0) {
    return(data.frame(gene_name = character(), pathway_count = integer(),
                      direction = character(), stringsAsFactors = FALSE))
  }
  
  ora_long <- do.call(rbind, ora_gene_list)
  
  ## Count KEGG pathway memberships per gene per direction
  ora_counts <- ora_long %>%
    dplyr::group_by(gene_name, direction) %>%
    dplyr::summarise(pathway_count = dplyr::n(), .groups = "drop")
  
  return(ora_counts)
}

## --- Helper: select ORA-informed genes for one direction ---------------
##
## de_dir: DE data already filtered to one direction (Up or Down)
## ora_counts_dir: ora_counts filtered to that direction
## top_n: how many genes to return
## min_pathway_count: minimum ORA pathway membership
## fallback_to_de: if no ORA genes found, fall back to DE ranking

select_genes_one_direction <- function(de_dir, ora_counts_dir, top_n,
                                       min_pathway_count, fallback_to_de) {
  if (nrow(de_dir) == 0) return(de_dir[0, ])
  
  ## Drop predicted/unannotated genes (Gm*, Rik, etc.) - not informative on plots
  de_dir <- de_dir %>%
    dplyr::filter(!grepl("^Gm\\d+$", gene_name),
                  !grepl("Rik$", gene_name))
  
  ## Merge DE with ORA counts
  merged <- dplyr::left_join(de_dir, ora_counts_dir, by = "gene_name")
  merged$pathway_count[is.na(merged$pathway_count)] <- 0
  
  ## Filter to ORA-supported genes
  ora_supported <- merged %>%
    dplyr::filter(pathway_count >= min_pathway_count)
  
  if (nrow(ora_supported) < 3 && fallback_to_de) {
    ## Fallback: rank by DE evidence alone
    cat("    [INFO] <3 ORA-supported genes; using DE-only ranking\n")
    ## Balanced fallback: half by significance, half by fold change, so the
    ## labels are not all clustered at the top of the volcano.
    n_half <- max(1, floor(top_n / 2))
    return(dplyr::bind_rows(
        de_dir %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = n_half),
        de_dir %>% dplyr::arrange(dplyr::desc(abs(logFC))) %>% dplyr::slice_head(n = top_n - n_half)
      ) %>% dplyr::distinct(gene_name, .keep_all = TRUE))
  }
  
  ## Compute composite score.
  ## PREVIOUSLY z-scored: scale(-log10(padj)) has an enormous right tail here
  ## (padj reaches 1e-27), so z_sig dwarfed z_fc and z_pw and the labels
  ## collapsed onto the most-significant genes at the top of the volcano --
  ## large fold-change genes were never labelled. Percentile RANKS put the
  ## three criteria on the same 0-1 footing so no single one can dominate.
  n <- nrow(ora_supported)
  scored <- ora_supported %>%
    dplyr::mutate(
      r_pw   = rank(pathway_count,  ties.method = "average") / n,
      r_fc   = rank(abs(logFC),     ties.method = "average") / n,
      r_sig  = rank(-log10(padj),   ties.method = "average") / n,
      composite_score = r_pw + r_fc + r_sig
    ) %>%
    dplyr::arrange(dplyr::desc(composite_score))

  ## GUARANTEE BOTH EXTREMES ARE REPRESENTED. Even a balanced score can drift,
  ## so reserve a share of the labels for the largest fold changes outright --
  ## this is the behaviour of the older volcanoes (top-N by FDR + top-N by FC).
  .fcfrac <- if (exists("volcano_fc_label_fraction")) volcano_fc_label_fraction else 0.5
  n_fc   <- max(1, floor(top_n * .fcfrac))
  by_fc  <- ora_supported %>% dplyr::arrange(dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = n_fc)
  by_comp <- scored %>% dplyr::slice_head(n = top_n)
  out <- dplyr::bind_rows(by_comp, by_fc) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::slice_head(n = top_n + n_fc)

  return(out)
}


if (generate_volcano_plots) {
  cat("\n=== SECTION B: VOLCANO PLOTS (ORA-informed labeling) ===\n")
  de_files <- .unique_de_files(list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE))
  
  vs <- volcano_gene_selection  ## shorthand
  
  for (de_file in de_files) {
    filename <- basename(de_file)
    contrast <- str_match(filename, "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
    if (is.na(contrast)) next
    
    cat("--- Volcano: ", contrast, " ---\n", sep = "")
    
    ## Tissue label for plot titles only. Read the depot off the front of the
    ## contrast rather than trying to strip a suffix -- suffixes vary
    ## ("_KAT8KD_vs_CTL", "_KD_vs_CTL") and a failed strip used to leak the whole
    ## contrast into downstream file matching.
    tissue <- sub("_.*$", "", contrast)

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
    
    ## ----- ORA-informed gene selection -----
    ora_counts <- build_ora_gene_counts(tables_dir, contrast, vs$ora_padj_cutoff)
    if (nrow(ora_counts) == 0)
      cat("  [WARN] no ORA tables matched contrast '", contrast,
          "' - labels will use DE evidence only\n", sep = "")
    
    tt_sig <- tt_plot %>% dplyr::filter(sig_cat != "NS")
    
    ## Select top genes per direction
    top_up <- select_genes_one_direction(
      de_dir           = tt_sig %>% dplyr::filter(sig_cat == "Up"),
      ora_counts_dir   = ora_counts %>% dplyr::filter(direction == "Up"),
      top_n            = vs$top_n_per_direction,
      min_pathway_count = vs$min_pathway_count,
      fallback_to_de   = vs$fallback_to_de
    )
    top_down <- select_genes_one_direction(
      de_dir           = tt_sig %>% dplyr::filter(sig_cat == "Down"),
      ora_counts_dir   = ora_counts %>% dplyr::filter(direction == "Down"),
      top_n            = vs$top_n_per_direction,
      min_pathway_count = vs$min_pathway_count,
      fallback_to_de   = vs$fallback_to_de
    )
    
    label_genes <- dplyr::bind_rows(top_up, top_down) %>%
      dplyr::distinct(gene_name, .keep_all = TRUE)
    
    ## Log selection summary
    n_up_ora <- sum(top_up$pathway_count > 0, na.rm = TRUE)
    n_dn_ora <- sum(top_down$pathway_count > 0, na.rm = TRUE)
    cat("  Labels: ", nrow(top_up), " Up (", n_up_ora, " ORA-supported), ",
        nrow(top_down), " Down (", n_dn_ora, " ORA-supported)\n", sep = "")
    
    ## Mandatory genes (biological anchors) - add if they pass DE thresholds
    mandatory <- vs$mandatory_genes
    mandatory_data <- tt_plot %>%
      dplyr::filter(gene_name %in% mandatory, sig_cat != "NS")
    if (nrow(mandatory_data) > 0) {
      label_genes <- dplyr::bind_rows(label_genes, mandatory_data) %>%
        dplyr::distinct(gene_name, .keep_all = TRUE)
    }
    
    ## Mark mandatory genes for bold styling
    label_genes <- label_genes %>%
      dplyr::mutate(is_highlight = gene_name %in% mandatory)
    
    highlight_data <- label_genes %>% dplyr::filter(is_highlight)
    
    p <- ggplot(tt_plot, aes(x = logFC, y = negLogFDR, color = sig_cat)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(
        values = c("Up" = volcano_colors$up, "Down" = volcano_colors$down, "NS" = volcano_colors$ns),
        name = "Significance"
      ) +
      geom_point(data = label_genes, size = 3) +
      ## Mandatory genes get special highlighting (larger point with yellow fill)
      geom_point(data = highlight_data, size = 5, shape = 21,
                 fill = "yellow", color = "black", stroke = 1.5) +
      ## Labels for standard ORA-selected genes
      geom_text_repel(data = label_genes %>% dplyr::filter(!is_highlight),
                      aes(label = gene_name), color = "black",
                      size = 3.5, box.padding = 0.25, max.overlaps = Inf) +
      ## Bold labels for mandatory/highlighted genes
      geom_text_repel(data = label_genes %>% dplyr::filter(is_highlight),
                      aes(label = gene_name), color = "black",
                      size = 4, box.padding = 0.5, max.overlaps = Inf,
                      fontface = "bold", min.segment.length = 0) +
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
  
  get_padded_limits <- function(values, pad_fraction = 0.02) {
    vals <- values[is.finite(values)]
    if (length(vals) == 0) return(NULL)
    rng <- range(vals)
    pad <- diff(rng) * pad_fraction
    c(rng[1] - pad, rng[2] + pad)
  }
  
  ora_global_ratio_max <- 0
  for (ora_file in ora_files) {
    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) next
    ratios <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    ora_global_ratio_max <- max(ora_global_ratio_max, max(ratios, na.rm = TRUE))
  }
  ora_xlim <- get_padded_limits(c(0, ora_global_ratio_max), pad_fraction = 0.05)
  
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
      dplyr::arrange(p.adjust, pvalue) %>%
      dplyr::slice_head(n = ora_params$top_n)
    
    if (nrow(plot_data) == 0) next
    
    plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    plot_data$neg_log10_padj <- -log10(plot_data$p.adjust)
    plot_data <- plot_data %>%
      dplyr::mutate(pathway_name = str_wrap(pathway_name, width = barplot_params$label_wrap_width),
                    pathway_name = factor(pathway_name, levels = rev(pathway_name)))
    
    direction_label <- if (direction == "Up") "Upregulated" else "Downregulated"
    
    p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = pathway_name, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      scale_fill_gradientn(colors = enrichment_colors$gradient, name = expression(-log[10]*"(FDR)")) +
      scale_x_continuous(limits = ora_xlim, expand = expansion(mult = c(0.02, 0.02))) +
      labs(title = paste0(db, " - ", direction_label, "\n", contrast), x = "Gene Ratio", y = NULL) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.margin = margin(4, 6, 4, 6))
    
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
  
  get_padded_limits <- function(values, pad_fraction = 0.05) {
    vals <- values[is.finite(values)]
    if (length(vals) == 0) return(NULL)
    rng <- range(vals)
    pad <- diff(rng) * pad_fraction
    c(rng[1] - pad, rng[2] + pad)
  }
  
  fgsea_global_abs_nes <- 0
  for (fgsea_file in fgsea_files) {
    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
    if (nrow(fgsea_data) == 0) next
    fgsea_global_abs_nes <- max(fgsea_global_abs_nes, max(abs(fgsea_data$NES), na.rm = TRUE))
  }
  fgsea_xlim <- get_padded_limits(c(-fgsea_global_abs_nes, fgsea_global_abs_nes), pad_fraction = 0.02)
  
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
        pathway_label = str_wrap(pathway_label, width = gsea_plot_params$label_wrap_width),
        pathway_label = factor(pathway_label, levels = pathway_label),
        neg_log10_padj = -log10(pmax(padj, .Machine$double.eps))
      )
    
    if (nrow(plot_data) == 0) next
    
    p <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = neg_log10_padj)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.75) +
      scale_fill_gradientn(colors = enrichment_colors$gradient, name = expression(-log[10]*"(FDR)")) +
      scale_x_continuous(limits = fgsea_xlim, expand = expansion(mult = c(0.01, 0.01))) +
      labs(title = paste0("GSEA ", db, "\n", contrast), x = "NES", y = NULL) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.margin = margin(4, 6, 4, 6)) +
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
  
  get_padded_limits <- function(values, extra_steps = 1) {
    vals <- values[is.finite(values)]
    if (length(vals) == 0) return(NULL)
    pretty_vals <- pretty(range(vals), n = 4)
    step <- if (length(pretty_vals) > 1) diff(pretty_vals)[1] else 1
    c(min(pretty_vals) - extra_steps * step, max(pretty_vals) + extra_steps * step)
  }
  
  dot_global_ratio_max <- 0
  for (ora_file in ora_files) {
    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    if (nrow(ora_data) == 0) next
    ratios <- sapply(strsplit(ora_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
    dot_global_ratio_max <- max(dot_global_ratio_max, max(ratios, na.rm = TRUE))
  }
  dot_xlim <- get_padded_limits(c(0, dot_global_ratio_max), extra_steps = 1)
  
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
      dplyr::arrange(p.adjust, pvalue) %>%
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
      scale_x_continuous(limits = dot_xlim, expand = expansion(mult = c(0.02, 0.02))) +
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
    
    de_files <- .unique_de_files(list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE))
    
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
        anno_list <- list(Genotype = sample_info_tissue$Genotype)
        col_list <- list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02"))
        
        ## Add Sex annotation if available
        if ("Sex" %in% colnames(sample_info_tissue)) {
          anno_list$Sex <- sample_info_tissue$Sex
          col_list$Sex <- c("F" = "#E7298A", "M" = "#7570B3")
        }
        
        col_anno <- HeatmapAnnotation(
          df = as.data.frame(anno_list),
          col = col_list,
          simple_anno_size = unit(4, "mm")
        )
      }
      
      ## Simplify sample names (JS_01 -> 1, JS_20 -> 20)
      simple_names <- gsub("^[A-Za-z]+_?0*", "", colnames(mat))
      colnames(mat) <- simple_names
      
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
    
    de_files <- .unique_de_files(list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE))
    
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
        anno_list <- list(Genotype = sample_info_tissue$Genotype)
        col_list <- list(Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02"))
        
        ## Add Sex annotation if available
        if ("Sex" %in% colnames(sample_info_tissue)) {
          anno_list$Sex <- sample_info_tissue$Sex
          col_list$Sex <- c("F" = "#E7298A", "M" = "#7570B3")
        }
        
        col_anno <- HeatmapAnnotation(
          df = as.data.frame(anno_list),
          col = col_list,
          simple_anno_size = unit(4, "mm")
        )
      }
      
      ## Simplify sample names (JS_01 -> 1, JS_20 -> 20)
      simple_names <- gsub("^[A-Za-z]+_?0*", "", colnames(mat))
      colnames(mat) <- simple_names
      
      ht <- Heatmap(
        mat, col = col_fun, name = "Z-score",
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
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
  
  ## D.3) AVERAGED 4-COLUMN HEATMAPS (genotype x sex group means)
  ## Shows z-scored VST expression averaged per group:
  ##   CTL-F | CTL-M | KAT8KD-F | KAT8KD-M
  ## Same pathway-grouped genes as the individual heatmap (D.2),
  ## but collapsed to 4 columns for a cleaner summary view.
  if (generate_avg4col_heatmaps) {
    cat("\n--- Creating Averaged 4-Column Heatmaps (genotype x sex) ---\n")
    
    ## Reuse gene_categories from D.2 (already built above)
    if (!exists("gene_categories") || length(gene_categories) == 0) {
      cat("[WARN] gene_categories not available, skipping avg4col heatmaps\n")
    } else {
      
      de_files <- .unique_de_files(
        list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE))

      for (de_file in de_files) {
        filename <- basename(de_file)
        contrast <- str_match(filename,
                              "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
        if (is.na(contrast)) next
        
        cat("  Avg 4-col heatmap: ", contrast, "\n", sep = "")
        
        tissue <- sub("_.*", "", contrast)
        sample_info <- vst_data$sample_info
        vst_mat <- vst_data$vst_mat
        
        if ("Depot" %in% colnames(sample_info)) {
          tissue_samples <- rownames(sample_info)[sample_info$Depot == tissue]
        } else {
          tissue_samples <- colnames(vst_mat)[
            grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
        }
        if (length(tissue_samples) == 0) next
        
        vst_tissue <- vst_mat[, colnames(vst_mat) %in% tissue_samples,
                              drop = FALSE]
        sample_info_tissue <- sample_info[colnames(vst_tissue), , drop = FALSE]
        
        ## Need both Genotype and Sex columns
        if (!all(c("Genotype", "Sex") %in% colnames(sample_info_tissue))) {
          cat("[WARN] Need Genotype and Sex columns, skipping\n")
          next
        }
        
        ## Get pathway genes
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
        
        ## Order genes by category
        gene_order <- c()
        for (cat_name in names(gene_categories)) {
          gene_order <- c(gene_order,
                          names(gene_to_cat[gene_to_cat == cat_name]))
        }
        gene_to_cat <- gene_to_cat[gene_order]
        
        ## Build group labels: CTL-F, CTL-M, KAT8KD-F, KAT8KD-M
        group_labels <- paste0(sample_info_tissue$Genotype, "-",
                               sample_info_tissue$Sex)
        group_levels <- c("CTL-F", "CTL-M", "KAT8KD-F", "KAT8KD-M")
        group_labels <- factor(group_labels, levels = group_levels)
        
        ## Average VST expression per group
        vst_sub <- vst_tissue[gene_order, , drop = FALSE]
        avg_mat <- sapply(group_levels, function(grp) {
          cols <- which(group_labels == grp)
          if (length(cols) == 0) return(rep(NA, nrow(vst_sub)))
          if (length(cols) == 1) return(vst_sub[, cols])
          rowMeans(vst_sub[, cols, drop = FALSE])
        })
        rownames(avg_mat) <- gene_order
        colnames(avg_mat) <- group_levels
        
        ## Drop any group with no samples
        keep <- !apply(avg_mat, 2, function(x) all(is.na(x)))
        avg_mat <- avg_mat[, keep, drop = FALSE]
        if (ncol(avg_mat) < 2) next
        
        ## Z-score across the 4 group means (per gene)
        avg_mat <- t(scale(t(avg_mat)))
        z_clip <- heatmap_params$lfc_clip
        avg_mat[avg_mat > z_clip] <- z_clip
        avg_mat[avg_mat < -z_clip] <- -z_clip
        
        col_fun <- colorRamp2(seq(-z_clip, z_clip, length.out = 5),
                              heatmap_params$custom_gradient)
        
        ## Column annotations matching original: Genotype + Sex bars
        col_geno <- sub("-.*", "", colnames(avg_mat))
        col_sex  <- sub(".*-", "", colnames(avg_mat))
        
        col_anno <- HeatmapAnnotation(
          Genotype = col_geno,
          Sex = col_sex,
          col = list(
            Genotype = c("CTL" = "#1B9E77", "KAT8KD" = "#D95F02"),
            Sex = c("F" = "#E7298A", "M" = "#7570B3")
          ),
          simple_anno_size = unit(4, "mm"),
          show_annotation_name = TRUE,
          annotation_name_gp = gpar(fontsize = 8)
        )
        
        ht <- Heatmap(
          avg_mat,
          col = col_fun,
          name = "Z-score",
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = FALSE,
          top_annotation = col_anno,
          row_split = factor(gene_to_cat, levels = names(gene_categories)),
          row_gap = unit(1.5, "mm"),
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 9, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "italic"),
          column_title = paste0(contrast, "\n(averaged per genotype x sex)"),
          column_title_gp = gpar(fontsize = 11, fontface = "bold"),
          border = TRUE,
          border_gp = gpar(col = "black", lwd = 0.8),
          rect_gp = gpar(col = "grey80", lwd = 0.3),
          width = unit(3.5, "cm"),
          heatmap_legend_param = list(
            title = "Z-score",
            title_gp = gpar(fontsize = 9)
          )
        )
        
        n_genes <- length(genes_to_plot)
        plot_file <- paste0("Heatmap_avg4col_", contrast, "_",
                            run_tag, ".png")
        png(file.path(plots_dir, plot_file),
            width = 600,
            height = max(400, n_genes * 18 + 100),
            res = 150)
        draw(ht, padding = unit(c(5, 12, 10, 5), "mm"))
        dev.off()
        cat("[OK] Saved: ", plot_file, "\n", sep = "")
      }
    }
  }
}

## ============================================================
## SECTION E: GENE EXPRESSION BOX PLOTS
## ============================================================
## Per-gene expression for GENES_OF_INTEREST, one panel per gene, CTL vs
## KAT8KD within each depot.
##
## Heatmaps and volcanoes show a gene's effect relative to everything else.
## A box plot shows the actual numbers: where each sample sits, how much the
## groups overlap, and whether one animal is carrying the result. That is what
## you need to look at before putting a gene in a figure or a talk, and it is
## the plot that translates directly to a qPCR follow-up.
##
## Individual sample points are drawn ON TOP of the box, always. With n = 9-10
## per group a box plot alone hides exactly the thing worth seeing -- a
## "difference" produced by one outlier animal.

if (generate_gene_boxplots) {
  cat("\n=== SECTION E: GENE EXPRESSION BOX PLOTS ===\n")

  ## Box plots get their own folder inside the run, alongside plots/panels/.
  box_dir <- file.path(plots_dir, "boxplots")
  dir.create(box_dir, recursive = TRUE, showWarnings = FALSE)

  goi <- if (exists("GENES_OF_INTEREST")) unique(GENES_OF_INTEREST) else character(0)
  if (!length(goi)) {
    cat("[SKIP] GENES_OF_INTEREST is empty (set it in parameters.R)\n")
  } else if (is.null(vst_data)) {
    cat("[SKIP] no VST matrix loaded\n")
  } else {
    vst_mat     <- vst_data$vst_mat
    sample_info <- vst_data$sample_info
    de_files    <- .unique_de_files(list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE))

    present <- intersect(goi, rownames(vst_mat))
    missing <- setdiff(goi, rownames(vst_mat))
    cat("[INFO] ", length(present), " of ", length(goi),
        " genes of interest are in the filtered matrix\n", sep = "")
    if (length(missing))
      cat("[INFO] not measured / filtered out: ", paste(missing, collapse = ", "), "\n", sep = "")

    for (de_file in de_files) {
      contrast <- str_match(basename(de_file), "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (is.na(contrast)) next
      tissue <- sub("_.*$", "", contrast)

      if ("Depot" %in% colnames(sample_info)) {
        keep <- rownames(sample_info)[sample_info$Depot == tissue]
      } else {
        keep <- colnames(vst_mat)[grepl(tissue, colnames(vst_mat), ignore.case = TRUE)]
      }
      keep <- intersect(keep, colnames(vst_mat))
      if (!length(keep) || !length(present)) next

      si_t <- sample_info[keep, , drop = FALSE]
      if (!"Genotype" %in% colnames(si_t)) {
        cat("[SKIP] ", contrast, ": no Genotype column\n", sep = ""); next
      }

      ## Long format: one row per gene per sample.
      bx <- do.call(rbind, lapply(present, function(g) {
        data.frame(gene     = g,
                   Sample   = keep,
                   expr     = as.numeric(vst_mat[g, keep]),
                   Genotype = si_t$Genotype,
                   Sex      = if ("Sex" %in% colnames(si_t)) si_t$Sex else NA_character_,
                   stringsAsFactors = FALSE)
      }))

      ## Annotate each panel with this gene's DE result for THIS depot, so the
      ## picture and the statistic can never disagree.
      tt <- read.csv(de_file, stringsAsFactors = FALSE)
      gcol <- if ("gene_name" %in% names(tt)) "gene_name" else names(tt)[1]
      lab <- vapply(present, function(g) {
        r <- tt[tt[[gcol]] == g, , drop = FALSE]
        if (!nrow(r) || !is.finite(r$padj[1])) return(paste0(g, "  (not tested)"))
        sprintf("%s  (log2FC %+.2f, FDR %s)", g, r$logFC[1],
                ifelse(r$padj[1] < 0.001, format(r$padj[1], digits = 2, scientific = TRUE),
                       sprintf("%.3f", r$padj[1])))
      }, character(1))
      bx$panel <- factor(lab[bx$gene], levels = lab[present])

      gt_levels <- c("CTL", "KAT8KD")
      bx$Genotype <- factor(bx$Genotype, levels = intersect(gt_levels, unique(bx$Genotype)))

      p_box <- ggplot(bx, aes(x = Genotype, y = expr, fill = Genotype)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.55, linewidth = 0.3, colour = "black") +
        geom_jitter(aes(shape = Sex), width = 0.16, height = 0, size = 1.5,
                    alpha = 0.9, colour = "grey20") +
        facet_wrap(~ panel, scales = "free_y",
                   ncol = min(5, max(1, ceiling(sqrt(length(present)))))) +
        scale_fill_manual(values = c(CTL = "#0072B2", KAT8KD = "#D55E00")) +
        labs(title = paste0(tissue, ": genes of interest, CTL vs KAT8KD"),
             subtitle = "VST-normalised expression; every sample is plotted. FDR is this depot's DESeq2 result.",
             x = NULL, y = "VST expression") +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              strip.text = element_text(size = 7.5, face = "bold"),
              legend.position = "bottom")

      ncol_used <- min(5, max(1, ceiling(sqrt(length(present)))))
      nrow_used <- ceiling(length(present) / ncol_used)
      box_file  <- paste0("Boxplot_GOI_", contrast, "_", run_tag, ".png")
      ggsave(file.path(box_dir, box_file), p_box,
             width  = max(6, 2.3 * ncol_used),
             height = max(4, 2.1 * nrow_used),
             dpi = 300, bg = "white", limitsize = FALSE)
      cat("[OK] Saved: ", box_file, " (", length(present), " genes, ",
          length(keep), " samples)\n", sep = "")
    }
  }
}

## ============================================================
## SECTION F: MARKER PANEL FOREST PLOTS
## ============================================================
## log2 fold change with a confidence interval, one row per gene, grouped into
## labelled blocks. Three figure types come out of the same builder:
##   Forest_<panel>_<contrast>      one depot, one panel
##   Combined_<name>_<depot>        two panels side by side, own axes each
##   CrossDepot_<panel>             the same gene in both depots, one row
##
## This is the plot that answers "how big is the effect, and how sure are we"
## at a glance. A volcano ranks thousands of genes and a heatmap shows relative
## pattern; neither says that Pparg fell 0.4 with an interval that just clears
## zero while Slc2a4 fell 0.67 with room to spare. The interval is the point --
## a wide bar crossing zero means "not resolved here", not "no effect".
##
## Panels, gene order, significance tiers and colours all come from
## parameters.R, so the figures can be re-cut without touching this script.

if (generate_marker_forest_plots) {
  cat("\n=== SECTION F: MARKER PANEL FOREST PLOTS ===\n")

  panels <- if (exists("MARKER_PANELS")) MARKER_PANELS else list()
  conf   <- if (exists("marker_panel_conf_level")) marker_panel_conf_level else 0.95
  zmult  <- stats::qnorm(1 - (1 - conf) / 2)
  tiers  <- if (exists("marker_sig_tiers")) marker_sig_tiers else
              list(list(cutoff = 0.001, label = "***"),
                   list(cutoff = 0.01,  label = "**"),
                   list(cutoff = 0.05,  label = "*"))
  dep_col <- if (exists("marker_depot_colours")) marker_depot_colours else
               c(iWAT = "#4FA5D8", gWAT = "#C4749E")

  ## Panel figures get their own folder inside the run, so a plots/ directory
  ## holding 40+ files stays navigable.
  panels_dir <- file.path(plots_dir, "panels")
  dir.create(panels_dir, recursive = TRUE, showWarnings = FALSE)

  ## Stars from the configured tiers. Ordered most- to least-stringent, so the
  ## first cutoff a gene clears is the one it gets.
  .stars <- function(padj) {
    out <- rep("", length(padj))
    for (k in rev(seq_along(tiers))) {
      tk <- tiers[[k]]
      out[!is.na(padj) & padj < tk$cutoff] <- tk$label
    }
    out
  }
  .sig_key <- paste(vapply(rev(tiers), function(t)
                      sprintf("%s FDR<%g", t$label, t$cutoff), character(1)),
                    collapse = "   ")

  ## Build the plotting frame for one panel in one depot.
  .panel_data <- function(tt, gcol, spec) {
    rows <- list(); absent <- character(0)
    for (gname in names(spec$groups)) {
      for (g in spec$groups[[gname]]) {
        r <- tt[tt[[gcol]] == g, , drop = FALSE]
        if (!nrow(r) || !is.finite(r$logFC[1]) || !is.finite(r$lfcSE[1])) {
          absent <- c(absent, g); next
        }
        rows[[length(rows) + 1]] <- data.frame(
          gene = g, group = gname, logFC = r$logFC[1], lfcSE = r$lfcSE[1],
          padj = ifelse(is.finite(r$padj[1]), r$padj[1], NA_real_),
          stringsAsFactors = FALSE)
      }
    }
    if (!length(rows)) return(NULL)
    fp <- do.call(rbind, rows)
    fp$lo <- fp$logFC - zmult * fp$lfcSE
    fp$hi <- fp$logFC + zmult * fp$lfcSE
    fp$negLogFDR <- ifelse(is.na(fp$padj), 0, -log10(pmax(fp$padj, .Machine$double.xmin)))
    fp$stars <- .stars(fp$padj)
    attr(fp, "absent") <- absent
    fp
  }

  ## One forest panel as a ggplot object, so single figures and the
  ## side-by-side combined figures are guaranteed to look identical.
  .forest_plot <- function(fp, spec, title, thr, compact = FALSE) {
    absent <- attr(fp, "absent")
    fp$gene  <- factor(fp$gene,  levels = rev(unique(fp$gene)))
    fp$group <- factor(as.character(fp$group), levels = names(spec$groups))

    ## Significance counted on FDR ALONE. The |log2FC| gate is a rule for
    ## genome-wide DEG lists; on a hand-picked panel it would hide exactly the
    ## modest real effects this plot exists to show. The gate is still DRAWN
    ## (dotted) so both criteria stay visible.
    n_sig  <- sum(!is.na(fp$padj) & fp$padj < thr$fdr_cut)
    n_up   <- sum(!is.na(fp$padj) & fp$padj < thr$fdr_cut & fp$logFC > 0)
    cap <- sprintf(
      "Tissue thresholds: FDR<%s, |log2FC|>%s (dotted) | n=%d | Up: %d | Down: %d | ns: %d\n%s   |   bars = %d%% CI",
      thr$fdr_cut, thr$logFC_cut, nrow(fp), n_up, n_sig - n_up, nrow(fp) - n_sig,
      .sig_key, round(100 * conf))
    if (length(absent))
      cap <- paste0(cap, "\nnot measured: ", paste(absent, collapse = ", "))

    rng   <- max(abs(c(fp$lo, fp$hi, thr$logFC_cut)), na.rm = TRUE) * 1.08
    nudge <- rng * 0.035

    ggplot(fp, aes(x = logFC, y = gene)) +
      geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
      geom_vline(xintercept = c(-thr$logFC_cut, thr$logFC_cut),
                 colour = "grey75", linetype = "dotted", linewidth = 0.4) +
      geom_errorbar(aes(xmin = lo, xmax = hi, colour = negLogFDR),
                    orientation = "y", width = 0, linewidth = 0.7) +
      geom_point(aes(colour = negLogFDR), size = if (compact) 2.6 else 3.2) +
      geom_text(aes(x = hi + nudge, label = stars), hjust = 0, vjust = 0.78,
                size = if (compact) 3.6 else 4.2, fontface = "bold", na.rm = TRUE) +
      scale_colour_gradientn(
        colours = c("grey78", "#F5E3C8", "#FDBB6E", "#E8850C", "#B33000"),
        name = expression(-log[10]~FDR), limits = c(0, NA)) +
      facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
      scale_x_continuous(limits = c(-rng, rng * 1.1)) +
      labs(title = title, subtitle = spec$subtitle,
           x = expression(log[2]~"Fold Change (KAT8KD/CTL)"), y = NULL,
           caption = cap) +
      theme_bw(base_size = if (compact) 10 else 12) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank(),
            axis.text.y  = element_text(face = "italic", size = if (compact) 9 else 11),
            strip.placement = "outside",
            strip.text.y.left = element_text(angle = 0, face = "bold"),
            strip.background = element_rect(fill = "grey95", colour = "grey70"),
            plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
            plot.caption  = element_text(hjust = 0.5, colour = "grey45",
                                         size = if (compact) 7 else 8))
  }

  if (!length(panels)) {
    cat("[SKIP] MARKER_PANELS not defined in parameters.R\n")
  } else {
    de_files <- .unique_de_files(
      .unique_de_files(list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)))

    cross <- list()            ## every depot's rows, for the cross-depot figures
    built <- list()            ## panel plots per depot, for the combined figures

    for (de_file in de_files) {
      contrast <- str_match(basename(de_file), "DE_tissue_(.+)_[0-9]+_[0-9]+\\.csv")[2]
      if (is.na(contrast)) next
      tissue <- sub("_.*$", "", contrast)

      tt   <- read.csv(de_file, stringsAsFactors = FALSE)
      gcol <- if ("gene_name" %in% names(tt)) "gene_name" else names(tt)[1]
      if (!"lfcSE" %in% names(tt)) {
        cat("[SKIP] ", contrast, ": DE table has no lfcSE column, cannot draw CIs\n", sep = "")
        next
      }
      thr <- get_tissue_thresholds(contrast)

      for (pname in names(panels)) {
        spec <- panels[[pname]]
        fp   <- .panel_data(tt, gcol, spec)
        if (is.null(fp)) {
          cat("[SKIP] ", tissue, " / ", pname, ": no genes measured\n", sep = ""); next
        }

        p_f <- .forest_plot(fp, spec, paste0(tissue, ": ", pname), thr)
        safe <- gsub("[^A-Za-z0-9]+", "_", pname)
        fp_file <- paste0("Forest_", safe, "_", contrast, "_", run_tag, ".png")
        ggsave(file.path(panels_dir, fp_file), p_f,
               width = 9, height = max(4, 0.42 * nrow(fp) + 2.4),
               dpi = 300, bg = "white", limitsize = FALSE)
        n_sig <- sum(!is.na(fp$padj) & fp$padj < thr$fdr_cut)
        cat("[OK] ", fp_file, "  (", nrow(fp), " genes, ", n_sig, " significant)\n", sep = "")

        built[[paste(tissue, pname, sep = "||")]] <-
          list(plot = .forest_plot(fp, spec, paste0(tissue, ": ", pname), thr,
                                   compact = TRUE), n = nrow(fp))
        fp$Depot <- tissue
        cross[[length(cross) + 1]] <- transform(fp, panel_name = pname)
      }

      ## ---- side-by-side combined figures, per depot --------------------
      combos <- if (exists("COMBINED_PANELS")) COMBINED_PANELS else list()
      for (cname in names(combos)) {
        want <- combos[[cname]]
        keys <- paste(tissue, want, sep = "||")
        if (!all(keys %in% names(built))) next
        pl <- lapply(keys, function(k) built[[k]]$plot)
        nn <- max(vapply(keys, function(k) built[[k]]$n, numeric(1)))
        cf <- paste0("Combined_", gsub("[^A-Za-z0-9]+", "_", cname),
                     "_", contrast, "_", run_tag, ".png")
        ok <- .save_side_by_side(
          pl, file.path(panels_dir, cf),
          title = paste0(tissue, ": ", cname),
          subtitle = sprintf("KAT8KD vs CTL | tissue thresholds (FDR<%s, |logFC|>%s)",
                             thr$fdr_cut, thr$logFC_cut),
          width = 15, height = max(5.5, 0.40 * nn + 2.6))
        if (ok) cat("[OK] ", cf, "\n", sep = "")
      }
    }

    ## ---- cross-depot: the same gene in both depots, on one row ----------
    ## Two figures side by side make that comparison from memory. One row makes
    ## it a direct read: overlapping intervals mean the depots are NOT resolved
    ## apart, however different the two panels looked.
    want_cd <- if (exists("CROSS_DEPOT_PANELS")) CROSS_DEPOT_PANELS else names(panels)
    cd <- if (length(cross)) dplyr::bind_rows(cross) else NULL
    if (!is.null(cd) && dplyr::n_distinct(cd$Depot) > 1 && length(want_cd)) {
      for (pname in intersect(want_cd, unique(cd$panel_name))) {
        sub <- cd[cd$panel_name == pname, , drop = FALSE]
        ## Genes must be measured in EVERY depot: a missing row would read as
        ## "no effect" rather than "not measured".
        keep_g  <- names(which(table(sub$gene) == dplyr::n_distinct(sub$Depot)))
        dropped <- setdiff(as.character(unique(sub$gene)), keep_g)
        sub <- sub[sub$gene %in% keep_g, , drop = FALSE]
        if (!nrow(sub)) next

        ord <- unlist(panels[[pname]]$groups, use.names = FALSE)
        ord <- ord[ord %in% keep_g]
        sub$gene  <- factor(as.character(sub$gene), levels = rev(ord))
        sub$Depot <- factor(sub$Depot, levels = intersect(names(dep_col), unique(sub$Depot)))

        lfc_cut <- get_tissue_thresholds("iWAT")$logFC_cut
        rng2 <- max(abs(c(sub$lo, sub$hi, lfc_cut)), na.rm = TRUE) * 1.08
        dg   <- ggplot2::position_dodge(width = 0.62)
        cap2 <- paste0(.sig_key, "   (per depot)   |   bars = ", round(100 * conf), "% CI",
                       if (length(dropped))
                         paste0("\nnot measured in both depots, so excluded: ",
                                paste(dropped, collapse = ", ")) else "")

        p_c <- ggplot(sub, aes(x = logFC, y = gene, colour = Depot,
                               shape = Depot, group = Depot)) +
          geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
          geom_vline(xintercept = c(-lfc_cut, lfc_cut), colour = "grey75",
                     linetype = "dotted", linewidth = 0.4) +
          geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                        width = 0, linewidth = 0.7, position = dg) +
          geom_point(size = 3, position = dg) +
          geom_text(aes(x = hi + rng2 * 0.03, label = stars), position = dg,
                    hjust = 0, vjust = 0.78, size = 3.8, fontface = "bold",
                    show.legend = FALSE, na.rm = TRUE) +
          scale_colour_manual(values = dep_col) +
          scale_shape_manual(values = c(iWAT = 16, gWAT = 17)) +
          scale_x_continuous(limits = c(-rng2, rng2 * 1.12)) +
          labs(title = paste0(pname, ": iWAT vs gWAT KAT8KD"),
               subtitle = sprintf("Cross-depot consistency check | tissue thresholds (FDR<%s, |logFC|>%s)",
                                  get_tissue_thresholds("iWAT")$fdr_cut, lfc_cut),
               x = expression(log[2]~"Fold Change (KAT8KD/CTL)"), y = NULL,
               caption = cap2) +
          theme_bw(base_size = 12) +
          theme(panel.grid.major.y = element_line(colour = "grey93"),
                panel.grid.minor   = element_blank(),
                axis.text.y  = element_text(face = "italic", size = 11),
                plot.title    = element_text(face = "bold", hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
                plot.caption  = element_text(hjust = 0.5, colour = "grey45", size = 8),
                legend.position = "right")

        safe <- gsub("[^A-Za-z0-9]+", "_", pname)
        cf <- paste0("CrossDepot_", safe, "_iWAT_vs_gWAT_", run_tag, ".png")
        ggsave(file.path(panels_dir, cf), p_c,
               width = 9.5, height = max(4.5, 0.46 * length(keep_g) + 2.6),
               dpi = 300, bg = "white", limitsize = FALSE)
        cat("[OK] ", cf, "  (", length(keep_g), " genes in both depots)\n", sep = "")
      }
    }
    cat("[INFO] panel figures -> ", panels_dir, "\n", sep = "")
  }
}

## ============================================================
## COMPLETE
## ============================================================

cat("\n======================================\n")
cat("=== PART 4 (ALL VISUALIZATIONS) COMPLETE ===\n")
cat("======================================\n")
cat("All plots saved to: ", plots_dir, "\n\n", sep = "")