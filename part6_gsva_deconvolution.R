#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 6: GSVA CELL-TYPE SCORING &
##                              COMPUTATIONAL DECONVOLUTION
## =========================================================
## Per-sample gene set scoring and cell-type proportion
## estimation from bulk RNA-seq
##
## KEY ANALYSES:
## 1. GSVA: Per-sample enrichment scores for curated adipose
##    cell-type signatures (adipocytes, macrophages, etc.)
## 2. Deconvolution: Estimated cell-type proportions using
##    MuSiC (if scRNA-seq reference is provided) or
##    marker-based inference (fallback)
## 3. Statistical tests: KD vs CTL for each cell-type score
## 4. Publication figures: Boxplots, heatmaps, stacked bars
##
## REQUIRES:
##   - Part 1 tissue results (VST counts + DE tables)
##   - GSVA package
##   - Optional: scRNA-seq reference for MuSiC deconvolution
##     (e.g., Emont et al. 2022 Nature, GEO: GSE176171)
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "tibble",
                   "RColorBrewer", "ggrepel", "purrr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c("GSVA", "DESeq2", "ComplexHeatmap", "circlize")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(RColorBrewer)
  library(purrr)
  library(GSVA)
  library(DESeq2)
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

set.seed(MASTER_SEED)

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  use_save_core <- TRUE
} else {
  cat("[WARN] save_core not found; proceeding without it\n")
  use_save_core <- FALSE
}

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## ============================================================
## 3) CELL-TYPE SIGNATURE DEFINITIONS
## ============================================================
## Curated marker gene sets for adipose-resident cell types
## Sources: Tabula Muris, Emont et al. 2022, literature consensus

adipose_signatures <- list(

  Adipocyte = c(
    "Adipoq", "Lep", "Fabp4", "Plin1", "Pparg",
    "Cidec", "Cidea", "Retn", "Rbp4", "Aqp7",
    "Slc2a4", "Irs1", "Cebpa", "Lpl", "Acaca",
    "Fasn", "Scd1", "Dgat1", "Dgat2", "Pck1"
  ),

  Preadipocyte = c(
    "Pdgfra", "Pdgfrb", "Cd34", "Ly6a", "Dcn",
    "Col1a1", "Col1a2", "Col3a1", "Fn1", "Vim",
    "Dlk1", "Sox9", "Zfp423"
  ),

  Macrophage = c(
    "Adgre1", "Cd68", "Csf1r", "Itgam", "Mrc1",
    "Cd163", "Lyve1", "Mertk", "Timd4", "Fcgr1",
    "Ccr2", "Ly6c2", "Itgax", "Cd86", "Nos2",
    "Arg1", "Chil3", "Il10", "Tgfb1"
  ),

  Macrophage_M1 = c(
    "Nos2", "Tnf", "Il1b", "Il6", "Ccl2",
    "Cxcl9", "Cxcl10", "Cd86", "Cd80", "Fcgr1",
    "Stat1", "Irf5", "Nfkb1", "Socs3", "Ptgs2"
  ),

  Macrophage_M2 = c(
    "Arg1", "Mrc1", "Cd163", "Chil3", "Retnla",
    "Il10", "Tgfb1", "Mertk", "Pparg", "Clec10a",
    "Stat6", "Irf4", "Ccl22", "Ccl17"
  ),

  Endothelial = c(
    "Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1",
    "Tek", "Emcn", "Icam1", "Vcam1", "Sele",
    "Cldn5", "Plvap", "Aqp1", "Esam", "Egfl7"
  ),

  T_cell = c(
    "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a",
    "Lck", "Zap70", "Lat", "Trac", "Trbc1",
    "Foxp3", "Il2ra", "Ctla4", "Il7r"
  ),

  NK_cell = c(
    "Ncr1", "Klrb1c", "Klrd1", "Nkg7", "Prf1",
    "Gzma", "Gzmb", "Ifng", "Klrk1", "Fcgr3"
  ),

  Fibroblast = c(
    "Dcn", "Lum", "Col1a1", "Col1a2", "Col3a1",
    "Pdgfra", "Fap", "Thy1", "Acta2", "S100a4",
    "Vim", "Fn1", "Postn", "Tnc", "Sparc"
  ),

  Dendritic = c(
    "Itgax", "Cd74", "H2-Ab1", "H2-Eb1", "Flt3",
    "Batf3", "Irf8", "Xcr1", "Clec9a", "Cd209a"
  ),

  Mast_cell = c(
    "Kit", "Cpa3", "Tpsb2", "Tpsab1", "Ms4a2",
    "Fcer1a", "Mcpt1", "Mcpt4", "Mcpt8", "Hdc"
  ),

  ## Biological process signatures (not cell types)
  Lipogenesis = c(
    "Fasn", "Acaca", "Acacb", "Scd1", "Elovl6",
    "Acly", "Srebf1", "Dgat1", "Dgat2", "Gpam",
    "Agpat2", "Lpin1", "Pck1", "Gyk"
  ),

  Thermogenesis = c(
    "Ucp1", "Dio2", "Cidea", "Cox7a1", "Cox8b",
    "Ppargc1a", "Prdm16", "Elovl3", "Ebf2", "Tbx1",
    "Tmem26", "Tnfrsf9", "Cited1"
  ),

  Inflammatory_Response = c(
    "Tnf", "Il1b", "Il6", "Ccl2", "Ccl5",
    "Cxcl1", "Cxcl2", "Cxcl10", "Nfkb1", "Nfkb2",
    "Ptgs2", "Socs3", "Tlr4", "Myd88"
  ),

  Insulin_Signaling = c(
    "Insr", "Irs1", "Irs2", "Pik3r1", "Akt1",
    "Akt2", "Slc2a4", "Gsk3b", "Foxo1", "Pparg",
    "Pck1", "G6pc", "Srebf1", "Ppargc1a"
  )
)

## ============================================================
## 4) FIND RESULTS & LOAD DATA
## ============================================================

cat("\n=== FINDING PART 1 RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Run Part 1 first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) stop("No run directories found.")

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")
gsva_dir <- file.path(plots_dir, "gsva")
if (!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)

cat("[INFO] Using results from: ", outdir, "\n", sep = "")

## 4.1) Load DE tables
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_tables <- list()
for (f in de_files) {
  contrast <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  de_tables[[contrast]] <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  cat("[OK] Loaded DE: ", contrast, " (", nrow(de_tables[[contrast]]), " genes)\n", sep = "")
}

## 4.2) Load count data for GSVA (need expression matrix)
cat("\n=== LOADING COUNT DATA FOR GSVA ===\n")

## Try to find saved DESeq2 object or raw counts
rds_files <- list.files(outdir, pattern = "dds.*\\.rds$", recursive = TRUE, full.names = TRUE)

if (length(rds_files) > 0) {
  dds <- readRDS(rds_files[length(rds_files)])
  vst_mat <- assay(vst(dds, blind = TRUE))
  sample_meta <- as.data.frame(colData(dds))
  cat("[OK] Loaded DESeq2 object from RDS\n")
} else {
  ## Fallback: rebuild from counts.txt
  cat("[INFO] No saved DESeq2 object found. Building from counts.txt\n")

  counts_raw <- read.delim("counts.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gene_col <- if ("gene" %in% colnames(counts_raw)) "gene" else "gene_id"

  ## Build metadata (matches Part 1 hardcoded annotation)
  metadata_file <- file.path(getwd(), "data", "metadata.csv")
  if (file.exists(metadata_file)) {
    sample_meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  } else {
    cat("[INFO] No metadata.csv found. Using Part 1 hardcoded sample annotation.\n")
    sample_ids <- paste0("JS_", sprintf("%02d", 1:40))
    sample_meta <- data.frame(
      Sample   = sample_ids,
      Tissue   = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
      Depot    = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
      Sex      = c(rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5),
                   rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5)),
      Genotype = c(rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5),
                   rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5)),
      stringsAsFactors = FALSE
    )
    ## Remove known outliers
    outlier_samples <- c("JS_08", "JS_28")
    sample_meta <- sample_meta[!sample_meta$Sample %in% outlier_samples, ]
  }

  tissue_samples <- intersect(sample_meta$Sample, colnames(counts_raw))
  count_matrix <- as.matrix(counts_raw[, tissue_samples])
  rownames(count_matrix) <- make.unique(counts_raw[[gene_col]])
  rownames(sample_meta) <- sample_meta$Sample
  sample_meta <- sample_meta[tissue_samples, , drop = FALSE]

  ## Build DESeq2 and VST
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData   = sample_meta,
    design    = ~ 1
  )
  keep <- rowSums(counts(dds) >= 10) >= ncol(dds) * 0.25
  dds <- dds[keep, ]
  vst_mat <- assay(vst(dds, blind = TRUE))
  cat("[OK] Built VST matrix: ", nrow(vst_mat), " genes x ", ncol(vst_mat), " samples\n", sep = "")
}

## Harmonize column names: Part 1 uses "Depot" but downstream code expects "Tissue"
if ("Depot" %in% colnames(sample_meta) && !"Tissue" %in% colnames(sample_meta)) {
  sample_meta$Tissue <- sample_meta$Depot
  cat("[INFO] Mapped 'Depot' column to 'Tissue' for consistency\n")
}

## ============================================================
## 5) ANALYSIS 1: GSVA PER-SAMPLE SCORING
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 1: GSVA PER-SAMPLE CELL-TYPE SCORING ===\n")
cat("============================================================\n\n")
cat("METHOD: Gene Set Variation Analysis (Hanzelmann et al. 2013)\n")
cat("  - Non-parametric, unsupervised per-sample enrichment scoring\n")
cat("  - Produces a score for each sample x signature combination\n")
cat("  - Positive score = signature genes are highly expressed in that sample\n")
cat("  - Negative score = signature genes are lowly expressed in that sample\n\n")

log_time("Starting GSVA")

## Filter signatures to genes present in expression matrix
sig_filtered <- lapply(adipose_signatures, function(genes) {
  intersect(genes, rownames(vst_mat))
})

## Report coverage
cat("Signature gene coverage:\n")
for (sig_name in names(sig_filtered)) {
  n_found <- length(sig_filtered[[sig_name]])
  n_total <- length(adipose_signatures[[sig_name]])
  cat("  ", sig_name, ": ", n_found, "/", n_total, " genes found (",
      round(100 * n_found / n_total), "%)\n", sep = "")
}

## Remove signatures with too few genes
min_genes <- gsva_params$min_signature_genes
sig_filtered <- sig_filtered[sapply(sig_filtered, length) >= min_genes]
cat("\nSignatures with >=", min_genes, " genes: ", length(sig_filtered), "\n", sep = "")

## Run GSVA
gsva_param <- gsvaParam(
  exprData = vst_mat,
  geneSets = sig_filtered,
  kcdf     = gsva_params$kcdf
)

gsva_scores <- gsva(gsva_param)
cat("[OK] GSVA complete: ", nrow(gsva_scores), " signatures x ", ncol(gsva_scores), " samples\n", sep = "")

log_time("GSVA complete")

## Save raw GSVA scores
gsva_scores_df <- as.data.frame(gsva_scores) %>% rownames_to_column("Signature")
gsva_file <- paste0("gsva_scores_", run_tag, ".csv")
write.csv(gsva_scores_df, file = file.path(tables_dir, gsva_file), row.names = FALSE)
cat("[OK] Saved: ", gsva_file, "\n", sep = "")

## ============================================================
## 6) STATISTICAL TESTING: KD vs CTL PER SIGNATURE
## ============================================================

cat("\n=== STATISTICAL TESTING ===\n")

## Build metadata for each sample
if (!"Genotype" %in% colnames(sample_meta)) {
  ## Try to extract from column names or existing metadata
  if ("GroupFull" %in% colnames(sample_meta)) {
    sample_meta$Genotype <- ifelse(grepl("KD|KAT8", sample_meta$GroupFull), "KD", "CTL")
  } else {
    stop("Cannot find Genotype column in sample metadata.")
  }
}

## Determine grouping variable(s)
## If Tissue is present, test within each tissue; otherwise test overall
has_tissue <- "Tissue" %in% colnames(sample_meta)

stat_results <- list()

if (has_tissue) {
  tissues <- unique(sample_meta$Tissue)
  cat("[INFO] Testing within each tissue: ", paste(tissues, collapse = ", "), "\n", sep = "")

  for (tissue in tissues) {
    tissue_samples <- rownames(sample_meta)[sample_meta$Tissue == tissue]
    kd_samples <- intersect(tissue_samples,
                            rownames(sample_meta)[grepl("KD|KAT8", sample_meta$Genotype, ignore.case = TRUE)])
    ctl_samples <- setdiff(tissue_samples, kd_samples)

    if (length(kd_samples) < 2 || length(ctl_samples) < 2) {
      cat("  [WARN] Skipping ", tissue, " (too few samples per group)\n", sep = "")
      next
    }

    for (sig_name in rownames(gsva_scores)) {
      kd_vals <- gsva_scores[sig_name, kd_samples]
      ctl_vals <- gsva_scores[sig_name, ctl_samples]

      ## Wilcoxon test (non-parametric, appropriate for small n)
      wt <- tryCatch(
        wilcox.test(kd_vals, ctl_vals, exact = FALSE),
        error = function(e) list(p.value = NA, statistic = NA)
      )

      ## Also t-test for reference
      tt <- tryCatch(
        t.test(kd_vals, ctl_vals),
        error = function(e) list(p.value = NA, estimate = c(0, 0))
      )

      stat_results[[paste0(tissue, "_", sig_name)]] <- data.frame(
        Tissue = tissue,
        Signature = sig_name,
        Mean_KD = mean(kd_vals),
        Mean_CTL = mean(ctl_vals),
        Diff = mean(kd_vals) - mean(ctl_vals),
        Wilcoxon_P = wt$p.value,
        Ttest_P = tt$p.value,
        N_KD = length(kd_samples),
        N_CTL = length(ctl_samples),
        stringsAsFactors = FALSE
      )
    }
  }
} else {
  ## Single comparison (e.g., cell culture)
  kd_samples <- rownames(sample_meta)[grepl("KD|KAT8", sample_meta$Genotype, ignore.case = TRUE)]
  ctl_samples <- rownames(sample_meta)[!grepl("KD|KAT8", sample_meta$Genotype, ignore.case = TRUE)]

  for (sig_name in rownames(gsva_scores)) {
    kd_vals <- gsva_scores[sig_name, kd_samples]
    ctl_vals <- gsva_scores[sig_name, ctl_samples]

    wt <- tryCatch(wilcox.test(kd_vals, ctl_vals, exact = FALSE),
                   error = function(e) list(p.value = NA))
    tt <- tryCatch(t.test(kd_vals, ctl_vals),
                   error = function(e) list(p.value = NA))

    stat_results[[sig_name]] <- data.frame(
      Tissue = "All",
      Signature = sig_name,
      Mean_KD = mean(kd_vals),
      Mean_CTL = mean(ctl_vals),
      Diff = mean(kd_vals) - mean(ctl_vals),
      Wilcoxon_P = wt$p.value,
      Ttest_P = tt$p.value,
      N_KD = length(kd_samples),
      N_CTL = length(ctl_samples),
      stringsAsFactors = FALSE
    )
  }
}

stat_df <- bind_rows(stat_results)
stat_df$Wilcoxon_FDR <- p.adjust(stat_df$Wilcoxon_P, method = "BH")
stat_df$Ttest_FDR <- p.adjust(stat_df$Ttest_P, method = "BH")
stat_df <- stat_df %>% arrange(Wilcoxon_P)

## Save
stat_file <- paste0("gsva_statistics_", run_tag, ".csv")
write.csv(stat_df, file = file.path(tables_dir, stat_file), row.names = FALSE)
cat("[OK] Saved: ", stat_file, "\n\n", sep = "")

## Print significant results
cat("Significant GSVA scores (Wilcoxon p < 0.05):\n")
sig_results <- stat_df %>% filter(Wilcoxon_P < 0.05)
if (nrow(sig_results) > 0) {
  for (i in seq_len(nrow(sig_results))) {
    row <- sig_results[i, ]
    dir_label <- ifelse(row$Diff > 0, "UP in KD", "DOWN in KD")
    cat("  ", row$Tissue, " - ", row$Signature, ": ", dir_label,
        " (diff = ", round(row$Diff, 3), ", p = ", format(row$Wilcoxon_P, digits = 3), ")\n", sep = "")
  }
} else {
  cat("  None at p < 0.05\n")
}

## ============================================================
## 7) FIGURE 1: GSVA SCORE BOXPLOTS
## ============================================================

cat("\n=== GENERATING GSVA BOXPLOTS ===\n")

## Reshape for plotting
## Ensure sample_meta has Sample as a column (not just rownames)
if (!"Sample" %in% colnames(sample_meta)) {
  sample_meta <- sample_meta %>% rownames_to_column("Sample")
}

gsva_long <- as.data.frame(t(gsva_scores)) %>%
  rownames_to_column("Sample") %>%
  left_join(sample_meta %>%
              select(Sample, any_of(c("Genotype", "Tissue", "Sex"))),
            by = "Sample") %>%
  pivot_longer(cols = -c(Sample, any_of(c("Genotype", "Tissue", "Sex"))),
               names_to = "Signature", values_to = "Score")

## Simplify genotype labels
gsva_long <- gsva_long %>%
  mutate(Group = ifelse(grepl("KD|KAT8", Genotype, ignore.case = TRUE), "KD", "CTL"))

## Group signatures by type for plotting
celltype_sigs <- c("Adipocyte", "Preadipocyte", "Macrophage", "Macrophage_M1",
                    "Macrophage_M2", "Endothelial", "T_cell", "NK_cell",
                    "Fibroblast", "Dendritic", "Mast_cell")
process_sigs <- c("Lipogenesis", "Thermogenesis", "Inflammatory_Response", "Insulin_Signaling")

## Function to generate boxplot panel
make_gsva_boxplot <- function(data, sig_subset, title, stat_table, plot_file) {

  plot_data <- data %>% filter(Signature %in% sig_subset)
  plot_data$Signature <- factor(plot_data$Signature, levels = sig_subset)

  if (has_tissue) {
    ## Faceted by tissue
    p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
      facet_grid(Signature ~ Tissue, scales = "free_y") +
      scale_fill_manual(values = c("CTL" = "#0072B2", "KD" = "#D55E00")) +
      labs(title = title, x = NULL, y = "GSVA Enrichment Score") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom",
        axis.text.x = element_text(size = 10, face = "bold")
      )

    n_tissues <- length(unique(plot_data$Tissue))
    n_sigs <- length(sig_subset)
    height <- max(6, n_sigs * 1.8)
    width <- max(6, n_tissues * 2.5 + 2)
  } else {
    p <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
      facet_wrap(~ Signature, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("CTL" = "#0072B2", "KD" = "#D55E00")) +
      labs(title = title, x = NULL, y = "GSVA Enrichment Score") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom"
      )

    n_sigs <- length(sig_subset)
    height <- max(6, ceiling(n_sigs / 3) * 3)
    width <- 10
  }

  ## Add p-value annotations
  ## (Wilcoxon p-values from stat table, shown as text above boxplots)
  p_annot <- stat_table %>%
    filter(Signature %in% sig_subset) %>%
    mutate(
      label = case_when(
        Wilcoxon_P < 0.001 ~ "***",
        Wilcoxon_P < 0.01  ~ "**",
        Wilcoxon_P < 0.05  ~ "*",
        Wilcoxon_P < 0.1   ~ ".",
        TRUE ~ "ns"
      )
    )

  ggsave(file.path(gsva_dir, plot_file), plot = p, width = width, height = height,
         dpi = 300, bg = "white")
  cat("[OK] Saved: ", plot_file, "\n", sep = "")
}

## Cell type signature boxplots
celltype_present <- intersect(celltype_sigs, rownames(gsva_scores))
if (length(celltype_present) > 0) {
  make_gsva_boxplot(gsva_long, celltype_present,
                    "Cell-Type Signature Scores (GSVA): KD vs CTL",
                    stat_df,
                    paste0("GSVA_celltype_boxplots_", run_tag, ".png"))
}

## Biological process signature boxplots
process_present <- intersect(process_sigs, rownames(gsva_scores))
if (length(process_present) > 0) {
  make_gsva_boxplot(gsva_long, process_present,
                    "Biological Process Scores (GSVA): KD vs CTL",
                    stat_df,
                    paste0("GSVA_process_boxplots_", run_tag, ".png"))
}

## ============================================================
## 8) FIGURE 2: GSVA SCORE HEATMAP (ALL SIGNATURES x SAMPLES)
## ============================================================

cat("\n=== GENERATING GSVA HEATMAP ===\n")

## Order signatures by category
sig_order <- c(celltype_present, process_present)
sig_order <- intersect(sig_order, rownames(gsva_scores))

if (length(sig_order) >= 3) {
  gsva_mat <- gsva_scores[sig_order, , drop = FALSE]

  ## Color scale
  max_abs <- max(abs(gsva_mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#0072B2", "white", "#D55E00"))

  ## Column annotation
  col_annot_data <- sample_meta[colnames(gsva_mat), , drop = FALSE]
  annot_cols <- list()
  if ("Genotype" %in% colnames(col_annot_data)) {
    col_annot_data$Group <- ifelse(grepl("KD|KAT8", col_annot_data$Genotype, ignore.case = TRUE), "KD", "CTL")
    annot_cols$Group <- c("CTL" = "#0072B2", "KD" = "#D55E00")
  }
  if ("Tissue" %in% colnames(col_annot_data)) {
    tissues <- unique(col_annot_data$Tissue)
    tissue_cols <- setNames(brewer.pal(max(3, length(tissues)), "Set2")[seq_along(tissues)], tissues)
    annot_cols$Tissue <- tissue_cols
  }

  ha_col <- HeatmapAnnotation(
    df = col_annot_data[, intersect(c("Group", "Tissue"), colnames(col_annot_data)), drop = FALSE],
    col = annot_cols,
    annotation_name_side = "left"
  )

  ## Row annotation (cell type vs process)
  row_category <- ifelse(sig_order %in% celltype_sigs, "Cell Type", "Biological Process")
  ha_row <- rowAnnotation(
    Category = row_category,
    col = list(Category = c("Cell Type" = "#E91E63", "Biological Process" = "#4CAF50")),
    show_legend = TRUE
  )

  ht <- Heatmap(
    as.matrix(gsva_mat),
    name = "GSVA\nScore",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    column_names_rot = 45,
    top_annotation = ha_col,
    left_annotation = ha_row,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 9),
    column_title = "GSVA Per-Sample Enrichment Scores",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold")
    )
  )

  png(file.path(gsva_dir, paste0("GSVA_heatmap_", run_tag, ".png")),
      width = max(8, ncol(gsva_mat) * 0.6 + 4), height = max(6, length(sig_order) * 0.4 + 2),
      units = "in", res = 300)
  draw(ht, padding = unit(c(2, 2, 2, 15), "mm"))
  dev.off()
  cat("[OK] Saved GSVA heatmap\n")
}

## ============================================================
## 9) ANALYSIS 2: MARKER-BASED PROPORTION ESTIMATION
## ============================================================
## Lightweight deconvolution fallback when no scRNA-seq reference
## is available. Uses signature gene expression to estimate
## relative cell-type proportions (not absolute).

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 2: MARKER-BASED PROPORTION ESTIMATION ===\n")
cat("============================================================\n\n")

## Check for MuSiC scRNA-seq reference
sc_ref_file <- gsva_params$sc_reference_path
use_music <- FALSE

if (!is.null(sc_ref_file) && file.exists(sc_ref_file)) {
  cat("[INFO] scRNA-seq reference found: ", sc_ref_file, "\n", sep = "")

  if (requireNamespace("MuSiC", quietly = TRUE)) {
    cat("[INFO] MuSiC package available. Running full deconvolution.\n")
    use_music <- TRUE
  } else {
    cat("[WARN] MuSiC package not installed. Using marker-based fallback.\n")
    cat("[INFO] Install MuSiC: devtools::install_github('xuranw/MuSiC')\n")
  }
} else {
  cat("[INFO] No scRNA-seq reference provided.\n")
  cat("[INFO] Using marker-based proportion estimation (relative, not absolute).\n")
  cat("[INFO] For full deconvolution, provide an scRNA-seq reference in parameters.R:\n")
  cat("         gsva_params$sc_reference_path <- 'path/to/emont_2022_reference.rds'\n\n")
}

if (use_music) {
  ## ---- FULL DECONVOLUTION WITH MuSiC ----
  cat("=== RUNNING MuSiC DECONVOLUTION ===\n")

  library(MuSiC)

  sc_ref <- readRDS(sc_ref_file)

  ## Build bulk ExpressionSet
  bulk_counts <- counts(dds)
  bulk_phenoData <- new("AnnotatedDataFrame", data = sample_meta)
  bulk_eset <- ExpressionSet(assayData = as.matrix(bulk_counts),
                              phenoData = bulk_phenoData)

  ## Run MuSiC
  music_result <- tryCatch({
    music_prop(
      bulk.mtx   = exprs(bulk_eset),
      sc.sce     = sc_ref,
      clusters   = gsva_params$sc_celltype_col,
      samples    = gsva_params$sc_subject_col
    )
  }, error = function(e) {
    cat("[ERROR] MuSiC deconvolution failed: ", conditionMessage(e), "\n", sep = "")
    NULL
  })

  if (!is.null(music_result)) {
    proportions <- as.data.frame(music_result$Est.prop.weighted)
    proportions$Sample <- rownames(proportions)

    prop_file <- paste0("deconvolution_proportions_", run_tag, ".csv")
    write.csv(proportions, file = file.path(tables_dir, prop_file), row.names = FALSE)
    cat("[OK] Saved deconvolution results: ", prop_file, "\n", sep = "")
  }

} else {
  ## ---- MARKER-BASED PROPORTION ESTIMATION (FALLBACK) ----
  ## Method: Mean z-scored expression of signature genes per cell type
  ## Then softmax to get relative proportions (sum to 1 per sample)

  ## Use only cell-type signatures (not process signatures)
  deconv_sigs <- sig_filtered[intersect(names(sig_filtered), celltype_sigs)]

  if (length(deconv_sigs) >= 3) {

    ## Z-score the expression matrix
    vst_z <- t(scale(t(vst_mat)))

    ## Mean z-score per signature per sample
    mean_z_scores <- sapply(deconv_sigs, function(genes) {
      genes_present <- intersect(genes, rownames(vst_z))
      if (length(genes_present) < 3) return(rep(NA, ncol(vst_z)))
      colMeans(vst_z[genes_present, , drop = FALSE], na.rm = TRUE)
    })

    ## Softmax to get relative proportions
    ## (exponential normalization so rows sum to 1)
    mean_z_scores[is.na(mean_z_scores)] <- 0
    exp_scores <- exp(mean_z_scores)
    proportions <- exp_scores / rowSums(exp_scores)
    proportions <- as.data.frame(proportions)
    proportions$Sample <- rownames(proportions)

    cat("[OK] Estimated relative proportions for ", ncol(proportions) - 1, " cell types\n", sep = "")

    prop_file <- paste0("marker_proportions_", run_tag, ".csv")
    write.csv(proportions, file = file.path(tables_dir, prop_file), row.names = FALSE)
    cat("[OK] Saved: ", prop_file, "\n", sep = "")

    ## Merge with metadata
    prop_long <- proportions %>%
      pivot_longer(cols = -Sample, names_to = "CellType", values_to = "Proportion") %>%
      left_join(sample_meta %>%
                  select(Sample, any_of(c("Genotype", "Tissue"))),
                by = "Sample") %>%
      mutate(Group = ifelse(grepl("KD|KAT8", Genotype, ignore.case = TRUE), "KD", "CTL"))

    ## Statistical testing
    prop_stats <- list()
    for (ct in setdiff(colnames(proportions), "Sample")) {
      if (has_tissue) {
        for (tissue in unique(prop_long$Tissue)) {
          subset_data <- prop_long %>% filter(CellType == ct, Tissue == tissue)
          kd_vals <- subset_data$Proportion[subset_data$Group == "KD"]
          ctl_vals <- subset_data$Proportion[subset_data$Group == "CTL"]

          if (length(kd_vals) >= 2 && length(ctl_vals) >= 2) {
            wt <- tryCatch(wilcox.test(kd_vals, ctl_vals, exact = FALSE),
                           error = function(e) list(p.value = NA))
            prop_stats[[paste0(tissue, "_", ct)]] <- data.frame(
              Tissue = tissue, CellType = ct,
              Mean_KD = mean(kd_vals), Mean_CTL = mean(ctl_vals),
              Diff = mean(kd_vals) - mean(ctl_vals),
              P_value = wt$p.value,
              stringsAsFactors = FALSE
            )
          }
        }
      } else {
        kd_vals <- prop_long$Proportion[prop_long$CellType == ct & prop_long$Group == "KD"]
        ctl_vals <- prop_long$Proportion[prop_long$CellType == ct & prop_long$Group == "CTL"]

        wt <- tryCatch(wilcox.test(kd_vals, ctl_vals, exact = FALSE),
                       error = function(e) list(p.value = NA))
        prop_stats[[ct]] <- data.frame(
          Tissue = "All", CellType = ct,
          Mean_KD = mean(kd_vals), Mean_CTL = mean(ctl_vals),
          Diff = mean(kd_vals) - mean(ctl_vals),
          P_value = wt$p.value,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(prop_stats) > 0) {
      prop_stats_df <- bind_rows(prop_stats)
      prop_stats_df$FDR <- p.adjust(prop_stats_df$P_value, method = "BH")
      prop_stats_df <- prop_stats_df %>% arrange(P_value)

      prop_stat_file <- paste0("proportion_statistics_", run_tag, ".csv")
      write.csv(prop_stats_df, file = file.path(tables_dir, prop_stat_file), row.names = FALSE)
      cat("[OK] Saved proportion statistics: ", prop_stat_file, "\n", sep = "")

      cat("\nSignificant proportion changes (p < 0.05):\n")
      sig_props <- prop_stats_df %>% filter(P_value < 0.05)
      if (nrow(sig_props) > 0) {
        for (i in seq_len(nrow(sig_props))) {
          row <- sig_props[i, ]
          dir_label <- ifelse(row$Diff > 0, "HIGHER in KD", "LOWER in KD")
          cat("  ", row$Tissue, " - ", row$CellType, ": ", dir_label,
              " (p = ", format(row$P_value, digits = 3), ")\n", sep = "")
        }
      } else {
        cat("  None at p < 0.05\n")
      }
    }

    ## ---- Stacked bar plot ----
    cat("\n=== GENERATING PROPORTION PLOTS ===\n")

    ## Order cell types for consistent plotting
    ct_order <- intersect(celltype_sigs, colnames(proportions))
    prop_long$CellType <- factor(prop_long$CellType, levels = rev(ct_order))

    ## Color palette for cell types
    n_ct <- length(ct_order)
    ct_colors <- setNames(
      colorRampPalette(brewer.pal(min(n_ct, 12), "Set3"))(n_ct),
      ct_order
    )

    p_stacked <- ggplot(prop_long, aes(x = Sample, y = Proportion, fill = CellType)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
      scale_fill_manual(values = ct_colors, name = "Cell Type") +
      labs(
        title = "Estimated Cell-Type Proportions (Marker-Based)",
        subtitle = "Relative proportions from mean z-scored signature gene expression",
        x = NULL, y = "Estimated Proportion"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right"
      )

    if (has_tissue) {
      p_stacked <- p_stacked + facet_wrap(~ Tissue, scales = "free_x")
    }

    stacked_file <- paste0("CellType_proportions_stacked_", run_tag, ".png")
    ggsave(file.path(gsva_dir, stacked_file), plot = p_stacked,
           width = max(8, ncol(vst_mat) * 0.6 + 3), height = 6, dpi = 300, bg = "white")
    cat("[OK] Saved: ", stacked_file, "\n", sep = "")

    ## ---- Proportion boxplots per cell type ----
    p_prop_box <- ggplot(prop_long, aes(x = Group, y = Proportion, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
      facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c("CTL" = "#0072B2", "KD" = "#D55E00")) +
      labs(
        title = "Cell-Type Proportion Estimates: KD vs CTL",
        x = NULL, y = "Estimated Proportion"
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom"
      )

    prop_box_file <- paste0("CellType_proportion_boxplots_", run_tag, ".png")
    ggsave(file.path(gsva_dir, prop_box_file), plot = p_prop_box,
           width = 12, height = max(6, ceiling(n_ct / 4) * 3), dpi = 300, bg = "white")
    cat("[OK] Saved: ", prop_box_file, "\n", sep = "")
  }
}

## ============================================================
## 10) INTERPRETATION NOTES
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== INTERPRETATION NOTES ===\n")
cat("============================================================\n\n")

cat("GSVA SCORES:\n")
cat("  - Positive score = signature genes are enriched (highly expressed) in that sample\n")
cat("  - Negative score = signature genes are depleted (lowly expressed) in that sample\n")
cat("  - Scores are relative within each sample (not absolute abundances)\n\n")

cat("MARKER-BASED PROPORTIONS:\n")
cat("  - These are RELATIVE estimates, not absolute cell counts\n")
cat("  - They indicate which signatures dominate in each sample\n")
cat("  - For absolute proportions, use MuSiC with an scRNA-seq reference\n\n")

cat("TO UPGRADE TO FULL DECONVOLUTION:\n")
cat("  1. Download Emont et al. 2022 scRNA-seq reference from GEO: GSE176171\n")
cat("  2. Process into a SingleCellExperiment object with cell type annotations\n")
cat("  3. Save as RDS: saveRDS(sce, 'emont_2022_mouse_wat_reference.rds')\n")
cat("  4. Set in parameters.R:\n")
cat("       gsva_params$sc_reference_path <- 'emont_2022_mouse_wat_reference.rds'\n")
cat("       gsva_params$sc_celltype_col <- 'cell_type'\n")
cat("       gsva_params$sc_subject_col <- 'subject'\n")
cat("  5. Install MuSiC: devtools::install_github('xuranw/MuSiC')\n")
cat("  6. Re-run this script\n\n")

## ============================================================
## 11) SESSION INFO
## ============================================================

session_file <- file.path(outdir, "logs", paste0("sessionInfo_gsva_", run_tag, ".txt"))
if (dir.exists(dirname(session_file))) {
  sink(session_file)
  cat("=== R SESSION INFORMATION - Part 6 (GSVA & Deconvolution) ===\n\n")
  print(sessionInfo())
  sink()
  cat("[INFO] Session info saved\n")
}

cat("\n=== Part 6 (GSVA Cell-Type Scoring & Deconvolution) COMPLETE ===\n")
cat("Output directory: ", gsva_dir, "\n", sep = "")
