#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 1b: QC DIAGNOSTICS
## =========================================================
## Runs after Part 1, before the pathway analyses. Everything here is a check
## you should look at BEFORE trusting a result -- these are the diagnostics
## that catch a broken model or a bad sample, and they were the main gap in
## this pipeline's QC.
##
##   A. p-value histograms          is the DE model calibrated?
##   B. independent filtering       how many genes got padj = NA, and why
##   C. unstable estimates          genes whose log2FC is driven by ~1 sample
##   D. contamination screen        foreign tissue in a fat sample, PER SAMPLE
##   E. PCA scree + PC-covariate    what actually drives the main components
##   F. sample-sample correlation   outlier detection independent of PCA
##   G. PCA convention comparison   top-variable genes vs all genes, side by side
##
## Reads Part 1's QC_data_*.rds and DE_tissue_*.csv from the current run folder.
## Writes plots + a QC_FLAGS_*.csv summarising everything that needs attention.
## =========================================================

required_pkgs <- c("dplyr", "ggplot2", "tidyr", "pheatmap", "RColorBrewer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
  library(pheatmap); library(RColorBrewer)
})

params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) source(params_file)

## ---- locate the run folder (same logic as Parts 2/3/4) ----
savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) stop("savepoints/ not found - run Part 1 first.", call. = FALSE)
outdir <- NULL
ptr <- file.path(savepoint_dir, "LATEST_RUN.txt")
if (file.exists(ptr)) {
  cand <- trimws(readLines(ptr, warn = FALSE)); cand <- cand[nzchar(cand)]
  if (length(cand) && dir.exists(cand[1])) outdir <- cand[1]
}
if (is.null(outdir)) {
  rd <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  rd <- rd[grepl("^RUN_", basename(rd))]
  if (!length(rd)) stop("No RUN_ folders found.", call. = FALSE)
  outdir <- rd[order(basename(rd), decreasing = TRUE)][1]
}
.clean_run_tag <- function(dir) {
  b <- basename(dir); if (b %in% c("tissue", "cells")) b <- basename(dirname(dir))
  b <- sub("^RUN_", "", b)
  m <- regmatches(b, regexpr("[0-9]{8}_[0-9]{6}$", b))
  if (length(m) == 1 && nzchar(m)) m else b
}
run_tag   <- .clean_run_tag(outdir)
tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
cat("[INFO] QC diagnostics for: ", outdir, "\n", sep = "")

FLAGS <- list()
flag <- function(check, value, status = "OK", action = "") {
  FLAGS[[length(FLAGS) + 1]] <<- data.frame(check = check, value = as.character(value),
                                            status = status, action = action,
                                            stringsAsFactors = FALSE)
  cat(sprintf("  [%-4s] %-44s %s\n", status, check, value))
}

## =========================================================
## A. p-VALUE HISTOGRAMS  -- the single best model diagnostic
## =========================================================
## A well-behaved contrast gives a FLAT histogram with a spike near 0.
##   flat + spike at 0      = healthy
##   hill in the middle     = test is mis-calibrated / model misspecified
##   spike at 1             = over-conservative filtering
## Anything other than the first pattern means the DE results are suspect,
## regardless of how good the biology looks.
cat("\n=== A. p-value histograms ===\n")
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
if (!length(de_files)) stop("No DE_tissue_*.csv found in ", tables_dir, call. = FALSE)

pv_list <- list()
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  pv <- d$pvalue[is.finite(d$pvalue)]
  pv_list[[cn]] <- data.frame(contrast = cn, pvalue = pv)

  ## Fraction in the flat region: for a calibrated test the density of p in
  ## [0.5,1] should be ~= the null density. A big hump in the middle shows up
  ## as the mid bins exceeding the high bins.
  h   <- hist(pv, breaks = seq(0, 1, 0.05), plot = FALSE)$counts
  mid <- mean(h[9:12]); high <- mean(h[17:20])
  ratio <- mid / max(high, 1)
  flag(paste0(cn, ": p-hist mid/high ratio"), sprintf("%.2f", ratio),
       ifelse(ratio > 1.5, "WARN", "OK"),
       ifelse(ratio > 1.5, "hump in the middle -> check model/covariates", ""))
  flag(paste0(cn, ": fraction p < 0.01"),
       sprintf("%.1f%%", 100 * mean(pv < 0.01)), "INFO")
}
pv_all <- dplyr::bind_rows(pv_list)
p <- ggplot(pv_all, aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.02), fill = "grey70", colour = "black", linewidth = .2) +
  facet_wrap(~ contrast, scales = "free_y") +
  labs(title = "p-value histograms (raw p, all tested genes)",
       subtitle = "healthy = flat with a spike at 0; a hump in the middle means the model is mis-specified",
       x = "raw p-value", y = "genes") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(plots_dir, paste0("QC_pvalue_histograms_", run_tag, ".png")),
       p, width = 9, height = 4.2, dpi = 300, bg = "white")

## =========================================================
## B. INDEPENDENT FILTERING
## =========================================================
cat("\n=== B. independent filtering ===\n")
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  n_na <- sum(is.na(d$padj))
  flag(paste0(cn, ": genes with padj = NA"),
       sprintf("%d / %d (%.1f%%)", n_na, nrow(d), 100 * n_na / nrow(d)),
       ifelse(n_na / nrow(d) > 0.5, "WARN", "OK"),
       ifelse(n_na / nrow(d) > 0.5, "most genes filtered out - check prefilter", ""))
}

## =========================================================
## C. UNSTABLE ESTIMATES  -- log2FC driven by very few samples
## =========================================================
## A huge |log2FC| paired with a huge lfcSE is not a big biological effect; it
## means the gene is ~zero in one group and nonzero in a couple of samples of
## the other. These dominate volcano plots and look impressive but are usually
## contamination or a single outlier library.
cat("\n=== C. unstable log2FC estimates ===\n")
unstable_all <- list()
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  d$gene <- if ("gene_name" %in% names(d)) d$gene_name else d[[1]]
  u <- d %>% dplyr::filter(!is.na(padj), padj < 0.05, lfcSE > 2) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, padj) %>%
    dplyr::mutate(contrast = cn)
  unstable_all[[cn]] <- u
  flag(paste0(cn, ": significant genes with lfcSE > 2"), nrow(u),
       ifelse(nrow(u) > 0, "WARN", "OK"),
       ifelse(nrow(u) > 0, "inspect per-sample counts before reporting these", ""))
  if (nrow(u)) cat("        ", paste(utils::head(u$gene, 10), collapse = ", "), "\n")
}
unstable <- dplyr::bind_rows(unstable_all)
if (nrow(unstable))
  write.csv(unstable, file.path(tables_dir, paste0("QC_unstable_estimates_", run_tag, ".csv")),
            row.names = FALSE)

## =========================================================
## D. CONTAMINATION SCREEN  -- per sample
## =========================================================
## Fat pads are dissected next to other tissues. A single sample carrying
## epididymis / testis / muscle can create "significant" DEGs on its own at
## n = 5 per group. This scores every sample against marker panels and, more
## importantly, tests whether contamination is CONFOUNDED WITH GENOTYPE --
## which is the case that actually biases results.
cat("\n=== D. contamination screen (per sample) ===\n")
qc_files <- list.files(tables_dir, pattern = "^QC_data_.*\\.rds$", full.names = TRUE)
if (!length(qc_files)) {
  flag("QC_data_*.rds", "NOT FOUND - contamination screen skipped", "WARN",
       "re-run Part 1 to regenerate")
} else {
  qc <- readRDS(qc_files[order(basename(qc_files), decreasing = TRUE)][1])
  vst <- qc$vst_mat_before; if (is.null(vst)) vst <- qc$vst_mat_qc
  si  <- qc$sample_info
  si  <- si[match(colnames(vst), si$Sample), , drop = FALSE]

  PANELS <- list(
    Sperm_testis = c("Prm1","Prm2","Tnp1","Tnp2","Smcp","Oaz3","Odf1","Akap4","Spata19","Izumo1"),
    Blood        = c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Alas2"),
    Muscle       = c("Acta1","Myh1","Myh2","Ckm","Tnnt3","Mylpf"),
    Liver        = c("Alb","Apoa1","Ttr","Serpina1a"),
    LymphNode    = c("Cd19","Ms4a1","Cd3e","Ccl21a","Ccl19"))

  score_rows <- list()
  for (nm in names(PANELS)) {
    g <- intersect(PANELS[[nm]], rownames(vst))
    if (length(g) < 2) { flag(paste0(nm, " panel"), "too few genes measured", "INFO"); next }
    sc <- colMeans(vst[g, , drop = FALSE])
    z  <- (sc - median(sc)) / (mad(sc) + 1e-9)          ## robust z, outlier-resistant
    score_rows[[nm]] <- data.frame(panel = nm, Sample = names(sc), score = as.numeric(sc),
                                   robust_z = as.numeric(z),
                                   Genotype = si$Genotype, Depot = si$Depot, Sex = si$Sex,
                                   stringsAsFactors = FALSE)
    hot <- names(sc)[z > 5]
    flag(paste0(nm, ": samples with robust z > 5"),
         ifelse(length(hot), paste(hot, collapse = ", "), "none"),
         ifelse(length(hot), "WARN", "OK"),
         ifelse(length(hot), "inspect these libraries", ""))

    ## The dangerous case: contamination tracking with genotype.
    if (length(unique(si$Genotype)) == 2) {
      pv <- tryCatch(suppressWarnings(wilcox.test(sc ~ si$Genotype)$p.value), error = function(e) NA)
      flag(paste0(nm, ": contamination vs GENOTYPE"), sprintf("p = %.3g", pv),
           ifelse(!is.na(pv) && pv < 0.05, "WARN", "OK"),
           ifelse(!is.na(pv) && pv < 0.05,
                  "CONFOUNDED with genotype -> will create false DEGs", ""))
    }
  }
  cont <- dplyr::bind_rows(score_rows)
  if (nrow(cont)) {
    write.csv(cont, file.path(tables_dir, paste0("QC_contamination_scores_", run_tag, ".csv")),
              row.names = FALSE)
    p <- ggplot(cont, aes(x = Sample, y = score, fill = Genotype)) +
      geom_col(colour = "black", linewidth = .2) +
      facet_wrap(~ panel, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c(CTL = "#0072B2", KAT8KD = "#D55E00")) +
      labs(title = "Contamination screen: foreign-tissue marker score per sample",
           subtitle = "one sample standing far above the rest = likely dissection carry-over",
           x = NULL, y = "mean VST of panel") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
            panel.grid.minor = element_blank())
    ggsave(file.path(plots_dir, paste0("QC_contamination_", run_tag, ".png")),
           p, width = 11, height = 12, dpi = 300, bg = "white")
  }

  ## =========================================================
  ## E. PCA SCREE + WHAT DRIVES EACH COMPONENT
  ## =========================================================
  ## A PCA plot shows structure but not its CAUSE. This tests every PC against
  ## the known covariates, so "PC1 is depot" or "PC2 is library size" becomes a
  ## number rather than an impression.
  cat("\n=== E. PCA scree + PC-covariate association ===\n")
  rv <- apply(vst, 1, stats::var)
  idx <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(vst[idx, , drop = FALSE]), scale. = FALSE)
  ve  <- summary(pca)$importance[2, ] * 100
  scree <- data.frame(PC = factor(paste0("PC", 1:min(10, length(ve))),
                                  levels = paste0("PC", 1:min(10, length(ve)))),
                      variance = as.numeric(ve[1:min(10, length(ve))]))
  p <- ggplot(scree, aes(PC, variance)) +
    geom_col(fill = "grey70", colour = "black", linewidth = .2) +
    geom_text(aes(label = sprintf("%.1f%%", variance)), vjust = -0.4, size = 3) +
    labs(title = "PCA scree", x = NULL, y = "% variance explained") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("QC_PCA_scree_", run_tag, ".png")),
         p, width = 7, height = 4, dpi = 300, bg = "white")

  libsize <- if (!is.null(qc$lib_sizes)) qc$lib_sizes[colnames(vst)] else NULL
  covs <- list(Genotype = si$Genotype, Depot = si$Depot, Sex = si$Sex)
  if (!is.null(libsize)) covs$LibrarySize <- as.numeric(libsize)
  for (pc in 1:min(4, ncol(pca$x))) {
    for (cvn in names(covs)) {
      cv <- covs[[cvn]]
      pv <- tryCatch({
        if (is.numeric(cv)) cor.test(pca$x[, pc], cv)$p.value
        else summary(aov(pca$x[, pc] ~ factor(cv)))[[1]][["Pr(>F)"]][1]
      }, error = function(e) NA_real_)
      if (!is.na(pv) && pv < 0.05)
        flag(sprintf("PC%d ~ %s", pc, cvn), sprintf("p = %.2g", pv), "INFO",
             "this covariate explains part of this component")
    }
  }
  if (is.null(libsize))
    flag("PC ~ LibrarySize", "lib_sizes not saved by Part 1", "INFO",
         "re-run Part 1 to enable this test")

  ## =========================================================
  ## F. SAMPLE-SAMPLE CORRELATION  -- outliers PCA can miss
  ## =========================================================
  cat("\n=== F. sample-sample correlation ===\n")
  cm <- cor(vst[idx, , drop = FALSE], method = "spearman")
  ann <- data.frame(Genotype = si$Genotype, Depot = si$Depot, row.names = si$Sample)
  png(file.path(plots_dir, paste0("QC_sample_correlation_", run_tag, ".png")),
      width = 1600, height = 1500, res = 180)
  pheatmap(cm, annotation_col = ann, annotation_row = ann,
           color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
           main = "Sample-sample Spearman correlation (top 500 variable genes)",
           fontsize = 7)
  dev.off()
  mean_cor <- (rowSums(cm) - 1) / (ncol(cm) - 1)
  z <- (mean_cor - median(mean_cor)) / (mad(mean_cor) + 1e-9)
  low <- names(mean_cor)[z < -5]
  flag("samples poorly correlated with the rest (z < -5)",
       ifelse(length(low), paste(low, collapse = ", "), "none"),
       ifelse(length(low), "WARN", "OK"),
       ifelse(length(low), "candidate outlier libraries", ""))

  ## =========================================================
  ## G. PCA CONVENTION COMPARISON
  ## =========================================================
  ## The two conventions genuinely look different and neither is wrong, so
  ## both are drawn side by side rather than one being chosen silently.
  cat("\n=== G. PCA convention comparison ===\n")
  mk <- function(mat, scale_it, label) {
    pr <- prcomp(t(mat), scale. = scale_it)
    v  <- summary(pr)$importance[2, 1:2] * 100
    data.frame(Sample = colnames(mat), PC1 = pr$x[, 1], PC2 = pr$x[, 2],
               Genotype = si$Genotype, Depot = si$Depot,
               panel = sprintf("%s  (PC1 %.1f%%, PC2 %.1f%%)", label, v[1], v[2]))
  }
  cmp <- rbind(mk(vst[idx, , drop = FALSE], FALSE, "top 500 var, unscaled (DESeq2 default)"),
               mk(vst,                      TRUE,  "all genes, unit-scaled (original)"))
  p <- ggplot(cmp, aes(PC1, PC2, colour = Depot, shape = Genotype)) +
    geom_point(size = 2.6) + facet_wrap(~ panel, scales = "free") +
    labs(title = "PCA: the two conventions compared",
         subtitle = "same data, different gene selection/scaling. PC sign is arbitrary; a mirrored panel is not an error.") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("QC_PCA_convention_comparison_", run_tag, ".png")),
         p, width = 11, height = 4.6, dpi = 300, bg = "white")
}

## ---- summary ----
fl <- dplyr::bind_rows(FLAGS)
write.csv(fl, file.path(tables_dir, paste0("QC_FLAGS_", run_tag, ".csv")), row.names = FALSE)
nw <- sum(fl$status == "WARN")
cat("\n==========================================================\n")
cat("QC DIAGNOSTICS: ", nrow(fl), " checks | ", nw, " need attention\n", sep = "")
if (nw > 0) print(fl[fl$status == "WARN", c("check", "value", "action")], row.names = FALSE)
cat("==========================================================\n")
cat("[DONE] Part 1b. Plots in ", plots_dir, "\n", sep = "")
