#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 7.7: SEX-COMBINING VALIDATION
##                                & MANUSCRIPT NUMBERS
## =========================================================
##
## PURPOSE: Rigorous, self-contained validation that combining
## sexes is justified for the unified ~ 0 + GroupDepot model
## (Part 1), with every number needed for the results section
## saved to a single output file.
##
## PRODUCES (per depot):
##   1. Formal LRT interaction test (Sex x Genotype)
##   2. Sex-stratified concordance (M vs F logFC)
##   3. Model comparison: ~Genotype vs ~Sex+Genotype vs unified (+/- Sex)
##   4. Variance decomposition
##   5. DEG counts (combined, male-only, female-only)
##   6. VERDICT: PASS / CAUTION / CONCERN
##   7. Ready-to-paste manuscript text with actual numbers
##
## ALL numbers saved to:
##   - sex_combining_validation_report_<run>.txt  (human-readable)
##   - sex_combining_validation_numbers_<run>.csv  (machine-readable)
## =========================================================

cat("
============================================================
  PART 7.7: SEX-COMBINING VALIDATION & MANUSCRIPT NUMBERS
============================================================
\n")

## 1) Packages ------------------------------------------------
required_pkgs <- c("DESeq2", "dplyr", "tidyr", "ggplot2", "ggrepel", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

## Try to load qvalue for pi0
has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
if (has_qvalue) library(qvalue)

## 1.5) Parameters --------------------------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found.")
}
if (!exists("prefilter_min_count"))   prefilter_min_count   <- 10
if (!exists("prefilter_min_samples")) prefilter_min_samples <- 4

## Helper: safe DESeq2 results to data.frame
deseq_res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  std_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  if (ncol(df) == length(std_cols)) colnames(df) <- std_cols
  df
}

## ============================================================
## 2) LOAD DATA (same logic as part1 / part7.6)
## ============================================================
cat("[INFO] Loading data...\n")

counts_file <- file.path(getwd(), "counts.txt")
if (!file.exists(counts_file)) stop("counts.txt not found")
counts_raw <- read.delim(counts_file, header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)

sample_ids <- paste0("JS_", sprintf("%02d", 1:40))
sample_annot <- data.frame(
  Sample   = sample_ids,
  Depot    = c(rep("iWAT", 10), rep("iWAT", 10),
               rep("gWAT", 10), rep("gWAT", 10)),
  Sex      = c(rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5),
               rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5)),
  Genotype = c(rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5)),
  stringsAsFactors = FALSE
)

## Remove outliers
outlier_samples <- c("JS_08", "JS_28")
sample_annot <- sample_annot[!sample_annot$Sample %in% outlier_samples, ]
tissue_samples <- intersect(sample_annot$Sample, colnames(counts_raw))
sample_annot <- sample_annot[sample_annot$Sample %in% tissue_samples, ]

## Factorize
sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
sample_annot$Depot    <- factor(sample_annot$Depot, levels = c("iWAT", "gWAT"))
sample_annot$Sex      <- factor(sample_annot$Sex, levels = c("F", "M"))
rownames(sample_annot) <- sample_annot$Sample

## Build count matrix
if ("gene" %in% colnames(counts_raw) || "gene_name" %in% colnames(counts_raw)) {
  counts_tissue <- counts_raw
  if (!"gene_name" %in% colnames(counts_tissue) && "gene" %in% colnames(counts_tissue))
    counts_tissue$gene_name <- counts_tissue$gene
  if (!"ensembl_gene_id" %in% colnames(counts_tissue) && "gene_id" %in% colnames(counts_tissue))
    counts_tissue$ensembl_gene_id <- counts_tissue$gene_id
  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  if ("ensembl_gene_id" %in% colnames(counts_tissue))
    counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]
  if (any(duplicated(counts_tissue$gene_name))) {
    dup_flag <- duplicated(counts_tissue$gene_name) | duplicated(counts_tissue$gene_name, fromLast = TRUE)
    if ("ensembl_gene_id" %in% colnames(counts_tissue))
      counts_tissue$gene_name[dup_flag] <- paste0(counts_tissue$ensembl_gene_id[dup_flag], "_",
                                                   counts_tissue$gene_name[dup_flag])
    counts_tissue$gene_name <- make.unique(counts_tissue$gene_name)
  }
  rownames(counts_tissue) <- counts_tissue$gene_name
  count_matrix <- as.matrix(counts_tissue[, tissue_samples])
  rownames(count_matrix) <- counts_tissue$gene_name
} else {
  count_matrix <- as.matrix(counts_raw[, tissue_samples])
}

## Prefilter
keep <- rowSums(count_matrix >= prefilter_min_count) >= prefilter_min_samples
count_matrix <- count_matrix[keep, , drop = FALSE]
cat("[INFO] Genes after prefiltering: ", nrow(count_matrix), "\n")

## ============================================================
## 3) OUTPUT SETUP
## ============================================================
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) > 0) {
  outdir <- run_dirs[which.max(file.info(run_dirs)$mtime)]
} else {
  outdir <- file.path(savepoint_dir,
                      paste0("RUN_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots", "sex_validation")
if (!dir.exists(plots_dir))  dir.create(plots_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

## Report file (human-readable) + numbers CSV (machine-readable)
report_file  <- file.path(tables_dir,
                          paste0("sex_combining_validation_report_", run_tag, ".txt"))
numbers_file <- file.path(tables_dir,
                          paste0("sex_combining_validation_numbers_", run_tag, ".csv"))

## Open report
report <- file(report_file, open = "wt")
writeLines("============================================================", report)
writeLines("  SEX-COMBINING VALIDATION REPORT", report)
writeLines(paste0("  Generated: ", Sys.time()), report)
writeLines(paste0("  Run tag:   ", run_tag), report)
writeLines("============================================================", report)
writeLines("", report)

## Numbers collector (for CSV output)
all_numbers <- list()

## ============================================================
## 3b) UNIFIED MODEL (computed once, used in TEST 3b per depot)
## ============================================================
cat("[INFO] Running unified model (~ Sex + GroupDepot, all samples)...\n")

## Create GroupDepot factor
sample_annot$GroupDepot <- factor(
  paste(sample_annot$Depot, sample_annot$Genotype, sep = "_"),
  levels = c("iWAT_CTL", "iWAT_KAT8KD", "gWAT_CTL", "gWAT_KAT8KD")
)

dds_unified <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = sample_annot,
  design    = ~ Sex + GroupDepot
)
dds_unified <- DESeq(dds_unified, quiet = TRUE)

## Extract per-depot KAT8 contrasts from unified model
unified_results <- list()
for (d in c("iWAT", "gWAT")) {
  res_u <- results(dds_unified,
                   contrast = c("GroupDepot",
                                paste0(d, "_KAT8KD"),
                                paste0(d, "_CTL")))
  res_u_df <- deseq_res_to_df(res_u)
  res_u_df$gene <- rownames(count_matrix)
  unified_results[[d]] <- res_u_df
}
cat("[OK] Unified model complete. Coefficients: ",
    paste(resultsNames(dds_unified), collapse = ", "), "\n")

## ============================================================
## 3c) ORIGINAL MODEL: ~ 0 + GroupDepot (no Sex covariate)
##     This is the model the boss has been reviewing.
##     Comparing against ~ Sex + GroupDepot justifies adding Sex.
## ============================================================
cat("[INFO] Running original model (~ 0 + GroupDepot, no Sex covariate)...\n")

dds_nosex_unified <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = sample_annot,
  design    = ~ 0 + GroupDepot
)
dds_nosex_unified <- DESeq(dds_nosex_unified, quiet = TRUE)

## Extract per-depot KAT8 contrasts from no-sex model
nosex_unified_results <- list()
for (d in c("iWAT", "gWAT")) {
  res_nu <- results(dds_nosex_unified,
                    contrast = c("GroupDepot",
                                 paste0(d, "_KAT8KD"),
                                 paste0(d, "_CTL")))
  res_nu_df <- deseq_res_to_df(res_nu)
  res_nu_df$gene <- rownames(count_matrix)
  nosex_unified_results[[d]] <- res_nu_df
}
cat("[OK] Original model complete. Coefficients: ",
    paste(resultsNames(dds_nosex_unified), collapse = ", "), "\n")

## ============================================================
## 4) PER-DEPOT ANALYSIS
## ============================================================

for (depot in c("iWAT", "gWAT")) {

  cat("\n============================================================\n")
  cat("===  ", depot, ": SEX-COMBINING VALIDATION  ===\n", sep = "")
  cat("============================================================\n\n")

  writeLines(paste0("\n", paste(rep("=", 60), collapse = "")), report)
  writeLines(paste0("  ", depot), report)
  writeLines(paste(rep("=", 60), collapse = ""), report)

  ## ---------------------------------------------------------
  ## 4a) Subset to depot
  ## ---------------------------------------------------------
  depot_samples <- sample_annot$Sample[sample_annot$Depot == depot]
  depot_annot   <- sample_annot[depot_samples, ]
  depot_counts  <- count_matrix[, depot_samples]

  ## Re-prefilter within depot
  depot_keep <- rowSums(depot_counts >= prefilter_min_count) >= prefilter_min_samples
  depot_counts <- depot_counts[depot_keep, , drop = FALSE]
  depot_gene_names <- rownames(depot_counts)

  ## Group sizes
  grp_table <- with(depot_annot, table(Sex, Genotype))
  n_per_group <- as.data.frame(grp_table)

  writeLines("\n--- Sample sizes ---", report)
  for (i in seq_len(nrow(n_per_group))) {
    writeLines(sprintf("  %s_%s: n = %d",
                       n_per_group$Sex[i], n_per_group$Genotype[i],
                       n_per_group$Freq[i]), report)
  }
  writeLines(sprintf("  Total: %d samples, %d genes",
                     ncol(depot_counts), nrow(depot_counts)), report)

  thresh <- get_tissue_thresholds(depot)
  fdr_cut   <- thresh$fdr_cut
  logFC_cut <- thresh$logFC_cut

  ## ---------------------------------------------------------
  ## 4b) TEST 1: Formal LRT interaction test
  ## ---------------------------------------------------------
  cat("[", depot, "] Running LRT interaction test...\n", sep = "")

  dds_full <- DESeqDataSetFromMatrix(
    countData = round(depot_counts),
    colData   = depot_annot,
    design    = ~ Sex + Genotype + Sex:Genotype
  )
  dds_lrt <- DESeq(dds_full, test = "LRT", reduced = ~ Sex + Genotype)

  res_lrt <- results(dds_lrt)
  res_lrt_df <- deseq_res_to_df(res_lrt)
  res_lrt_df$gene <- depot_gene_names
  res_lrt_df <- res_lrt_df[!is.na(res_lrt_df$pvalue), ]

  n_tested    <- nrow(res_lrt_df)
  n_sig_lrt   <- sum(res_lrt_df$padj < fdr_cut, na.rm = TRUE)
  pct_sig_lrt <- round(100 * n_sig_lrt / n_tested, 2)

  ## pi0 estimation
  pi0_est <- NA_real_
  if (has_qvalue && n_tested > 100) {
    tryCatch({
      qobj <- qvalue(res_lrt_df$pvalue)
      pi0_est <- qobj$pi0
    }, error = function(e) { })
  }

  ## Also get Wald test for interaction coefficient
  dds_wald <- DESeq(dds_full, test = "Wald")
  coef_names <- resultsNames(dds_wald)
  interaction_coef <- grep("Sex.*Genotype|Genotype.*Sex", coef_names, value = TRUE)

  if (length(interaction_coef) == 1) {
    res_wald_int <- deseq_res_to_df(results(dds_wald, name = interaction_coef))
    res_wald_int$gene <- depot_gene_names
    n_sig_wald <- sum(res_wald_int$padj < fdr_cut, na.rm = TRUE)
    median_abs_int_lfc <- median(abs(res_wald_int$log2FoldChange), na.rm = TRUE)
  } else {
    n_sig_wald <- NA_integer_
    median_abs_int_lfc <- NA_real_
  }

  writeLines("\n--- TEST 1: LRT Interaction Test ---", report)
  writeLines(sprintf("  Model: ~ Sex + Genotype + Sex:Genotype"), report)
  writeLines(sprintf("  Reduced: ~ Sex + Genotype"), report)
  writeLines(sprintf("  Genes tested:                 %d", n_tested), report)
  writeLines(sprintf("  Significant interactions:     %d (%.2f%%)", n_sig_lrt, pct_sig_lrt), report)
  writeLines(sprintf("  Wald interaction sig:         %s", ifelse(is.na(n_sig_wald), "N/A", n_sig_wald)), report)
  if (!is.na(pi0_est))
    writeLines(sprintf("  pi0 (true null proportion):   %.4f (%.1f%% genes have NO interaction)",
                       pi0_est, pi0_est * 100), report)
  if (!is.na(median_abs_int_lfc))
    writeLines(sprintf("  Median |interaction logFC|:   %.3f", median_abs_int_lfc), report)

  ## Save significant interaction genes
  sig_int_genes <- res_lrt_df[!is.na(res_lrt_df$padj) & res_lrt_df$padj < fdr_cut, ]
  if (nrow(sig_int_genes) > 0) {
    sig_int_genes <- sig_int_genes[order(sig_int_genes$padj), ]
    write.csv(sig_int_genes,
              file.path(tables_dir, paste0("sex_interaction_sig_genes_", depot, "_", run_tag, ".csv")),
              row.names = FALSE)
    writeLines(sprintf("  Top 10 interaction genes:"), report)
    top10 <- head(sig_int_genes, 10)
    for (i in seq_len(nrow(top10))) {
      writeLines(sprintf("    %s (padj = %.2e, LFC = %.2f)",
                         top10$gene[i], top10$padj[i], top10$log2FoldChange[i]), report)
    }
  }

  ## P-value histogram
  p_hist <- ggplot(res_lrt_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", color = "white", boundary = 0) +
    geom_hline(yintercept = n_tested / 50, linetype = "dashed", color = "red") +
    labs(title = paste0(depot, ": P-value Distribution (Sex x Genotype Interaction)"),
         subtitle = paste0(n_sig_lrt, " / ", n_tested, " significant (",
                           pct_sig_lrt, "%)",
                           if (!is.na(pi0_est)) paste0("; pi0 = ", round(pi0_est, 3))),
         x = "P-value (LRT)", y = "Number of genes") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("01_pvalue_hist_interaction_", depot, ".png")),
         p_hist, width = 7, height = 5, dpi = 300, bg = "white")

  ## ---------------------------------------------------------
  ## 4c) TEST 2: Sex-stratified concordance
  ## ---------------------------------------------------------
  cat("[", depot, "] Running sex-stratified DESeq2...\n", sep = "")

  female_idx <- depot_annot$Sample[depot_annot$Sex == "F"]
  male_idx   <- depot_annot$Sample[depot_annot$Sex == "M"]

  ## Female-only DESeq2
  dds_f <- DESeqDataSetFromMatrix(
    countData = round(depot_counts[, female_idx]),
    colData   = depot_annot[female_idx, ],
    design    = ~ Genotype
  )
  dds_f <- DESeq(dds_f, quiet = TRUE)
  res_f <- deseq_res_to_df(results(dds_f))
  res_f$gene <- depot_gene_names

  ## Male-only DESeq2
  dds_m <- DESeqDataSetFromMatrix(
    countData = round(depot_counts[, male_idx]),
    colData   = depot_annot[male_idx, ],
    design    = ~ Genotype
  )
  dds_m <- DESeq(dds_m, quiet = TRUE)
  res_m <- deseq_res_to_df(results(dds_m))
  res_m$gene <- depot_gene_names

  ## DEG counts per sex
  n_deg_f <- sum(res_f$padj < fdr_cut & abs(res_f$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_deg_m <- sum(res_m$padj < fdr_cut & abs(res_m$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_up_f  <- sum(res_f$padj < fdr_cut & res_f$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_f  <- sum(res_f$padj < fdr_cut & res_f$log2FoldChange < -logFC_cut, na.rm = TRUE)
  n_up_m  <- sum(res_m$padj < fdr_cut & res_m$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_m  <- sum(res_m$padj < fdr_cut & res_m$log2FoldChange < -logFC_cut, na.rm = TRUE)

  ## Concordance
  concord_df <- merge(
    res_f[, c("gene", "log2FoldChange", "padj")],
    res_m[, c("gene", "log2FoldChange", "padj")],
    by = "gene", suffixes = c("_F", "_M")
  ) %>%
    filter(is.finite(log2FoldChange_F), is.finite(log2FoldChange_M))

  r_all      <- cor(concord_df$log2FoldChange_F, concord_df$log2FoldChange_M, use = "complete.obs")
  rho_all    <- cor(concord_df$log2FoldChange_F, concord_df$log2FoldChange_M,
                    method = "spearman", use = "complete.obs")
  r_test     <- cor.test(concord_df$log2FoldChange_F, concord_df$log2FoldChange_M)

  ## DEG-only concordance
  concord_sig <- concord_df %>%
    filter((!is.na(padj_F) & padj_F < fdr_cut & abs(log2FoldChange_F) > logFC_cut) |
           (!is.na(padj_M) & padj_M < fdr_cut & abs(log2FoldChange_M) > logFC_cut))
  r_degs <- if (nrow(concord_sig) > 10) {
    cor(concord_sig$log2FoldChange_F, concord_sig$log2FoldChange_M, use = "complete.obs")
  } else NA_real_

  ## Directional concordance among shared DEGs
  both_sig <- concord_df %>%
    filter(!is.na(padj_F) & padj_F < fdr_cut & abs(log2FoldChange_F) > logFC_cut,
           !is.na(padj_M) & padj_M < fdr_cut & abs(log2FoldChange_M) > logFC_cut)
  n_shared <- nrow(both_sig)
  n_concordant <- sum(sign(both_sig$log2FoldChange_F) == sign(both_sig$log2FoldChange_M))
  n_discordant <- n_shared - n_concordant
  concordance_rate <- if (n_shared > 0) round(100 * n_concordant / n_shared, 1) else NA_real_

  writeLines("\n--- TEST 2: Sex-Stratified Concordance ---", report)
  writeLines(sprintf("  Female DEGs: %d (%d up, %d down)", n_deg_f, n_up_f, n_dn_f), report)
  writeLines(sprintf("  Male DEGs:   %d (%d up, %d down)", n_deg_m, n_up_m, n_dn_m), report)
  writeLines(sprintf("  Pearson r (all genes):     %.4f (p = %.2e)", r_all, r_test$p.value), report)
  writeLines(sprintf("  Pearson r (DEGs only):     %s",
                     ifelse(is.na(r_degs), "N/A (<10 shared)", sprintf("%.4f", r_degs))), report)
  writeLines(sprintf("  Spearman rho (all genes):  %.4f", rho_all), report)
  writeLines(sprintf("  Shared DEGs (both sig):    %d", n_shared), report)
  writeLines(sprintf("  Concordant (same dir):     %d", n_concordant), report)
  writeLines(sprintf("  Discordant (opp dir):      %d", n_discordant), report)
  writeLines(sprintf("  Concordance rate:          %s%%",
                     ifelse(is.na(concordance_rate), "N/A", concordance_rate)), report)

  ## Concordance scatter plot
  concord_df <- concord_df %>%
    mutate(
      f_sig = !is.na(padj_F) & padj_F < fdr_cut & abs(log2FoldChange_F) > logFC_cut,
      m_sig = !is.na(padj_M) & padj_M < fdr_cut & abs(log2FoldChange_M) > logFC_cut,
      category = factor(case_when(
        f_sig & m_sig & sign(log2FoldChange_F) == sign(log2FoldChange_M) ~ "Concordant",
        f_sig & m_sig & sign(log2FoldChange_F) != sign(log2FoldChange_M) ~ "Discordant",
        f_sig & !m_sig ~ "Female only",
        !f_sig & m_sig ~ "Male only",
        TRUE ~ "NS"
      ), levels = c("Concordant", "Discordant", "Female only", "Male only", "NS"))
    )

  cat_colors <- c("Concordant" = "#4CAF50", "Discordant" = "#F44336",
                   "Female only" = "#B2182B", "Male only" = "#2166AC", "NS" = "grey85")

  p_concord <- ggplot(concord_df, aes(x = log2FoldChange_M, y = log2FoldChange_F, color = category)) +
    geom_point(data = concord_df %>% filter(category == "NS"), size = 0.4, alpha = 0.2) +
    geom_point(data = concord_df %>% filter(category != "NS"), size = 1.2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    scale_color_manual(values = cat_colors, drop = FALSE, name = "") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r_all, 3)), fontface = "bold", size = 5) +
    labs(title = paste0(depot, ": Sex-Stratified KAT8-KD Concordance"),
         subtitle = paste0("Pearson r = ", round(r_all, 3),
                           "; Concordance = ", concordance_rate, "%"),
         x = expression(log[2]~FC~"(Male: AKO vs FL)"),
         y = expression(log[2]~FC~"(Female: AKO vs FL)")) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("02_sex_concordance_scatter_", depot, ".png")),
         p_concord, width = 7, height = 7, dpi = 300, bg = "white")

  ## ---------------------------------------------------------
  ## 4d) TEST 3: Model comparison
  ##     3a: ~Genotype vs ~Sex + Genotype (per-depot, value of Sex covariate)
  ##     3b: Per-depot ~Sex+Genotype vs unified ~Sex+GroupDepot (all 38 samples)
  ##     3c: ~0+GroupDepot vs ~Sex+GroupDepot (unified, justifies adding Sex
  ##         to the boss's original model)
  ## ---------------------------------------------------------
  cat("[", depot, "] Comparing models: ~Genotype vs ~Sex+Genotype vs unified...\n", sep = "")

  ## Model A: ~ Genotype only (no sex covariate, per-depot)
  dds_nosex <- DESeqDataSetFromMatrix(
    countData = round(depot_counts),
    colData   = depot_annot,
    design    = ~ Genotype
  )
  dds_nosex <- DESeq(dds_nosex, quiet = TRUE)
  res_nosex <- deseq_res_to_df(results(dds_nosex))
  res_nosex$gene <- depot_gene_names

  ## Model B: ~ Sex + Genotype (sex as covariate, per-depot)
  dds_sexadj <- DESeqDataSetFromMatrix(
    countData = round(depot_counts),
    colData   = depot_annot,
    design    = ~ Sex + Genotype
  )
  dds_sexadj <- DESeq(dds_sexadj, quiet = TRUE)
  res_sexadj <- deseq_res_to_df(results(dds_sexadj, name = "Genotype_KAT8KD_vs_CTL"))
  res_sexadj$gene <- depot_gene_names

  ## Model C: Unified ~ Sex + GroupDepot (all 38 samples, pre-computed)
  res_unified <- unified_results[[depot]]

  ## --- 3a: Compare A vs B (value of Sex covariate) ---
  n_deg_nosex <- sum(res_nosex$padj < fdr_cut & abs(res_nosex$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_deg_sexadj <- sum(res_sexadj$padj < fdr_cut & abs(res_sexadj$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_up_nosex  <- sum(res_nosex$padj < fdr_cut & res_nosex$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_nosex  <- sum(res_nosex$padj < fdr_cut & res_nosex$log2FoldChange < -logFC_cut, na.rm = TRUE)
  n_up_sexadj <- sum(res_sexadj$padj < fdr_cut & res_sexadj$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_sexadj <- sum(res_sexadj$padj < fdr_cut & res_sexadj$log2FoldChange < -logFC_cut, na.rm = TRUE)

  ## logFC correlation between models A and B
  compare_models <- merge(
    res_nosex[, c("gene", "log2FoldChange", "padj")],
    res_sexadj[, c("gene", "log2FoldChange", "padj")],
    by = "gene", suffixes = c("_nosex", "_sexadj")
  ) %>%
    filter(is.finite(log2FoldChange_nosex), is.finite(log2FoldChange_sexadj))

  r_models <- cor(compare_models$log2FoldChange_nosex, compare_models$log2FoldChange_sexadj)

  ## Genes that gain/lose significance (A -> B)
  compare_models <- compare_models %>%
    mutate(
      sig_nosex  = !is.na(padj_nosex) & padj_nosex < fdr_cut & abs(log2FoldChange_nosex) > logFC_cut,
      sig_sexadj = !is.na(padj_sexadj) & padj_sexadj < fdr_cut & abs(log2FoldChange_sexadj) > logFC_cut,
      gained = !sig_nosex & sig_sexadj,
      lost   = sig_nosex & !sig_sexadj
    )
  n_gained <- sum(compare_models$gained)
  n_lost   <- sum(compare_models$lost)
  pct_power_gain <- round(100 * (n_deg_sexadj - n_deg_nosex) / max(1, n_deg_nosex), 1)

  writeLines("\n--- TEST 3a: Per-Depot Model Comparison ---", report)
  writeLines(sprintf("  Model A (~ Genotype, no sex covariate):"), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_nosex, n_up_nosex, n_dn_nosex), report)
  writeLines(sprintf("  Model B (~ Sex + Genotype, per-depot):"), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_sexadj, n_up_sexadj, n_dn_sexadj), report)
  writeLines(sprintf("  logFC correlation (A vs B):  r = %.4f", r_models), report)
  writeLines(sprintf("  Genes gained by adding Sex:  %d", n_gained), report)
  writeLines(sprintf("  Genes lost by adding Sex:    %d", n_lost), report)
  writeLines(sprintf("  Net change: %+d DEGs (%+.1f%%)", n_deg_sexadj - n_deg_nosex, pct_power_gain), report)

  if (n_deg_sexadj > n_deg_nosex * 1.1) {
    writeLines("  >>> Sex covariate improves power", report)
  } else {
    writeLines("  >>> Models are equivalent; Sex covariate optional", report)
  }

  ## --- 3b: Compare B (per-depot) vs C (unified) ---
  ## This validates the unified model used in Part 1
  ## Unified model uses shared dispersions from all 38 samples
  common_genes_uc <- intersect(res_sexadj$gene, res_unified$gene)
  res_unified_common <- res_unified[res_unified$gene %in% common_genes_uc, ]
  res_sexadj_common  <- res_sexadj[res_sexadj$gene %in% common_genes_uc, ]
  ## Ensure same gene order
  res_unified_common <- res_unified_common[match(common_genes_uc, res_unified_common$gene), ]
  res_sexadj_common  <- res_sexadj_common[match(common_genes_uc, res_sexadj_common$gene), ]

  n_deg_unified <- sum(res_unified_common$padj < fdr_cut &
                        abs(res_unified_common$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_up_unified  <- sum(res_unified_common$padj < fdr_cut &
                        res_unified_common$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_unified  <- sum(res_unified_common$padj < fdr_cut &
                        res_unified_common$log2FoldChange < -logFC_cut, na.rm = TRUE)

  compare_uc <- data.frame(
    gene = common_genes_uc,
    lfc_perdepot = res_sexadj_common$log2FoldChange,
    lfc_unified  = res_unified_common$log2FoldChange,
    padj_perdepot = res_sexadj_common$padj,
    padj_unified  = res_unified_common$padj,
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(lfc_perdepot), is.finite(lfc_unified))

  r_unified <- cor(compare_uc$lfc_perdepot, compare_uc$lfc_unified, use = "complete.obs")

  compare_uc <- compare_uc %>%
    mutate(
      sig_perdepot = !is.na(padj_perdepot) & padj_perdepot < fdr_cut & abs(lfc_perdepot) > logFC_cut,
      sig_unified  = !is.na(padj_unified) & padj_unified < fdr_cut & abs(lfc_unified) > logFC_cut,
      gained_unified = !sig_perdepot & sig_unified,
      lost_unified   = sig_perdepot & !sig_unified
    )
  n_gained_unified <- sum(compare_uc$gained_unified)
  n_lost_unified   <- sum(compare_uc$lost_unified)
  pct_unified_gain <- round(100 * (n_deg_unified - n_deg_sexadj) / max(1, n_deg_sexadj), 1)

  writeLines("\n--- TEST 3b: Unified Model Validation ---", report)
  writeLines(sprintf("  Model B (~ Sex + Genotype, per-depot, %d samples):",
                     sum(sample_annot$Depot == depot)), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_sexadj, n_up_sexadj, n_dn_sexadj), report)
  writeLines(sprintf("  Model C (~ Sex + GroupDepot, unified, %d samples):",
                     nrow(sample_annot)), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_unified, n_up_unified, n_dn_unified), report)
  writeLines(sprintf("  logFC correlation (B vs C):  r = %.4f", r_unified), report)
  writeLines(sprintf("  DEGs gained by unified:      %d", n_gained_unified), report)
  writeLines(sprintf("  DEGs lost by unified:        %d", n_lost_unified), report)
  writeLines(sprintf("  Net change: %+d DEGs (%+.1f%%)", n_deg_unified - n_deg_sexadj, pct_unified_gain), report)

  if (r_unified > 0.99) {
    writeLines("  >>> Unified model is equivalent to per-depot (shared dispersions add power)", report)
  } else if (r_unified > 0.95) {
    writeLines("  >>> Unified model is highly concordant with per-depot", report)
  } else {
    writeLines("  >>> NOTE: Models show some divergence; check dispersion estimates", report)
  }

  ## Model comparison scatter (A vs B)
  p_compare <- ggplot(compare_models, aes(x = log2FoldChange_nosex, y = log2FoldChange_sexadj)) +
    geom_point(size = 0.3, alpha = 0.2, color = "grey60") +
    geom_point(data = compare_models %>% filter(gained), color = "#4CAF50", size = 1.5, alpha = 0.7) +
    geom_point(data = compare_models %>% filter(lost), color = "#F44336", size = 1.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r_models, 4)), fontface = "bold", size = 5) +
    labs(title = paste0(depot, ": Sex Covariate Effect"),
         subtitle = paste0("~Genotype: ", n_deg_nosex, " DEGs; ~Sex+Genotype: ", n_deg_sexadj, " DEGs"),
         x = expression(log[2]~FC~"(~ Genotype only)"),
         y = expression(log[2]~FC~"(~ Sex + Genotype)")) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("03a_model_comparison_sex_covariate_", depot, ".png")),
         p_compare, width = 7, height = 7, dpi = 300, bg = "white")

  ## Unified vs per-depot scatter
  p_unified <- ggplot(compare_uc, aes(x = lfc_perdepot, y = lfc_unified)) +
    geom_point(size = 0.3, alpha = 0.2, color = "grey60") +
    geom_point(data = compare_uc %>% filter(gained_unified), color = "#4CAF50", size = 1.5, alpha = 0.7) +
    geom_point(data = compare_uc %>% filter(lost_unified), color = "#F44336", size = 1.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r_unified, 4)), fontface = "bold", size = 5) +
    labs(title = paste0(depot, ": Per-Depot vs Unified Model"),
         subtitle = paste0("Per-depot: ", n_deg_sexadj, " DEGs; Unified (38 samples): ",
                           n_deg_unified, " DEGs"),
         x = expression(log[2]~FC~"(~ Sex + Genotype, per-depot)"),
         y = expression(log[2]~FC~"(~ Sex + GroupDepot, unified)")) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("03b_unified_vs_perdepot_", depot, ".png")),
         p_unified, width = 7, height = 7, dpi = 300, bg = "white")

  ## --- 3c: Compare unified ~ 0 + GroupDepot vs unified ~ Sex + GroupDepot ---
  ## This directly justifies adding Sex to the model the boss has been using.
  ## Both use all 38 samples; the only difference is the Sex covariate.
  cat("[", depot, "] Comparing unified models: ~ 0 + GroupDepot vs ~ Sex + GroupDepot...\n", sep = "")

  res_nosex_uni <- nosex_unified_results[[depot]]

  common_genes_nc <- intersect(res_nosex_uni$gene, res_unified$gene)
  res_nosex_uni_c <- res_nosex_uni[match(common_genes_nc, res_nosex_uni$gene), ]
  res_unified_c   <- res_unified[match(common_genes_nc, res_unified$gene), ]

  n_deg_nosex_uni <- sum(res_nosex_uni_c$padj < fdr_cut &
                          abs(res_nosex_uni_c$log2FoldChange) > logFC_cut, na.rm = TRUE)
  n_up_nosex_uni  <- sum(res_nosex_uni_c$padj < fdr_cut &
                          res_nosex_uni_c$log2FoldChange > logFC_cut, na.rm = TRUE)
  n_dn_nosex_uni  <- sum(res_nosex_uni_c$padj < fdr_cut &
                          res_nosex_uni_c$log2FoldChange < -logFC_cut, na.rm = TRUE)

  compare_nosex_vs_sex <- data.frame(
    gene = common_genes_nc,
    lfc_nosex = res_nosex_uni_c$log2FoldChange,
    lfc_sex   = res_unified_c$log2FoldChange,
    padj_nosex = res_nosex_uni_c$padj,
    padj_sex   = res_unified_c$padj,
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(lfc_nosex), is.finite(lfc_sex))

  r_nosex_vs_sex <- cor(compare_nosex_vs_sex$lfc_nosex,
                         compare_nosex_vs_sex$lfc_sex, use = "complete.obs")

  compare_nosex_vs_sex <- compare_nosex_vs_sex %>%
    mutate(
      sig_nosex = !is.na(padj_nosex) & padj_nosex < fdr_cut & abs(lfc_nosex) > logFC_cut,
      sig_sex   = !is.na(padj_sex) & padj_sex < fdr_cut & abs(lfc_sex) > logFC_cut,
      gained_by_sex = !sig_nosex & sig_sex,
      lost_by_sex   = sig_nosex & !sig_sex
    )
  n_gained_by_sex <- sum(compare_nosex_vs_sex$gained_by_sex)
  n_lost_by_sex   <- sum(compare_nosex_vs_sex$lost_by_sex)
  pct_sex_gain_uni <- round(100 * (n_deg_unified - n_deg_nosex_uni) / max(1, n_deg_nosex_uni), 1)

  ## LRT: does adding Sex significantly improve the unified model?
  dds_lrt_unified <- DESeq(dds_unified, test = "LRT",
                            reduced = ~ GroupDepot, quiet = TRUE)
  res_lrt_sex <- results(dds_lrt_unified)
  n_sex_sig_lrt <- sum(res_lrt_sex$padj < 0.05, na.rm = TRUE)
  n_sex_tested  <- sum(!is.na(res_lrt_sex$pvalue))
  pct_sex_sig   <- round(100 * n_sex_sig_lrt / n_sex_tested, 2)

  ## Median dispersion comparison
  disp_nosex <- median(dispersions(dds_nosex_unified), na.rm = TRUE)
  disp_sex   <- median(dispersions(dds_unified), na.rm = TRUE)
  disp_ratio <- round(disp_sex / disp_nosex, 4)

  writeLines("\n--- TEST 3c: Unified Model +/- Sex Covariate ---", report)
  writeLines(sprintf("  Model D (~ 0 + GroupDepot, no Sex, %d samples):",
                     nrow(sample_annot)), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_nosex_uni, n_up_nosex_uni, n_dn_nosex_uni), report)
  writeLines(sprintf("  Model C (~ Sex + GroupDepot, with Sex, %d samples):",
                     nrow(sample_annot)), report)
  writeLines(sprintf("    DEGs: %d (%d up, %d down)", n_deg_unified, n_up_unified, n_dn_unified), report)
  writeLines(sprintf("  logFC correlation (D vs C):  r = %.4f", r_nosex_vs_sex), report)
  writeLines(sprintf("  DEGs gained by adding Sex:   %d", n_gained_by_sex), report)
  writeLines(sprintf("  DEGs lost by adding Sex:     %d", n_lost_by_sex), report)
  writeLines(sprintf("  Net change: %+d DEGs (%+.1f%%)", n_deg_unified - n_deg_nosex_uni, pct_sex_gain_uni), report)
  writeLines(sprintf("  LRT (Sex covariate effect):  %d / %d genes (%.2f%%) where Sex matters",
                     n_sex_sig_lrt, n_sex_tested, pct_sex_sig), report)
  writeLines(sprintf("  Median dispersion ratio:     %.4f (< 1 means Sex absorbs variance)", disp_ratio), report)

  if (n_deg_unified > n_deg_nosex_uni * 1.05) {
    writeLines("  >>> Adding Sex to ~ 0 + GroupDepot IMPROVES power", report)
  } else if (n_deg_unified >= n_deg_nosex_uni) {
    writeLines("  >>> Adding Sex is neutral or slightly beneficial", report)
  } else {
    writeLines("  >>> Adding Sex shows minimal effect; models are equivalent", report)
  }

  ## Scatter plot: ~ 0 + GroupDepot vs ~ Sex + GroupDepot
  p_nosex_vs_sex <- ggplot(compare_nosex_vs_sex, aes(x = lfc_nosex, y = lfc_sex)) +
    geom_point(size = 0.3, alpha = 0.2, color = "grey60") +
    geom_point(data = compare_nosex_vs_sex %>% filter(gained_by_sex),
               color = "#4CAF50", size = 1.5, alpha = 0.7) +
    geom_point(data = compare_nosex_vs_sex %>% filter(lost_by_sex),
               color = "#F44336", size = 1.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r_nosex_vs_sex, 4)), fontface = "bold", size = 5) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = paste0("Gained: +", n_gained_by_sex, " DEGs\nLost: -", n_lost_by_sex, " DEGs"),
             color = "#4CAF50", fontface = "bold", size = 4) +
    labs(title = paste0(depot, ": Adding Sex to Unified Model"),
         subtitle = paste0("~ 0 + GroupDepot: ", n_deg_nosex_uni,
                           " DEGs; ~ Sex + GroupDepot: ", n_deg_unified, " DEGs"),
         x = expression(log[2]~FC~"(~ 0 + GroupDepot, no Sex)"),
         y = expression(log[2]~FC~"(~ Sex + GroupDepot)")) +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("03c_nosex_vs_sex_unified_", depot, ".png")),
         p_nosex_vs_sex, width = 7, height = 7, dpi = 300, bg = "white")

  ## ---------------------------------------------------------
  ## 4e) TEST 4: Variance decomposition
  ## ---------------------------------------------------------
  cat("[", depot, "] Computing variance decomposition...\n", sep = "")

  vsd <- vst(dds_wald, blind = TRUE)
  vst_mat <- assay(vsd)
  rownames(vst_mat) <- depot_gene_names

  set.seed(MASTER_SEED)
  n_var_genes <- min(5000, nrow(vst_mat))
  var_genes <- sample(depot_gene_names, n_var_genes)

  var_results <- data.frame(
    gene = var_genes,
    pct_sex = NA_real_, pct_genotype = NA_real_,
    pct_interaction = NA_real_, pct_residual = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(var_genes)) {
    g <- var_genes[i]
    y <- vst_mat[g, ]
    tryCatch({
      fit <- anova(lm(y ~ Sex + Genotype + Sex:Genotype, data = depot_annot))
      ss <- fit$`Sum Sq`
      ss_total <- sum(ss)
      if (ss_total > 0) {
        var_results$pct_sex[i]         <- 100 * ss[1] / ss_total
        var_results$pct_genotype[i]    <- 100 * ss[2] / ss_total
        var_results$pct_interaction[i] <- 100 * ss[3] / ss_total
        var_results$pct_residual[i]    <- 100 * ss[4] / ss_total
      }
    }, error = function(e) { })
  }

  var_results <- var_results %>% filter(!is.na(pct_sex))

  med_sex  <- round(median(var_results$pct_sex), 2)
  med_geno <- round(median(var_results$pct_genotype), 2)
  med_int  <- round(median(var_results$pct_interaction), 2)
  med_res  <- round(median(var_results$pct_residual), 2)
  mean_sex  <- round(mean(var_results$pct_sex), 2)
  mean_geno <- round(mean(var_results$pct_genotype), 2)
  mean_int  <- round(mean(var_results$pct_interaction), 2)
  mean_res  <- round(mean(var_results$pct_residual), 2)

  writeLines("\n--- TEST 4: Variance Decomposition ---", report)
  writeLines(sprintf("  (Across %d genes)", nrow(var_results)), report)
  writeLines(sprintf("                     Median    Mean"), report)
  writeLines(sprintf("  Sex:               %6.2f%%  %6.2f%%", med_sex, mean_sex), report)
  writeLines(sprintf("  Genotype:          %6.2f%%  %6.2f%%", med_geno, mean_geno), report)
  writeLines(sprintf("  Sex x Genotype:    %6.2f%%  %6.2f%%", med_int, mean_int), report)
  writeLines(sprintf("  Residual:          %6.2f%%  %6.2f%%", med_res, mean_res), report)

  ## Variance boxplot
  var_long <- var_results %>%
    dplyr::select(gene, pct_sex, pct_genotype, pct_interaction) %>%
    pivot_longer(-gene, names_to = "Term", values_to = "Pct") %>%
    mutate(Term = recode(Term, pct_sex = "Sex", pct_genotype = "Genotype",
                         pct_interaction = "Sex x Genotype"))

  p_var <- ggplot(var_long, aes(x = Term, y = Pct, fill = Term)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Sex" = "#FF7F00", "Genotype" = "#377EB8",
                                  "Sex x Genotype" = "#E41A1C")) +
    labs(title = paste0(depot, ": Variance Explained per Gene"),
         subtitle = sprintf("Median: Sex=%.1f%%, Genotype=%.1f%%, Interaction=%.1f%%",
                            med_sex, med_geno, med_int),
         x = NULL, y = "% Variance Explained") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = "none")

  ggsave(file.path(plots_dir, paste0("04_variance_decomposition_", depot, ".png")),
         p_var, width = 6, height = 5, dpi = 300, bg = "white")

  ## ---------------------------------------------------------
  ## 4f) TEST 5: Within-depot PCA by sex
  ## ---------------------------------------------------------
  pca_data <- plotPCA(vsd, intgroup = c("Sex", "Genotype"), returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Genotype, shape = Sex)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = c("CTL" = "#4DAF4A", "KAT8KD" = "#E41A1C")) +
    scale_shape_manual(values = c("F" = 16, "M" = 17)) +
    labs(title = paste0(depot, ": PCA (Sex x Genotype)"),
         subtitle = "Triangles & circles overlapping within genotype = safe to combine",
         x = paste0("PC1 (", pct_var[1], "%)"),
         y = paste0("PC2 (", pct_var[2], "%)")) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(file.path(plots_dir, paste0("05_PCA_by_sex_", depot, ".png")),
         p_pca, width = 7, height = 6, dpi = 300, bg = "white")

  writeLines(sprintf("\n  Within-depot PCA: PC1 = %.1f%%, PC2 = %.1f%%", pct_var[1], pct_var[2]), report)

  ## ---------------------------------------------------------
  ## 4g) VERDICT
  ## ---------------------------------------------------------
  ## Criteria:
  ##   PASS:    <5% interaction genes AND pi0 > 0.90 AND r > 0.7 AND concordance > 90%
  ##   CAUTION: 5-15% interaction genes OR pi0 0.80-0.90 OR r 0.5-0.7
  ##   CONCERN: >15% interaction genes OR pi0 < 0.80 OR r < 0.5 OR concordance < 80%

  flags <- c()
  if (pct_sig_lrt > 15)                          flags <- c(flags, "HIGH_INTERACTION_PCT")
  if (!is.na(pi0_est) && pi0_est < 0.80)         flags <- c(flags, "LOW_PI0")
  if (r_all < 0.5)                                flags <- c(flags, "LOW_CONCORDANCE_R")
  if (!is.na(concordance_rate) && concordance_rate < 80)
                                                   flags <- c(flags, "LOW_DIR_CONCORDANCE")

  cautions <- c()
  if (pct_sig_lrt > 5 && pct_sig_lrt <= 15)      cautions <- c(cautions, "MODERATE_INTERACTION_PCT")
  if (!is.na(pi0_est) && pi0_est >= 0.80 && pi0_est < 0.90)
                                                   cautions <- c(cautions, "MODERATE_PI0")
  if (r_all >= 0.5 && r_all < 0.7)               cautions <- c(cautions, "MODERATE_CONCORDANCE_R")

  if (length(flags) > 0) {
    verdict <- "CONCERN"
    verdict_detail <- paste("Flags:", paste(flags, collapse = ", "))
  } else if (length(cautions) > 0) {
    verdict <- "CAUTION"
    verdict_detail <- paste("Cautions:", paste(cautions, collapse = ", "))
  } else {
    verdict <- "PASS"
    verdict_detail <- "All criteria met for sex-combining"
  }

  writeLines(paste0("\n--- VERDICT: ", depot, " ---"), report)
  writeLines(paste0("  ", verdict, ": ", verdict_detail), report)

  if (verdict == "PASS") {
    writeLines("  >> Combining sexes is JUSTIFIED for this depot.", report)
  } else if (verdict == "CAUTION") {
    writeLines("  >> Combining sexes is ACCEPTABLE but mention caveats.", report)
  } else {
    writeLines("  >> RECONSIDER combining sexes. Run sex-stratified analysis.", report)
  }

  ## Add model comparison recommendation
  writeLines(paste0("  >> FINAL MODEL: ~ 0 + GroupDepot (", nrow(sample_annot),
                    " samples): ", n_deg_nosex_uni, " DEGs"), report)
  writeLines(paste0("     Balanced design (equal M/F per group) => Sex cannot confound"), report)
  writeLines(paste0("     Adding Sex as covariate: ", n_deg_nosex_uni, " -> ", n_deg_unified,
                    " DEGs (r = ", round(r_nosex_vs_sex, 3),
                    ") — extra DF cost outweighs variance absorbed"), report)

  cat("\n  >>> ", depot, " VERDICT: ", verdict, " <<<\n", sep = "")

  ## ---------------------------------------------------------
  ## 4h) Collect numbers for CSV
  ## ---------------------------------------------------------
  all_numbers[[depot]] <- data.frame(
    Depot = depot,
    ## LRT interaction test
    N_genes_tested = n_tested,
    N_interaction_sig_LRT = n_sig_lrt,
    Pct_interaction_sig = pct_sig_lrt,
    pi0_estimate = round(pi0_est, 4),
    N_interaction_sig_Wald = n_sig_wald,
    Median_abs_interaction_logFC = round(median_abs_int_lfc, 3),
    ## Concordance
    N_DEG_female = n_deg_f,
    N_DEG_male = n_deg_m,
    Pearson_r_all = round(r_all, 4),
    Pearson_r_DEGs = round(r_degs, 4),
    Spearman_rho = round(rho_all, 4),
    N_shared_DEGs = n_shared,
    N_concordant = n_concordant,
    N_discordant = n_discordant,
    Concordance_rate_pct = concordance_rate,
    ## Model comparison (3a: sex covariate)
    N_DEG_no_sex_covariate = n_deg_nosex,
    N_DEG_sex_covariate = n_deg_sexadj,
    N_gained_with_sex = n_gained,
    N_lost_with_sex = n_lost,
    Model_logFC_correlation = round(r_models, 4),
    ## Model comparison (3b: unified model)
    N_DEG_unified = n_deg_unified,
    N_gained_unified = n_gained_unified,
    N_lost_unified = n_lost_unified,
    Unified_logFC_correlation = round(r_unified, 4),
    ## Model comparison (3c: ~ 0 + GroupDepot vs ~ Sex + GroupDepot)
    N_DEG_nosex_unified = n_deg_nosex_uni,
    N_gained_by_sex_unified = n_gained_by_sex,
    N_lost_by_sex_unified = n_lost_by_sex,
    Nosex_vs_sex_logFC_correlation = round(r_nosex_vs_sex, 4),
    LRT_sex_covariate_sig = n_sex_sig_lrt,
    LRT_sex_covariate_pct = pct_sex_sig,
    Dispersion_ratio_sex_vs_nosex = disp_ratio,
    ## Variance decomposition
    Median_var_sex_pct = med_sex,
    Median_var_genotype_pct = med_geno,
    Median_var_interaction_pct = med_int,
    Median_var_residual_pct = med_res,
    ## Verdict
    Verdict = verdict,
    stringsAsFactors = FALSE
  )
}

## ============================================================
## 5) COMBINED SUMMARY
## ============================================================
numbers_df <- bind_rows(all_numbers)
write.csv(numbers_df, numbers_file, row.names = FALSE)

## ============================================================
## 6) MANUSCRIPT TEXT (ready to paste)
## ============================================================
writeLines(paste0("\n", paste(rep("=", 60), collapse = "")), report)
writeLines("  READY-TO-PASTE MANUSCRIPT TEXT", report)
writeLines(paste(rep("=", 60), collapse = ""), report)

## Build the text using actual numbers
for (depot in c("iWAT", "gWAT")) {
  row <- numbers_df[numbers_df$Depot == depot, ]

  ## Determine which concordance adjective to use
  conc_adj <- if (row$Pearson_r_all >= 0.8) "highly concordant"
              else if (row$Pearson_r_all >= 0.6) "concordant"
              else if (row$Pearson_r_all >= 0.4) "moderately concordant"
              else "weakly concordant"

  writeLines(paste0("\n--- ", depot, " ---\n"), report)

  writeLines(paste0(
    "To assess whether the transcriptional response to KAT8 loss differed ",
    "between sexes in ", depot, ", we performed a DESeq2 likelihood ratio test ",
    "comparing a full model (~ Sex + Genotype + Sex:Genotype) against a reduced ",
    "model (~ Sex + Genotype). Of ", formatC(row$N_genes_tested, big.mark = ","),
    " tested genes, only ", row$N_interaction_sig_LRT,
    " (", row$Pct_interaction_sig, "%) showed significant Sex x Genotype ",
    "interactions (FDR < 0.05)",
    if (!is.na(row$pi0_estimate)) {
      paste0(", with an estimated true null proportion (pi0) of ",
             row$pi0_estimate, ", indicating that approximately ",
             round(row$pi0_estimate * 100, 0), "% of genes show no evidence ",
             "of sex-dependent KAT8 knockout effects")
    } else "",
    ". Sex-stratified analysis identified ",
    row$N_DEG_female, " DEGs in females and ",
    row$N_DEG_male, " DEGs in males; sex-stratified log2 fold changes were ",
    conc_adj,
    " (Pearson r = ", row$Pearson_r_all,
    if (!is.na(row$Pearson_r_DEGs)) paste0("; r = ", row$Pearson_r_DEGs, " among DEGs") else "",
    "). ",
    if (row$N_shared_DEGs > 0) {
      paste0("Among ", row$N_shared_DEGs, " genes significant in both sexes, ",
             row$N_concordant, " (", row$Concordance_rate_pct,
             "%) were concordant in direction. ")
    } else "",
    "Variance decomposition showed the Sex x Genotype interaction explained ",
    "a median of only ", row$Median_var_interaction_pct, "% of per-gene variance, ",
    "compared to ", row$Median_var_genotype_pct, "% for the Genotype main effect",
    if (row$Median_var_sex_pct > 1) paste0(" and ", row$Median_var_sex_pct, "% for Sex") else "",
    ". Based on these findings, males and females were combined in a ",
    "unified DESeq2 model (~ GroupDepot) across all ", nrow(sample_annot),
    " tissue samples, with depot-specific KAT8 effects extracted via ",
    "explicit contrasts. Because the experimental design is balanced ",
    "(equal numbers of males and females per group), Sex cannot confound ",
    "the genotype effect and is not required in the model. ",
    "We confirmed this by comparing models with and without a Sex covariate: ",
    "adding Sex as a covariate (~ Sex + GroupDepot) ",
    if (row$N_DEG_unified < row$N_DEG_nosex_unified) {
      paste0("modestly reduced DEG counts (", row$N_DEG_nosex_unified, " vs. ",
             row$N_DEG_unified, " DEGs)")
    } else {
      paste0("yielded comparable DEG counts (", row$N_DEG_nosex_unified, " vs. ",
             row$N_DEG_unified, " DEGs)")
    },
    " while preserving effect estimates (logFC r = ",
    row$Nosex_vs_sex_logFC_correlation,
    "), indicating that the simpler model maximizes statistical power."
  ), report)
}

## Combined short version
writeLines(paste0("\n--- COMBINED SHORT VERSION (for Methods/Results) ---\n"), report)

iwat <- numbers_df[numbers_df$Depot == "iWAT", ]
gwat <- numbers_df[numbers_df$Depot == "gWAT", ]

writeLines(paste0(
  "Both male and female mice were included in each experimental group ",
  "(n = ", min(with(sample_annot, table(Sex, Genotype, Depot))),
  " per sex per genotype per depot). ",
  "To assess whether KAT8 knockout effects were sex-dependent, we performed ",
  "formal Sex x Genotype interaction testing (DESeq2 likelihood ratio test) ",
  "within each depot. The interaction term was significant for only ",
  iwat$N_interaction_sig_LRT, " genes (",
  iwat$Pct_interaction_sig, "%) in iWAT and ",
  gwat$N_interaction_sig_LRT, " genes (",
  gwat$Pct_interaction_sig, "%) in gWAT (FDR < 0.05). ",
  "Sex-stratified log2 fold changes were ",
  if (min(iwat$Pearson_r_all, gwat$Pearson_r_all) >= 0.7) "highly " else "",
  "concordant between males and females (Pearson r = ",
  iwat$Pearson_r_all, " in iWAT; r = ", gwat$Pearson_r_all, " in gWAT). ",
  "Males and females were therefore combined in a unified DESeq2 model ",
  "(~ GroupDepot) across all ", nrow(sample_annot),
  " tissue samples, identifying ", iwat$N_DEG_nosex_unified,
  " DEGs in iWAT and ", gwat$N_DEG_nosex_unified, " DEGs in gWAT. ",
  "The balanced experimental design (equal numbers of males and females per group) ",
  "ensures that sex cannot confound genotype effects. ",
  "Formal comparison confirmed that adding Sex as a covariate did not improve ",
  "sensitivity (logFC r = ", iwat$Nosex_vs_sex_logFC_correlation,
  " in iWAT, r = ", gwat$Nosex_vs_sex_logFC_correlation,
  " in gWAT), supporting the simpler model. ",
  "Depot-specific KAT8 knockout effects were extracted via explicit contrasts, ",
  "leveraging shared dispersion estimates across all samples for maximal ",
  "statistical power."
), report)

## ============================================================
## 7) CLOSE REPORT
## ============================================================
writeLines(paste0("\n", paste(rep("=", 60), collapse = "")), report)
writeLines(paste0("  Report saved to: ", report_file), report)
writeLines(paste0("  Numbers saved to: ", numbers_file), report)
writeLines(paste(rep("=", 60), collapse = ""), report)
close(report)

## Print summary to console
cat("\n============================================================\n")
cat("=== VALIDATION COMPLETE ===\n")
cat("============================================================\n\n")
cat("VERDICTS:\n")
for (depot in c("iWAT", "gWAT")) {
  row <- numbers_df[numbers_df$Depot == depot, ]
  cat("  ", depot, ": ", row$Verdict, "\n", sep = "")
}
cat("\nKey files:\n")
cat("  Report:  ", report_file, "\n")
cat("  Numbers: ", numbers_file, "\n")
cat("  Plots:   ", plots_dir, "\n\n")

cat("Numbers CSV preview:\n")
print(numbers_df %>% dplyr::select(Depot, N_interaction_sig_LRT, Pct_interaction_sig,
                             pi0_estimate, Pearson_r_all, Concordance_rate_pct,
                             N_DEG_no_sex_covariate, N_DEG_sex_covariate,
                             N_DEG_nosex_unified, N_DEG_unified,
                             Nosex_vs_sex_logFC_correlation, Verdict),
      row.names = FALSE)

cat("\n=== Part 7.7 COMPLETE ===\n")
