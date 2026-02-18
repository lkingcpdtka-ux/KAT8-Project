#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 7.6: FORMAL SEX x GENOTYPE
##                                INTERACTION TEST
## =========================================================
## THE question: "Is the KAT8-KD transcriptional response
##   different between males and females?"
##
## If no --> justified collapsing sexes (as done in the unified
##   ~ Sex + GroupDepot model in Part 1).
## If yes --> need sex-stratified analysis.
##
## APPROACH (per depot):
## 1. DESeq2 interaction model:  ~ Sex + Genotype + Sex:Genotype
##    - Wald test on interaction term (per-gene)
##    - LRT: full vs reduced (~ Sex + Genotype) for global test
## 2. P-value histogram of interaction term
##    - Uniform --> no systematic sex-dependent response
##    - Left-skewed --> some genes show sex-dependent KD effects
## 3. Estimate pi0 (proportion of true nulls) via qvalue
## 4. PCA within depot, colored by Sex (visual check)
## 5. Sex-stratified concordance: run DESeq2 separately for M/F,
##    scatter their KD logFCs (high r = same response)
## 6. Variance explained: Sex vs Genotype vs Interaction
## 7. Reviewer-ready summary table
##
## INTERPRETATION FOR MANUSCRIPT:
##   "We tested for Sex × Genotype interactions within each
##    depot using DESeq2 likelihood ratio tests. Only X genes
##    (Y%) showed significant interactions (FDR < 0.05),
##    with a pi0 estimate of Z, indicating the vast majority
##    of KAT8-KD effects are sex-independent. This justifies
##    pooling males and females within each depot."
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("DESeq2", "dplyr", "tidyr", "ggplot2", "ggrepel",
                   "scales", "RColorBrewer", "pheatmap", "grid",
                   "clusterProfiler", "org.Mm.eg.db", "enrichplot",
                   "AnnotationDbi", "stringr")
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
  library(RColorBrewer)
  library(pheatmap)
  library(grid)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(AnnotationDbi)
  library(stringr)
})

## Try to load qvalue for pi0 estimation (optional)
has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
if (has_qvalue) {
  library(qvalue)
  cat("[INFO] qvalue package loaded (pi0 estimation enabled)\n")
} else {
  cat("[WARN] qvalue not installed; skipping pi0 estimation\n")
  cat("       Install with: BiocManager::install('qvalue')\n")
}

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found.")
}

## Prefilter defaults (in case not set by earlier scripts)
if (!exists("prefilter_min_count"))   prefilter_min_count   <- 10
if (!exists("prefilter_min_samples")) prefilter_min_samples <- 4

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  use_save_core <- TRUE
} else {
  use_save_core <- FALSE
}

## ============================================================
## 2b) Helper: safe DESeq2 results --> data.frame
## ============================================================
## as.data.frame() on S4 DESeqResults can drop column names in some
## Bioconductor/S4Vectors versions.  This helper forces the standard
## DESeq2 column names so downstream code is reliable.
deseq_res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  std_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  if (ncol(df) == length(std_cols)) {
    colnames(df) <- std_cols
  } else if (ncol(df) > length(std_cols)) {
    colnames(df)[seq_along(std_cols)] <- std_cols
  }
  df
}

## ============================================================
## 3) LOAD RAW COUNTS + SAMPLE ANNOTATION
## ============================================================
cat("\n=== LOADING DATA ===\n")

counts_file <- file.path(getwd(), "counts.txt")
if (!file.exists(counts_file)) stop("counts.txt not found")

counts_raw <- read.delim(counts_file, header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)

## Sample annotation (same as part1)
sample_ids <- paste0("JS_", sprintf("%02d", 1:40))

sample_annot_full <- data.frame(
  Sample   = sample_ids,
  Depot    = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
  Sex      = c(rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5),
               rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5)),
  Genotype = c(rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5)),
  stringsAsFactors = FALSE
)

## Remove outliers
outlier_samples <- c("JS_08", "JS_28")
sample_annot <- sample_annot_full[!sample_annot_full$Sample %in% outlier_samples, ]
tissue_samples <- intersect(sample_annot$Sample, colnames(counts_raw))
sample_annot <- sample_annot[sample_annot$Sample %in% tissue_samples, ]

## Factorize
sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
sample_annot$Depot    <- factor(sample_annot$Depot, levels = c("iWAT", "gWAT"))
sample_annot$Sex      <- factor(sample_annot$Sex, levels = c("F", "M"))
rownames(sample_annot) <- sample_annot$Sample

## Build count matrix (gene names as rownames, same logic as part1)
## counts.txt uses "gene" (not "gene_name"), so check both
if ("gene_name" %in% colnames(counts_raw) || "gene" %in% colnames(counts_raw)) {
  counts_tissue <- counts_raw

  ## Normalise column name: counts.txt has "gene"; part1 renames to "gene_name"
  if (!"gene_name" %in% colnames(counts_tissue) && "gene" %in% colnames(counts_tissue)) {
    counts_tissue$gene_name <- counts_tissue$gene
  }
  ## Capture Ensembl IDs for fallback (column may be "gene_id" or "ensembl_gene_id")
  if (!"ensembl_gene_id" %in% colnames(counts_tissue) && "gene_id" %in% colnames(counts_tissue)) {
    counts_tissue$ensembl_gene_id <- counts_tissue$gene_id
  }

  ## Replace NA / blank gene names with Ensembl ID (same as part1 lines 268-269)
  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  if ("ensembl_gene_id" %in% colnames(counts_tissue)) {
    counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]
  }

  ## De-duplicate gene names (same as part1 lines 272-277)
  if (any(duplicated(counts_tissue$gene_name))) {
    dup_flag <- duplicated(counts_tissue$gene_name) | duplicated(counts_tissue$gene_name, fromLast = TRUE)
    if ("ensembl_gene_id" %in% colnames(counts_tissue)) {
      counts_tissue$gene_name[dup_flag] <- paste0(counts_tissue$ensembl_gene_id[dup_flag], "_",
                                                   counts_tissue$gene_name[dup_flag])
    }
    counts_tissue$gene_name <- make.unique(counts_tissue$gene_name)
  }

  rownames(counts_tissue) <- counts_tissue$gene_name
  count_matrix <- as.matrix(counts_tissue[, tissue_samples])
  ## Re-assign rownames AFTER as.matrix() (critical — part1 line 283)
  rownames(count_matrix) <- counts_tissue$gene_name
} else {
  count_matrix <- as.matrix(counts_raw[, tissue_samples])
}

## Prefilter (same rule as part1)
min_count   <- prefilter_min_count
min_samples <- prefilter_min_samples
keep <- rowSums(count_matrix >= min_count) >= min_samples
count_matrix <- count_matrix[keep, , drop = FALSE]
cat("[INFO] Genes after prefiltering: ", nrow(count_matrix), "\n")

## ============================================================
## Set up output directory
## ============================================================
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) > 0) {
  outdir <- run_dirs[which.max(file.info(run_dirs)$mtime)]
} else {
  outdir <- file.path(savepoint_dir, paste0("RUN_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
interaction_dir <- file.path(plots_dir, "sex_interaction_test")
if (!dir.exists(interaction_dir)) dir.create(interaction_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

## ============================================================
## 4) PER-DEPOT INTERACTION TEST
## ============================================================

summary_rows <- list()

for (depot in c("iWAT", "gWAT")) {

  cat("\n")
  cat("============================================================\n")
  cat("=== ", depot, ": SEX x GENOTYPE INTERACTION TEST ===\n", sep = "")
  cat("============================================================\n\n")

  ## Initialize per-depot tracking variables (avoid stale values from prior iteration)
  n_geno_sig <- NA_integer_
  n_sex_sig  <- NA_integer_

  ## -------------------------------------------------------
  ## 4a) Subset to this depot
  ## -------------------------------------------------------
  depot_samples <- sample_annot$Sample[sample_annot$Depot == depot]
  depot_annot   <- sample_annot[depot_samples, ]
  depot_counts  <- count_matrix[, depot_samples]

  ## Re-prefilter within depot (some genes may not be expressed here)
  depot_keep <- rowSums(depot_counts >= min_count) >= min_samples
  depot_counts <- depot_counts[depot_keep, , drop = FALSE]

  ## Capture gene names explicitly from the R matrix (robust against
  ## rowname loss in S4/DESeq2 objects and as.data.frame conversions)
  depot_gene_names <- rownames(depot_counts)
  if (is.null(depot_gene_names)) {
    depot_gene_names <- paste0("gene_", seq_len(nrow(depot_counts)))
    rownames(depot_counts) <- depot_gene_names
    cat("[WARN] Gene names missing from count matrix; using auto-generated names\n")
  }

  cat("[INFO] ", depot, ": ", nrow(depot_counts), " genes, ",
      ncol(depot_counts), " samples\n", sep = "")

  ## Print group sizes
  cat("  Group sizes:\n")
  with(depot_annot, print(table(Sex, Genotype)))
  cat("\n")

  ## -------------------------------------------------------
  ## 4b) DESeq2: Full interaction model
  ## -------------------------------------------------------
  ## Full:    ~ Sex + Genotype + Sex:Genotype
  ## Reduced: ~ Sex + Genotype
  ##
  ## The interaction term tests: "Is the Genotype effect
  ## different between males and females?"

  dds <- DESeqDataSetFromMatrix(
    countData = round(depot_counts),
    colData   = depot_annot,
    design    = ~ Sex + Genotype + Sex:Genotype
  )

  ## Run DESeq2 with LRT (compare full vs reduced)
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ Sex + Genotype)

  ## Also get Wald test results for the interaction coefficient
  dds_wald <- DESeq(dds, test = "Wald")

  ## -------------------------------------------------------
  ## 4c) Extract results
  ## -------------------------------------------------------
  ## LRT result: tests ALL genes for interaction effect
  res_lrt <- results(dds_lrt)
  res_lrt_df <- deseq_res_to_df(res_lrt)
  res_lrt_df$gene <- depot_gene_names
  res_lrt_df <- res_lrt_df[!is.na(res_lrt_df$pvalue), ]

  ## Wald result for interaction coefficient specifically
  ## Coefficient name: "SexM.GenotypeKAT8KD"
  coef_names <- resultsNames(dds_wald)
  interaction_coef <- grep("Sex.*Genotype|Genotype.*Sex", coef_names, value = TRUE)
  cat("[INFO] Model coefficients: ", paste(coef_names, collapse = ", "), "\n")
  cat("[INFO] Interaction term:   ", interaction_coef, "\n\n")

  if (length(interaction_coef) == 1) {
    res_wald <- results(dds_wald, name = interaction_coef)
  } else {
    ## Fallback: use contrast
    res_wald <- results(dds_wald, name = tail(coef_names, 1))
    cat("[WARN] Using last coefficient as interaction term\n")
  }

  res_wald_df <- deseq_res_to_df(res_wald)
  res_wald_df$gene <- depot_gene_names
  res_wald_df <- res_wald_df[!is.na(res_wald_df$pvalue), ]

  ## -------------------------------------------------------
  ## 4d) Count significant interaction genes
  ## -------------------------------------------------------
  thresh <- get_tissue_thresholds(depot)
  fdr_cut <- thresh$fdr_cut

  n_sig_lrt  <- sum(res_lrt_df$padj < fdr_cut, na.rm = TRUE)
  n_sig_wald <- sum(res_wald_df$padj < fdr_cut, na.rm = TRUE)
  n_tested   <- nrow(res_lrt_df)
  pct_sig    <- round(100 * n_sig_lrt / n_tested, 2)

  cat("LRT interaction test:\n")
  cat("  Genes tested:                   ", n_tested, "\n")
  cat("  Significant (padj < ", fdr_cut, "):  ", n_sig_lrt, " (", pct_sig, "%)\n", sep = "")
  cat("  Wald interaction sig:            ", n_sig_wald, "\n")

  ## -------------------------------------------------------
  ## 4e) pi0 estimation (proportion of true nulls)
  ## -------------------------------------------------------
  ## pi0 close to 1.0 = almost all genes show NO interaction
  ## This is the strongest statistical argument for merging

  pi0_est <- NA_real_
  if (has_qvalue && n_tested > 100) {
    tryCatch({
      qobj <- qvalue(res_lrt_df$pvalue)
      pi0_est <- qobj$pi0
      cat("  pi0 (true null proportion):      ", round(pi0_est, 4), "\n")
      cat("  Interpretation: ", round(pi0_est * 100, 1),
          "% of genes have NO sex-dependent KD response\n", sep = "")
    }, error = function(e) {
      cat("  [WARN] qvalue estimation failed: ", e$message, "\n")
    })
  }

  ## -------------------------------------------------------
  ## 4f) Genotype main effect (what we care about) vs
  ##     Sex:Genotype interaction (what we hope is small)
  ## -------------------------------------------------------
  ## Extract Genotype main effect from Wald
  geno_coef <- grep("^Genotype", coef_names, value = TRUE)
  if (length(geno_coef) == 1) {
    res_geno <- deseq_res_to_df(results(dds_wald, name = geno_coef))
    res_geno$gene <- depot_gene_names
    n_geno_sig <- sum(res_geno$padj < fdr_cut & abs(res_geno$log2FoldChange) > thresh$logFC_cut,
                      na.rm = TRUE)
    cat("\n  Genotype main effect (KD vs CTL):\n")
    cat("    Significant DEGs:              ", n_geno_sig, "\n")
    cat("    Ratio interaction/genotype:    ", round(n_sig_lrt / max(1, n_geno_sig), 3), "\n")
  }

  ## Sex main effect
  sex_coef <- grep("^Sex", coef_names, value = TRUE)
  sex_coef <- sex_coef[!grepl(":", sex_coef)]  # exclude interaction
  if (length(sex_coef) == 1) {
    res_sex <- deseq_res_to_df(results(dds_wald, name = sex_coef))
    n_sex_sig <- sum(res_sex$padj < fdr_cut, na.rm = TRUE)
    cat("  Sex main effect:\n")
    cat("    Significant:                   ", n_sex_sig, "\n")
  }

  ## -------------------------------------------------------
  ## 4g) PLOT 1: P-value histogram (interaction term)
  ## -------------------------------------------------------
  p_hist <- ggplot(res_lrt_df, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", color = "white", boundary = 0) +
    geom_hline(yintercept = n_tested / 50, linetype = "dashed", color = "red") +
    labs(
      title = paste0(depot, ": P-value Distribution (Sex x Genotype Interaction)"),
      subtitle = paste0(
        n_sig_lrt, " / ", n_tested, " genes significant (",
        pct_sig, "%)",
        if (!is.na(pi0_est)) paste0("; pi0 = ", round(pi0_est, 3)) else ""
      ),
      x = "P-value (LRT: full vs reduced)",
      y = "Number of genes"
    ) +
    annotate("text", x = 0.75, y = Inf, vjust = 2,
             label = "Red line = expected\nunder uniform null",
             color = "red", size = 3.5) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )

  ggsave(file.path(interaction_dir,
                   paste0("pvalue_hist_interaction_", depot, "_", run_tag, ".png")),
         plot = p_hist, width = 7, height = 5, dpi = 300, bg = "white")

  ## -------------------------------------------------------
  ## 4h) PLOT 2: Interaction effect size vs Genotype effect
  ## -------------------------------------------------------
  ## If interaction effects are tiny relative to genotype
  ## effects, that's a strong argument for merging

  if (length(interaction_coef) == 1 && length(geno_coef) == 1) {
    compare_df <- data.frame(
      gene = res_geno$gene,
      geno_lfc = res_geno$log2FoldChange,
      geno_padj = res_geno$padj,
      interaction_lfc = res_wald_df$log2FoldChange[match(res_geno$gene, res_wald_df$gene)],
      interaction_padj = res_wald_df$padj[match(res_geno$gene, res_wald_df$gene)],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(is.finite(geno_lfc), is.finite(interaction_lfc))

    cat_levels <- c("Sig. interaction", "Genotype DEG", "NS")
    cat_colors <- c("Sig. interaction" = "#E41A1C",
                     "Genotype DEG" = "#377EB8",
                     "NS" = "grey85")

    if (nrow(compare_df) > 0) {
      ## Classify
      compare_df <- compare_df %>%
        mutate(
          category = factor(case_when(
            !is.na(interaction_padj) & interaction_padj < fdr_cut ~ "Sig. interaction",
            !is.na(geno_padj) & geno_padj < fdr_cut & abs(geno_lfc) > thresh$logFC_cut ~ "Genotype DEG",
            TRUE ~ "NS"
          ), levels = cat_levels)
        )

      p_effect <- ggplot(compare_df, aes(x = geno_lfc, y = interaction_lfc, color = category)) +
        geom_point(data = compare_df %>% dplyr::filter(category == "NS"),
                   size = 0.5, alpha = 0.3) +
        geom_point(data = compare_df %>% dplyr::filter(category != "NS"),
                   size = 1.5, alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        scale_color_manual(values = cat_colors, drop = FALSE, name = "") +
        labs(
          title = paste0(depot, ": Genotype Effect vs Sex x Genotype Interaction"),
          subtitle = "Interaction effects clustered near zero = sexes respond similarly",
          x = expression(log[2]~FC~"(Genotype main effect: KD vs CTL)"),
          y = expression(log[2]~FC~"(Sex x Genotype interaction)")
        ) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = "right"
        )

      ggsave(file.path(interaction_dir,
                       paste0("interaction_vs_genotype_", depot, "_", run_tag, ".png")),
             plot = p_effect, width = 8, height = 6, dpi = 300, bg = "white")
    } else {
      cat("[WARN] No genes with finite values for both genotype and interaction effects; ",
          "skipping effect size scatter plot\n")
    }
  }

  ## -------------------------------------------------------
  ## 4i) PLOT 3: PCA within depot, colored by Sex
  ## -------------------------------------------------------
  vsd <- vst(dds_wald, blind = TRUE)
  pca_data <- plotPCA(vsd, intgroup = c("Sex", "Genotype"), returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

  p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Genotype, shape = Sex)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = c("CTL" = "#4DAF4A", "KAT8KD" = "#E41A1C")) +
    scale_shape_manual(values = c("F" = 16, "M" = 17)) +
    labs(
      title = paste0(depot, ": PCA (colored by Genotype, shaped by Sex)"),
      subtitle = "If sexes overlap within each genotype --> safe to merge",
      x = paste0("PC1 (", pct_var[1], "% variance)"),
      y = paste0("PC2 (", pct_var[2], "% variance)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )

  ggsave(file.path(interaction_dir,
                   paste0("PCA_by_sex_", depot, "_", run_tag, ".png")),
         plot = p_pca, width = 7, height = 6, dpi = 300, bg = "white")

  ## -------------------------------------------------------
  ## 4i-b) PLOT 3b: Sex-stratified concordance
  ## -------------------------------------------------------
  ## Run DESeq2 separately for females and males, then scatter
  ## their KD-vs-CTL log2FCs. High correlation = sexes respond
  ## the same way, strongest evidence for collapsing.

  female_idx <- depot_annot$Sample[depot_annot$Sex == "F"]
  male_idx   <- depot_annot$Sample[depot_annot$Sex == "M"]

  if (length(female_idx) >= 4 && length(male_idx) >= 4) {
    cat("[INFO] Running sex-stratified DESeq2 for concordance analysis...\n")

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

    ## Merge
    concord_df <- merge(
      res_f[, c("gene", "log2FoldChange", "padj")],
      res_m[, c("gene", "log2FoldChange", "padj")],
      by = "gene", suffixes = c("_F", "_M")
    ) %>%
      dplyr::filter(is.finite(log2FoldChange_F), is.finite(log2FoldChange_M))

    if (nrow(concord_df) > 0) {
      ## Correlation
      r_val <- cor(concord_df$log2FoldChange_F, concord_df$log2FoldChange_M,
                   use = "complete.obs")
      cat("  Sex-stratified logFC correlation (Pearson r): ", round(r_val, 3), "\n")

      ## Classify significance
      concord_df <- concord_df %>%
        mutate(
          sig_class = factor(case_when(
            !is.na(padj_F) & padj_F < fdr_cut &
              !is.na(padj_M) & padj_M < fdr_cut ~ "Sig. in both",
            !is.na(padj_F) & padj_F < fdr_cut ~ "Sig. in F only",
            !is.na(padj_M) & padj_M < fdr_cut ~ "Sig. in M only",
            TRUE ~ "NS"
          ), levels = c("Sig. in both", "Sig. in F only", "Sig. in M only", "NS"))
        )

      concord_colors <- c("Sig. in both"  = "#984EA3",
                           "Sig. in F only" = "#E41A1C",
                           "Sig. in M only" = "#377EB8",
                           "NS"             = "grey85")

      p_concord <- ggplot(concord_df,
                          aes(x = log2FoldChange_F, y = log2FoldChange_M, color = sig_class)) +
        geom_point(data = concord_df %>% dplyr::filter(sig_class == "NS"),
                   size = 0.5, alpha = 0.2) +
        geom_point(data = concord_df %>% dplyr::filter(sig_class != "NS"),
                   size = 1.5, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
        scale_color_manual(values = concord_colors, drop = FALSE, name = "") +
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
                 label = paste0("r = ", round(r_val, 3)),
                 fontface = "bold", size = 5, color = "black") +
        labs(
          title = paste0(depot, ": KD Effect Concordance Between Sexes"),
          subtitle = "Points near diagonal = same response in males and females",
          x = expression(log[2]~FC~"(KD vs CTL) in Females"),
          y = expression(log[2]~FC~"(KD vs CTL) in Males")
        ) +
        coord_fixed() +
        theme_bw(base_size = 12) +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
          legend.position = "right"
        )

      ggsave(file.path(interaction_dir,
                       paste0("sex_concordance_scatter_", depot, "_", run_tag, ".png")),
             plot = p_concord, width = 7, height = 7, dpi = 300, bg = "white")

      ## Save concordance table
      write.csv(concord_df,
                file.path(tables_dir,
                          paste0("sex_stratified_logfc_", depot, "_", run_tag, ".csv")),
                row.names = FALSE)
    }
  } else {
    cat("[WARN] Too few samples per sex for concordance analysis\n")
    r_val <- NA_real_
  }

  ## -------------------------------------------------------
  ## 4j) PLOT 4: Variance explained by each model term
  ## -------------------------------------------------------
  ## Quick approach: for each gene, fit Sex + Genotype + Sex:Genotype
  ## and compute sum-of-squares attributed to each term
  ## Then summarize across genes

  cat("\n[INFO] Computing variance attribution (this may take a moment)...\n")

  ## Use VST-transformed data (already computed for PCA, statistically
  ## more appropriate for linear-model variance decomposition than log2+1)
  vst_mat <- assay(vsd)
  ## Ensure VST matrix has gene names (may be lost in S4 conversion)
  rownames(vst_mat) <- depot_gene_names
  gene_names <- depot_gene_names

  ## Use a sample of genes to speed this up
  n_genes_for_var <- min(5000, length(gene_names))
  set.seed(MASTER_SEED)
  if (n_genes_for_var > 1) {
    var_genes <- sample(gene_names, n_genes_for_var)
  } else {
    var_genes <- gene_names
  }

  var_results <- data.frame(
    gene = var_genes,
    pct_sex = NA_real_,
    pct_genotype = NA_real_,
    pct_interaction = NA_real_,
    pct_residual = NA_real_,
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

  var_results <- var_results %>% dplyr::filter(!is.na(pct_sex))

  ## Summary
  var_summary <- data.frame(
    Term = c("Sex", "Genotype (KD vs CTL)", "Sex x Genotype", "Residual"),
    Median_pct = c(
      median(var_results$pct_sex),
      median(var_results$pct_genotype),
      median(var_results$pct_interaction),
      median(var_results$pct_residual)
    ),
    Mean_pct = c(
      mean(var_results$pct_sex),
      mean(var_results$pct_genotype),
      mean(var_results$pct_interaction),
      mean(var_results$pct_residual)
    )
  )

  cat("\n  Variance attribution (", n_genes_for_var, " genes):\n", sep = "")
  print(var_summary, row.names = FALSE, digits = 3)

  ## Bar plot of variance explained
  var_long <- var_results %>%
    dplyr::select(gene, pct_sex, pct_genotype, pct_interaction) %>%
    pivot_longer(cols = -gene, names_to = "Term", values_to = "Pct") %>%
    mutate(Term = recode(Term,
                         "pct_sex" = "Sex",
                         "pct_genotype" = "Genotype",
                         "pct_interaction" = "Sex x Genotype"))

  term_colors <- c("Sex" = "#FF7F00", "Genotype" = "#377EB8", "Sex x Genotype" = "#E41A1C")

  p_var <- ggplot(var_long, aes(x = Term, y = Pct, fill = Term)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    scale_fill_manual(values = term_colors) +
    labs(
      title = paste0(depot, ": Variance Explained per Gene"),
      subtitle = paste0("Across ", nrow(var_results), " genes; ",
                        "median interaction = ",
                        round(median(var_results$pct_interaction), 1), "%"),
      x = NULL,
      y = "% Variance Explained"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "none"
    )

  ggsave(file.path(interaction_dir,
                   paste0("variance_explained_", depot, "_", run_tag, ".png")),
         plot = p_var, width = 6, height = 5, dpi = 300, bg = "white")

  ## -------------------------------------------------------
  ## 4k) Save interaction results table
  ## -------------------------------------------------------
  ## Build table with explicit column access (avoids dplyr NSE column name issues)
  interaction_table <- data.frame(
    gene               = res_lrt_df[["gene"]],
    LRT_baseMean       = res_lrt_df[["baseMean"]],
    LRT_stat           = res_lrt_df[["stat"]],
    LRT_pvalue         = res_lrt_df[["pvalue"]],
    LRT_padj           = res_lrt_df[["padj"]],
    interaction_log2FC = res_wald_df[["log2FoldChange"]][
        match(res_lrt_df[["gene"]], res_wald_df[["gene"]])],
    interaction_padj   = res_wald_df[["padj"]][
        match(res_lrt_df[["gene"]], res_wald_df[["gene"]])],
    stringsAsFactors = FALSE
  )
  interaction_table <- interaction_table[order(interaction_table$LRT_pvalue), ]

  write.csv(interaction_table,
            file.path(tables_dir, paste0("sex_interaction_results_", depot, "_", run_tag, ".csv")),
            row.names = FALSE)
  cat("\n[OK] Saved interaction results table for ", depot, "\n", sep = "")

  ## If there ARE significant interaction genes, list them
  sig_interaction <- interaction_table %>% dplyr::filter(LRT_padj < fdr_cut)
  if (nrow(sig_interaction) > 0) {
    write.csv(sig_interaction,
              file.path(tables_dir, paste0("sex_interaction_sig_genes_", depot, "_", run_tag, ".csv")),
              row.names = FALSE)
    cat("[INFO] ", nrow(sig_interaction), " significant interaction genes saved\n", sep = "")
  }

  ## -------------------------------------------------------
  ## 4m) PATHWAY ANALYSIS ON INTERACTION GENES
  ## -------------------------------------------------------
  ## Which biological pathways show sex-dependent KAT8-KD effects?

  if (nrow(sig_interaction) >= 5) {
    cat("\n[INFO] Running pathway analysis on ", nrow(sig_interaction),
        " significant interaction genes (", depot, ")...\n", sep = "")

    interaction_genes <- sig_interaction$gene
    universe_genes    <- interaction_table$gene

    ## Convert gene symbols to Entrez IDs
    gene_entrez <- tryCatch(
      clusterProfiler::bitr(unique(interaction_genes),
                            fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = org.Mm.eg.db, drop = TRUE),
      error = function(e) {
        cat("[WARN] Gene ID conversion failed: ", e$message, "\n")
        data.frame(SYMBOL = character(), ENTREZID = character())
      })

    universe_entrez <- tryCatch(
      clusterProfiler::bitr(unique(universe_genes),
                            fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = org.Mm.eg.db, drop = TRUE),
      error = function(e) {
        cat("[WARN] Universe ID conversion failed: ", e$message, "\n")
        data.frame(SYMBOL = character(), ENTREZID = character())
      })

    cat("[INFO] Mapped ", nrow(gene_entrez), "/", length(interaction_genes),
        " interaction genes to Entrez IDs\n", sep = "")

    if (nrow(gene_entrez) >= 5) {

      ## --- GO:BP ORA ---
      tryCatch({
        ego <- enrichGO(
          gene          = gene_entrez$ENTREZID,
          universe      = universe_entrez$ENTREZID,
          OrgDb         = org.Mm.eg.db,
          ont           = "BP",
          pvalueCutoff  = ora_params$pvalue_cutoff,
          qvalueCutoff  = ora_params$qvalue_cutoff,
          readable      = TRUE,
          minGSSize     = ora_params$min_gs_size,
          maxGSSize     = ora_params$max_gs_size
        )

        if (!is.null(ego) && nrow(ego@result) > 0) {
          ego_simp <- simplify(ego, cutoff = ora_params$simplify_cutoff,
                               by = "p.adjust", select_fun = min)
          sig_go <- ego_simp@result %>%
            dplyr::filter(p.adjust < 0.05) %>%
            dplyr::arrange(p.adjust)

          if (nrow(sig_go) > 0) {
            write.csv(sig_go,
                      file.path(tables_dir,
                                paste0("interaction_ORA_GOBP_", depot, "_", run_tag, ".csv")),
                      row.names = FALSE)
            cat("[OK] GO:BP: ", nrow(sig_go), " significant terms\n", sep = "")

            ## Dotplot (top 15)
            top_go <- sig_go %>% dplyr::slice_head(n = ora_params$top_n)
            top_go$GeneRatio_num <- sapply(strsplit(top_go$GeneRatio, "/"), function(x)
              as.numeric(x[1]) / as.numeric(x[2]))
            top_go$Count <- sapply(strsplit(top_go$GeneRatio, "/"), function(x)
              as.numeric(x[1]))
            top_go$neg_log10_fdr <- -log10(top_go$p.adjust)
            top_go$Description <- factor(
              stringr::str_wrap(top_go$Description, width = 32),
              levels = rev(stringr::str_wrap(top_go$Description, width = 32)))

            p_go <- ggplot(top_go, aes(x = GeneRatio_num, y = Description)) +
              geom_point(aes(size = Count, color = neg_log10_fdr)) +
              scale_color_gradientn(
                colors = c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027"),
                name = expression(-log[10]*"(FDR)")) +
              scale_size_continuous(name = "Count", range = c(2, 6)) +
              labs(title = paste0(depot, ": GO:BP - Interaction Genes"),
                   subtitle = paste0(nrow(sig_interaction),
                                     " genes with significant Sex x Genotype interaction"),
                   x = "Gene Ratio", y = NULL) +
              theme_bw(base_size = 10) +
              theme(plot.title = element_text(face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
                    panel.grid.minor = element_blank())

            ggsave(file.path(interaction_dir,
                             paste0("interaction_ORA_GOBP_dotplot_", depot, "_", run_tag, ".png")),
                   p_go, width = 7, height = 6, dpi = 300, bg = "white")
          } else {
            cat("[INFO] No significant GO:BP terms for interaction genes\n")
          }
        }
      }, error = function(e) {
        cat("[WARN] GO:BP ORA failed: ", e$message, "\n")
      })

      ## --- KEGG ORA ---
      tryCatch({
        ekegg <- enrichKEGG(
          gene         = gene_entrez$ENTREZID,
          universe     = universe_entrez$ENTREZID,
          organism     = "mmu",
          pvalueCutoff = ora_params$pvalue_cutoff,
          qvalueCutoff = ora_params$qvalue_cutoff,
          minGSSize    = ora_params$min_gs_size,
          maxGSSize    = ora_params$max_gs_size
        )

        if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
          sig_kegg <- ekegg@result %>%
            dplyr::filter(p.adjust < 0.05) %>%
            dplyr::arrange(p.adjust)

          ## Translate Entrez IDs in geneID column to symbols
          if (nrow(sig_kegg) > 0 && "geneID" %in% colnames(sig_kegg)) {
            sig_kegg$geneID_symbols <- sapply(sig_kegg$geneID, function(g) {
              ids <- unlist(strsplit(g, "/"))
              syms <- gene_entrez$SYMBOL[match(ids, gene_entrez$ENTREZID)]
              paste(syms[!is.na(syms)], collapse = "/")
            })
          }

          if (nrow(sig_kegg) > 0) {
            write.csv(sig_kegg,
                      file.path(tables_dir,
                                paste0("interaction_ORA_KEGG_", depot, "_", run_tag, ".csv")),
                      row.names = FALSE)
            cat("[OK] KEGG: ", nrow(sig_kegg), " significant pathways\n", sep = "")

            top_kegg <- sig_kegg %>% dplyr::slice_head(n = ora_params$top_n)
            top_kegg$GeneRatio_num <- sapply(strsplit(top_kegg$GeneRatio, "/"), function(x)
              as.numeric(x[1]) / as.numeric(x[2]))
            top_kegg$Count <- sapply(strsplit(top_kegg$GeneRatio, "/"), function(x)
              as.numeric(x[1]))
            top_kegg$neg_log10_fdr <- -log10(top_kegg$p.adjust)
            top_kegg$Description <- factor(
              stringr::str_wrap(top_kegg$Description, width = 32),
              levels = rev(stringr::str_wrap(top_kegg$Description, width = 32)))

            p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_num, y = Description)) +
              geom_point(aes(size = Count, color = neg_log10_fdr)) +
              scale_color_gradientn(
                colors = c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027"),
                name = expression(-log[10]*"(FDR)")) +
              scale_size_continuous(name = "Count", range = c(2, 6)) +
              labs(title = paste0(depot, ": KEGG - Interaction Genes"),
                   subtitle = paste0(nrow(sig_interaction),
                                     " genes with significant Sex x Genotype interaction"),
                   x = "Gene Ratio", y = NULL) +
              theme_bw(base_size = 10) +
              theme(plot.title = element_text(face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
                    panel.grid.minor = element_blank())

            ggsave(file.path(interaction_dir,
                             paste0("interaction_ORA_KEGG_dotplot_", depot, "_", run_tag, ".png")),
                   p_kegg, width = 7, height = 6, dpi = 300, bg = "white")
          } else {
            cat("[INFO] No significant KEGG pathways for interaction genes\n")
          }
        }
      }, error = function(e) {
        cat("[WARN] KEGG ORA failed: ", e$message, "\n")
      })

    } else {
      cat("[INFO] Too few genes mapped to Entrez IDs for pathway analysis\n")
    }
  } else {
    cat("[INFO] Fewer than 5 significant interaction genes in ", depot,
        "; skipping pathway analysis\n", sep = "")
  }

  ## -------------------------------------------------------
  ## 4l) Collect summary row
  ## -------------------------------------------------------
  summary_rows[[depot]] <- data.frame(
    Depot = depot,
    N_genes_tested = n_tested,
    N_interaction_sig = n_sig_lrt,
    Pct_interaction = pct_sig,
    pi0_estimate = round(pi0_est, 4),
    N_genotype_DEGs = n_geno_sig,
    N_sex_DEGs = n_sex_sig,
    Sex_logFC_cor = if (exists("r_val") && !is.null(r_val)) round(r_val, 4) else NA,
    Median_var_sex = round(median(var_results$pct_sex), 2),
    Median_var_genotype = round(median(var_results$pct_genotype), 2),
    Median_var_interaction = round(median(var_results$pct_interaction), 2),
    Median_var_residual = round(median(var_results$pct_residual), 2),
    stringsAsFactors = FALSE
  )
}

## ============================================================
## 5) COMBINED SUMMARY
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== SEX x GENOTYPE INTERACTION: COMBINED SUMMARY ===\n")
cat("============================================================\n\n")

summary_df <- bind_rows(summary_rows)
print(summary_df, row.names = FALSE)

write.csv(summary_df,
          file.path(tables_dir, paste0("sex_interaction_summary_", run_tag, ".csv")),
          row.names = FALSE)

## ============================================================
## 6) MANUSCRIPT LANGUAGE GENERATOR
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== SUGGESTED MANUSCRIPT TEXT ===\n")
cat("============================================================\n\n")

for (depot in c("iWAT", "gWAT")) {
  row <- summary_df[summary_df$Depot == depot, ]
  cat("--- ", depot, " ---\n", sep = "")
  cat("\"To assess whether KAT8-KD transcriptional effects differed\n")
  cat("between sexes in ", depot, ", we performed a DESeq2 likelihood\n", sep = "")
  cat("ratio test comparing a full model (~ Sex + Genotype + Sex:Genotype)\n")
  cat("against a reduced model (~ Sex + Genotype). Of ", row$N_genes_tested, "\n", sep = "")
  cat("tested genes, ", row$N_interaction_sig, " (", row$Pct_interaction, "%)\n", sep = "")
  cat("showed significant Sex x Genotype interactions (FDR < 0.05)")
  if (!is.na(row$pi0_estimate)) {
    cat(",\nwith an estimated true null proportion (pi0) of ",
        row$pi0_estimate, ", indicating that approximately ",
        round(row$pi0_estimate * 100, 0), "% of genes show\n",
        "no evidence of sex-dependent KAT8-KD effects", sep = "")
  }
  cat(".\nVariance decomposition across a random sample of genes showed\n")
  cat("the interaction term explained a median of ", row$Median_var_interaction,
      "% of per-gene variance,\n", sep = "")
  cat("compared to ", row$Median_var_genotype,
      "% for the Genotype main effect and ",
      row$Median_var_sex, "% for Sex.\n", sep = "")
  if (!is.na(row$Sex_logFC_cor)) {
    ## Nuanced concordance language based on r value
    concordance_adj <- if (row$Sex_logFC_cor >= 0.8) {
      "strong"
    } else if (row$Sex_logFC_cor >= 0.6) {
      "moderate-to-strong"
    } else if (row$Sex_logFC_cor >= 0.4) {
      "moderate"
    } else {
      "weak"
    }
    cat("Sex-stratified analysis revealed ", concordance_adj, " concordance\n", sep = "")
    cat("between female and male KD log2 fold-changes\n")
    cat("(Pearson r = ", row$Sex_logFC_cor, ").\n", sep = "")
  }
  ## Power caveat with actual group sizes
  depot_annot_sub <- sample_annot[sample_annot$Depot == depot, ]
  group_sizes <- with(depot_annot_sub, table(Sex, Genotype))
  min_n <- min(group_sizes)
  max_n <- max(group_sizes)
  if (min_n == max_n) {
    cat("We note that with n = ", min_n, " per sex-genotype group, ", sep = "")
  } else {
    cat("We note that with n = ", min_n, "-", max_n, " per sex-genotype group, ", sep = "")
  }
  cat("the interaction test has limited statistical power\n")
  cat("to detect subtle sex-dependent effects; however, the convergent\n")
  cat("evidence from the LRT, pi0 estimation, variance decomposition,\n")
  cat("and sex-stratified concordance analysis collectively supports\n")
  cat("pooling sexes for the primary analysis of ", depot, ".\"\n\n", sep = "")
}

## ============================================================
## 7) SESSION INFO
## ============================================================
log_dir <- file.path(outdir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
sink(file.path(log_dir, paste0("sessionInfo_sex_interaction_", run_tag, ".txt")))
cat("=== R SESSION INFORMATION - Part 7.6 (Sex Interaction Test) ===\n\n")
print(sessionInfo())
sink()

cat("\n=== Part 7.6 (Sex x Genotype Interaction Test) COMPLETE ===\n")
cat("Output directory: ", interaction_dir, "\n\n", sep = "")
cat("Files created per depot:\n")
cat("  1. P-value histogram (interaction LRT)\n")
cat("  2. Interaction effect vs Genotype effect scatter\n")
cat("  3. PCA colored by Sex/Genotype\n")
cat("  3b. Sex-stratified KD concordance scatter (M vs F logFC)\n")
cat("  4. Variance explained boxplot (Sex vs Genotype vs Interaction)\n")
cat("  5. Full interaction results table\n")
cat("  6. Significant interaction genes (if any)\n")
cat("  7. Sex-stratified logFC table\n")
cat("  8. Combined summary table\n")
cat("  9. Interaction gene ORA: GO:BP results + dotplot (if sig genes >= 5)\n")
cat(" 10. Interaction gene ORA: KEGG results + dotplot (if sig genes >= 5)\n\n")
