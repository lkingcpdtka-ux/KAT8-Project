#!/usr/bin/env Rscript

################################################################################
## KAT8 RNA-seq Diagnostic: Combined vs. Split Depot Analysis
##
## This script compares two analysis strategies:
##   1. COMBINED: All depots analyzed together (current approach)
##   2. SPLIT: iWAT and gWAT analyzed separately
##
## Purpose: Determine which approach is more appropriate for KAT8 chromatin
##          modifier study with depot-specific effects
##
## Output:
##   - DEG count comparisons
##   - Dispersion plot comparisons
##   - Overlap analysis
##   - Recommendation
################################################################################

## Packages
required_pkgs <- c("dplyr", "ggplot2", "gridExtra")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
})

cat("\n")
cat("================================================================================\n")
cat("=== DIAGNOSTIC: COMBINED VS. SPLIT DEPOT ANALYSIS ===\n")
cat("================================================================================\n")
cat("\n")

## Setup
outdir <- file.path(getwd(), "diagnostic_output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Load data
cat("Loading data...\n")
counts_raw <- read.delim("raw_counts.txt", row.names = 1, stringsAsFactors = FALSE)
metadata <- read.csv("metadata.csv", stringsAsFactors = FALSE)

## Remove gene_id column if present
if ("gene_id" %in% colnames(counts_raw)) {
  counts_raw <- counts_raw[, -which(colnames(counts_raw) == "gene_id"), drop = FALSE]
}

## Match samples
shared_samples <- intersect(colnames(counts_raw), metadata$Sample)
counts_matrix <- as.matrix(counts_raw[, shared_samples])
metadata <- metadata[match(shared_samples, metadata$Sample), ]

## Remove outliers (from part1)
outlier_samples <- c("JS_08", "JS_28")
keep_idx <- !metadata$Sample %in% outlier_samples
counts_matrix <- counts_matrix[, keep_idx, drop = FALSE]
metadata <- metadata[keep_idx, ]

cat("Samples after outlier removal: ", ncol(counts_matrix), "\n")
cat("iWAT samples: ", sum(metadata$Tissue == "iWAT"), "\n")
cat("gWAT samples: ", sum(metadata$Tissue == "gWAT"), "\n\n")

## Filtering parameters (global filtering as updated)
min_count <- 10
min_samples <- 4

################################################################################
## APPROACH 1: COMBINED ANALYSIS (Current)
################################################################################

cat("================================================================================\n")
cat("APPROACH 1: COMBINED ANALYSIS (All depots together)\n")
cat("================================================================================\n\n")

## Create GroupFull factor
metadata$GroupFull <- paste(metadata$Tissue, metadata$Sex, metadata$Genotype, sep = "_")
metadata$GroupFull <- factor(metadata$GroupFull)

## Filter genes
keep_combined <- rowSums(counts_matrix >= min_count) >= min_samples
counts_filt_combined <- counts_matrix[keep_combined, , drop = FALSE]

cat("Genes after filtering: ", nrow(counts_filt_combined), " (",
    round(100 * nrow(counts_filt_combined) / nrow(counts_matrix), 1), "%)\n\n")

## DESeq2 combined
dds_combined <- DESeqDataSetFromMatrix(
  countData = counts_filt_combined,
  colData = metadata,
  design = ~ 0 + GroupFull
)

cat("Running DESeq2 (combined)...\n")
dds_combined <- DESeq(dds_combined, quiet = TRUE)

## Extract results for each contrast
contrasts_combined <- list(
  iWAT_F = c("GroupFull", "iWAT_F_KAT8KD", "iWAT_F_CTL"),
  iWAT_M = c("GroupFull", "iWAT_M_KAT8KD", "iWAT_M_CTL"),
  gWAT_F = c("GroupFull", "gWAT_F_KAT8KD", "gWAT_F_CTL"),
  gWAT_M = c("GroupFull", "gWAT_M_KAT8KD", "gWAT_M_CTL")
)

results_combined <- list()
for (cn in names(contrasts_combined)) {
  res <- results(dds_combined, contrast = contrasts_combined[[cn]])
  res_df <- as.data.frame(res) %>%
    dplyr::filter(!is.na(padj))
  results_combined[[cn]] <- res_df
}

## Count DEGs
deg_counts_combined <- data.frame(
  Contrast = names(results_combined),
  Approach = "Combined",
  N_Tested = sapply(results_combined, nrow),
  DEG_FDR05 = sapply(results_combined, function(x) sum(x$padj < 0.05, na.rm = TRUE)),
  DEG_FDR05_FC1 = sapply(results_combined, function(x)
    sum(x$padj < 0.05 & abs(x$log2FoldChange) > 1, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

print(deg_counts_combined)
cat("\n")

################################################################################
## APPROACH 2: SPLIT ANALYSIS (Separate by depot)
################################################################################

cat("================================================================================\n")
cat("APPROACH 2: SPLIT ANALYSIS (Separate iWAT and gWAT)\n")
cat("================================================================================\n\n")

results_split <- list()
deg_counts_split_list <- list()

for (depot in c("iWAT", "gWAT")) {
  cat("--- Analyzing ", depot, " only ---\n", sep = "")

  ## Subset data
  depot_idx <- metadata$Tissue == depot
  counts_depot <- counts_matrix[, depot_idx, drop = FALSE]
  metadata_depot <- metadata[depot_idx, ]

  cat("Samples: ", ncol(counts_depot), "\n")

  ## Create Sex_Genotype factor
  metadata_depot$Sex_Genotype <- paste(metadata_depot$Sex, metadata_depot$Genotype, sep = "_")
  metadata_depot$Sex_Genotype <- factor(metadata_depot$Sex_Genotype)

  ## Filter genes (depot-specific)
  keep_depot <- rowSums(counts_depot >= min_count) >= min_samples
  counts_filt_depot <- counts_depot[keep_depot, , drop = FALSE]

  cat("Genes after filtering: ", nrow(counts_filt_depot), " (",
      round(100 * nrow(counts_filt_depot) / nrow(counts_depot), 1), "%)\n")

  ## DESeq2 split
  dds_depot <- DESeqDataSetFromMatrix(
    countData = counts_filt_depot,
    colData = metadata_depot,
    design = ~ 0 + Sex_Genotype
  )

  cat("Running DESeq2 (", depot, ")...\n", sep = "")
  dds_depot <- DESeq(dds_depot, quiet = TRUE)

  ## Extract results for each sex
  contrasts_depot <- list(
    F = c("Sex_Genotype", "F_KAT8KD", "F_CTL"),
    M = c("Sex_Genotype", "M_KAT8KD", "M_CTL")
  )

  for (sex in names(contrasts_depot)) {
    cn <- paste0(depot, "_", sex)
    res <- results(dds_depot, contrast = contrasts_depot[[sex]])
    res_df <- as.data.frame(res) %>%
      dplyr::filter(!is.na(padj))
    results_split[[cn]] <- res_df
  }

  cat("\n")
}

## Count DEGs
deg_counts_split <- data.frame(
  Contrast = names(results_split),
  Approach = "Split",
  N_Tested = sapply(results_split, nrow),
  DEG_FDR05 = sapply(results_split, function(x) sum(x$padj < 0.05, na.rm = TRUE)),
  DEG_FDR05_FC1 = sapply(results_split, function(x)
    sum(x$padj < 0.05 & abs(x$log2FoldChange) > 1, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

print(deg_counts_split)
cat("\n")

################################################################################
## COMPARISON
################################################################################

cat("================================================================================\n")
cat("COMPARISON: COMBINED VS. SPLIT\n")
cat("================================================================================\n\n")

## Combine results
deg_comparison <- rbind(deg_counts_combined, deg_counts_split)

## Calculate differences
cat("--- DEG Count Differences (Split - Combined) ---\n")
for (cn in names(results_combined)) {
  combined_count <- deg_counts_combined$DEG_FDR05_FC1[deg_counts_combined$Contrast == cn]
  split_count <- deg_counts_split$DEG_FDR05_FC1[deg_counts_split$Contrast == cn]
  diff <- split_count - combined_count
  pct_diff <- round(100 * diff / combined_count, 1)

  cat(sprintf("%-10s: %+5d DEGs (%+5.1f%%)\n", cn, diff, pct_diff))
}
cat("\n")

## Gene overlap analysis
cat("--- Gene Overlap Analysis (FDR<0.05, |FC|>1) ---\n")
for (cn in names(results_combined)) {
  genes_combined <- rownames(results_combined[[cn]])[
    results_combined[[cn]]$padj < 0.05 & abs(results_combined[[cn]]$log2FoldChange) > 1
  ]
  genes_split <- rownames(results_split[[cn]])[
    results_split[[cn]]$padj < 0.05 & abs(results_split[[cn]]$log2FoldChange) > 1
  ]

  n_overlap <- length(intersect(genes_combined, genes_split))
  n_combined_only <- length(setdiff(genes_combined, genes_split))
  n_split_only <- length(setdiff(genes_split, genes_combined))

  pct_overlap <- if (length(genes_combined) > 0) {
    round(100 * n_overlap / length(genes_combined), 1)
  } else 0

  cat(sprintf("%-10s: %4d overlap (%5.1f%%), %4d combined-only, %4d split-only\n",
              cn, n_overlap, pct_overlap, n_combined_only, n_split_only))
}
cat("\n")

## Save comparison table
write.csv(deg_comparison, file = file.path(outdir, "deg_comparison_combined_vs_split.csv"),
          row.names = FALSE)
cat("Saved comparison table to: ", file.path(outdir, "deg_comparison_combined_vs_split.csv"), "\n\n")

################################################################################
## VISUALIZATION
################################################################################

cat("Creating comparison plots...\n")

## Plot 1: DEG counts comparison
p1 <- ggplot(deg_comparison, aes(x = Contrast, y = DEG_FDR05_FC1, fill = Approach)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Combined" = "#E69F00", "Split" = "#56B4E9")) +
  labs(
    title = "DEG Counts: Combined vs. Split Analysis",
    subtitle = "FDR < 0.05, |log2FC| > 1",
    x = "Contrast",
    y = "Number of DEGs"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(file.path(outdir, "deg_counts_comparison.png"), plot = p1,
       width = 8, height = 6, dpi = 300)

## Plot 2: Genes tested comparison
p2 <- ggplot(deg_comparison, aes(x = Contrast, y = N_Tested, fill = Approach)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("Combined" = "#E69F00", "Split" = "#56B4E9")) +
  labs(
    title = "Genes Tested: Combined vs. Split Analysis",
    x = "Contrast",
    y = "Number of Genes Tested"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(file.path(outdir, "genes_tested_comparison.png"), plot = p2,
       width = 8, height = 6, dpi = 300)

cat("Saved plots to: ", outdir, "\n\n")

################################################################################
## RECOMMENDATION
################################################################################

cat("================================================================================\n")
cat("RECOMMENDATION\n")
cat("================================================================================\n\n")

## Calculate average differences
avg_diff_iwat <- mean(c(
  deg_counts_split$DEG_FDR05_FC1[deg_counts_split$Contrast == "iWAT_F"] -
    deg_counts_combined$DEG_FDR05_FC1[deg_counts_combined$Contrast == "iWAT_F"],
  deg_counts_split$DEG_FDR05_FC1[deg_counts_split$Contrast == "iWAT_M"] -
    deg_counts_combined$DEG_FDR05_FC1[deg_counts_combined$Contrast == "iWAT_M"]
))

avg_diff_gwat <- mean(c(
  deg_counts_split$DEG_FDR05_FC1[deg_counts_split$Contrast == "gWAT_F"] -
    deg_counts_combined$DEG_FDR05_FC1[deg_counts_combined$Contrast == "gWAT_F"],
  deg_counts_split$DEG_FDR05_FC1[deg_counts_split$Contrast == "gWAT_M"] -
    deg_counts_combined$DEG_FDR05_FC1[deg_counts_combined$Contrast == "gWAT_M"]
))

cat("Average DEG difference (Split - Combined):\n")
cat("  iWAT: ", round(avg_diff_iwat, 0), " DEGs\n", sep = "")
cat("  gWAT: ", round(avg_diff_gwat, 0), " DEGs\n\n", sep = "")

if (abs(avg_diff_iwat) < 100 && abs(avg_diff_gwat) < 100) {
  cat("** RECOMMENDATION: KEEP COMBINED ANALYSIS **\n\n")
  cat("Reasons:\n")
  cat("  1. DEG counts are similar between approaches (< 100 difference)\n")
  cat("  2. Combined approach has more statistical power (more samples)\n")
  cat("  3. Allows formal testing of depot x treatment interaction\n")
  cat("  4. Standard approach for depot comparison papers\n")
  cat("  5. Reviewers expect this for comparative studies\n\n")
  cat("Note: If you want to emphasize depot-specific biology, you can:\n")
  cat("  - Use combined analysis for main results\n")
  cat("  - Show split analysis in supplement as sensitivity analysis\n")
} else if (avg_diff_gwat > 100) {
  cat("** RECOMMENDATION: CONSIDER SPLIT ANALYSIS **\n\n")
  cat("Reasons:\n")
  cat("  1. Split analysis finds substantially more DEGs in gWAT (+", round(avg_diff_gwat, 0), ")\n", sep = "")
  cat("  2. Suggests depot-specific normalization improves gWAT sensitivity\n")
  cat("  3. Combined analysis may mask gWAT-specific effects\n")
  cat("  4. Appropriate for chromatin modifier with depot-specific action\n\n")
  cat("Suggested approach:\n")
  cat("  - Run split analysis for main results\n")
  cat("  - Analyze iWAT and gWAT as independent biological questions\n")
  cat("  - Compare qualitatively in discussion (not statistically)\n")
} else {
  cat("** RECOMMENDATION: EITHER APPROACH ACCEPTABLE **\n\n")
  cat("Reasons:\n")
  cat("  1. Results are similar between approaches\n")
  cat("  2. Choice depends on your biological question:\n")
  cat("     - 'Is KAT8 effect depot-dependent?' → Combined\n")
  cat("     - 'What does KAT8 do in each depot?' → Split\n\n")
  cat("Suggested approach:\n")
  cat("  - Use combined for main figures (standard)\n")
  cat("  - Use split for supplement (shows robustness)\n")
}

cat("\n")
cat("Diagnostic complete!\n")
cat("Output saved to: ", outdir, "\n")
cat("\n")
