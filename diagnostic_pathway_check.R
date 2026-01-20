#!/usr/bin/env Rscript

## =========================================================
## DIAGNOSTIC SCRIPT: Why few down-regulated pathways in GWAT?
## =========================================================
## This script checks:
## 1. DEG counts (up vs down) for all contrasts
## 2. Wald statistic distributions (to check for bias)
## 3. Gene overlap between iWAT and GWAT
## 4. Pathway enrichment parameters
## =========================================================

library(dplyr)
library(ggplot2)

## Find most recent Part 1 run
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))
tables_dir <- file.path(outdir, "tables")

cat("\n==========================================================\n")
cat("=== DIAGNOSTIC: GWAT DOWN-REGULATED PATHWAY ANALYSIS ===\n")
cat("==========================================================\n")
cat("Run directory: ", outdir, "\n", sep = "")
cat("Run tag: ", run_tag, "\n\n", sep = "")

## Parameters (match your analysis)
logFC_cut <- 1
fdr_cut <- 0.05

## 1) Load all DE tables and count DEGs
cat("=== STEP 1: DEG COUNTS (UP vs DOWN) ===\n\n")

de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("OVERALL", de_files)]  # Skip overall contrast

deg_summary <- data.frame()

for (de_file in de_files) {
  de_filename <- basename(de_file)
  contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

  de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Count DEGs
  n_up <- sum(de_table$adj.P.Val < fdr_cut & de_table$logFC > logFC_cut, na.rm = TRUE)
  n_down <- sum(de_table$adj.P.Val < fdr_cut & de_table$logFC < -logFC_cut, na.rm = TRUE)
  n_total <- n_up + n_down

  ## Count genes with stat column (for fgsea)
  n_valid_stat <- sum(is.finite(de_table$stat), na.rm = TRUE)
  n_pos_stat <- sum(de_table$stat > 0, na.rm = TRUE)
  n_neg_stat <- sum(de_table$stat < 0, na.rm = TRUE)

  deg_summary <- rbind(deg_summary, data.frame(
    Contrast = contrast_name,
    Up_DEGs = n_up,
    Down_DEGs = n_down,
    Total_DEGs = n_total,
    Pct_Up = round(100 * n_up / n_total, 1),
    Pct_Down = round(100 * n_down / n_total, 1),
    Valid_Stats = n_valid_stat,
    Pos_Stats = n_pos_stat,
    Neg_Stats = n_neg_stat,
    Pct_Pos_Stats = round(100 * n_pos_stat / n_valid_stat, 1),
    Pct_Neg_Stats = round(100 * n_neg_stat / n_valid_stat, 1),
    stringsAsFactors = FALSE
  ))
}

print(deg_summary)

cat("\n--- KEY OBSERVATIONS ---\n")
cat("Look for:\n")
cat("  - GWAT contrasts with <50 down-regulated DEGs (low power for pathways)\n")
cat("  - Asymmetry in Pct_Down vs Pct_Up (indicates directional bias)\n")
cat("  - Low Pct_Neg_Stats (means few genes with negative Wald stats)\n\n")

## 2) Check Wald statistic distributions
cat("\n=== STEP 2: WALD STATISTIC DISTRIBUTIONS ===\n\n")

for (de_file in de_files) {
  de_filename <- basename(de_file)
  contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

  de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Get stats
  stats <- de_table$stat[is.finite(de_table$stat)]

  ## Summary statistics
  cat("Contrast: ", contrast_name, "\n", sep = "")
  cat("  Wald stat range: [", round(min(stats), 2), ", ", round(max(stats), 2), "]\n", sep = "")
  cat("  Wald stat median: ", round(median(stats), 2), "\n", sep = "")
  cat("  Wald stat mean: ", round(mean(stats), 2), " (should be ~0 if balanced)\n", sep = "")

  if (abs(mean(stats)) > 0.5) {
    cat("  *** WARNING: Mean Wald stat is far from 0 - indicates directional bias! ***\n")
  }

  ## Count extreme negative stats (these drive down-regulated pathways in fgsea)
  n_very_neg <- sum(stats < -2, na.rm = TRUE)
  cat("  Genes with Wald stat < -2: ", n_very_neg, " (these drive down-regulated pathway enrichment)\n", sep = "")
  cat("\n")
}

## 3) Gene overlap analysis
cat("\n=== STEP 3: GENE OVERLAP (iWAT vs GWAT) ===\n\n")

## Load iWAT F down genes
iwat_f_file <- de_files[grepl("iWAT_F_KD_vs_CTL", de_files)]
iwat_f_de <- read.csv(iwat_f_file, row.names = 1, stringsAsFactors = FALSE)
iwat_f_down <- rownames(iwat_f_de)[iwat_f_de$adj.P.Val < fdr_cut & iwat_f_de$logFC < -logFC_cut]

## Load gWAT F down genes
gwat_f_file <- de_files[grepl("gWAT_F_KD_vs_CTL", de_files)]
gwat_f_de <- read.csv(gwat_f_file, row.names = 1, stringsAsFactors = FALSE)
gwat_f_down <- rownames(gwat_f_de)[gwat_f_de$adj.P.Val < fdr_cut & gwat_f_de$logFC < -logFC_cut]

## Overlap
overlap_f <- intersect(iwat_f_down, gwat_f_down)
cat("iWAT_F down-regulated DEGs: ", length(iwat_f_down), "\n", sep = "")
cat("gWAT_F down-regulated DEGs: ", length(gwat_f_down), "\n", sep = "")
cat("Overlap: ", length(overlap_f), " genes (",
    round(100 * length(overlap_f) / length(gwat_f_down), 1), "% of gWAT_F down DEGs)\n\n", sep = "")

if (length(overlap_f) > 0) {
  cat("Shared down-regulated genes (first 20):\n")
  print(head(overlap_f, 20))
  cat("\n")
}

## Same for males
iwat_m_file <- de_files[grepl("iWAT_M_KD_vs_CTL", de_files)]
iwat_m_de <- read.csv(iwat_m_file, row.names = 1, stringsAsFactors = FALSE)
iwat_m_down <- rownames(iwat_m_de)[iwat_m_de$adj.P.Val < fdr_cut & iwat_m_de$logFC < -logFC_cut]

gwat_m_file <- de_files[grepl("gWAT_M_KD_vs_CTL", de_files)]
gwat_m_de <- read.csv(gwat_m_file, row.names = 1, stringsAsFactors = FALSE)
gwat_m_down <- rownames(gwat_m_de)[gwat_m_de$adj.P.Val < fdr_cut & gwat_m_de$logFC < -logFC_cut]

overlap_m <- intersect(iwat_m_down, gwat_m_down)
cat("iWAT_M down-regulated DEGs: ", length(iwat_m_down), "\n", sep = "")
cat("gWAT_M down-regulated DEGs: ", length(gwat_m_down), "\n", sep = "")
cat("Overlap: ", length(overlap_m), " genes (",
    round(100 * length(overlap_m) / length(gwat_m_down), 1), "% of gWAT_M down DEGs)\n\n", sep = "")

if (length(overlap_m) > 0) {
  cat("Shared down-regulated genes (first 20):\n")
  print(head(overlap_m, 20))
  cat("\n")
}

## 4) Check pathway enrichment results
cat("\n=== STEP 4: PATHWAY ENRICHMENT SUMMARY ===\n\n")

## Check ORA results
ora_files <- list.files(tables_dir, pattern = "^ORA_.*gWAT.*\\.csv$", full.names = TRUE)
if (length(ora_files) > 0) {
  cat("--- ORA Results for GWAT ---\n")
  for (ora_file in ora_files) {
    ora_basename <- basename(ora_file)
    ora_data <- read.csv(ora_file, stringsAsFactors = FALSE)
    cat(ora_basename, ": ", nrow(ora_data), " significant pathways\n", sep = "")
  }
  cat("\n")
} else {
  cat("No ORA results found for GWAT\n\n")
}

## Check fgsea results
fgsea_files <- list.files(tables_dir, pattern = "^fgsea_.*gWAT.*\\.csv$", full.names = TRUE)
if (length(fgsea_files) > 0) {
  cat("--- fgsea Results for GWAT ---\n")
  for (fgsea_file in fgsea_files) {
    fgsea_basename <- basename(fgsea_file)
    fgsea_data <- read.csv(fgsea_file, stringsAsFactors = FALSE)
    n_up <- sum(fgsea_data$NES > 0, na.rm = TRUE)
    n_down <- sum(fgsea_data$NES < 0, na.rm = TRUE)
    cat(fgsea_basename, ": ", nrow(fgsea_data), " total (",
        n_up, " up, ", n_down, " down)\n", sep = "")
  }
  cat("\n")
} else {
  cat("No fgsea results found for GWAT\n\n")
}

## 5) Recommendations
cat("\n==========================================================\n")
cat("=== RECOMMENDATIONS ===\n")
cat("==========================================================\n\n")

## Check if GWAT has <50 down DEGs
gwat_f_down_count <- deg_summary$Down_DEGs[deg_summary$Contrast == "gWAT_F_KD_vs_CTL"]
gwat_m_down_count <- deg_summary$Down_DEGs[deg_summary$Contrast == "gWAT_M_KD_vs_CTL"]

if (gwat_f_down_count < 50 || gwat_m_down_count < 50) {
  cat("⚠️  ISSUE: GWAT has very few down-regulated DEGs\n")
  cat("    gWAT_F: ", gwat_f_down_count, " down DEGs\n", sep = "")
  cat("    gWAT_M: ", gwat_m_down_count, " down DEGs\n\n", sep = "")
  cat("This is likely WHY you're not seeing many down-regulated pathways.\n\n")
  cat("POSSIBLE REASONS:\n")
  cat("  1. Biology: KAT8 knockdown has less effect on gene repression in GWAT vs iWAT\n")
  cat("  2. Technical: Lower statistical power in GWAT (more variability?)\n")
  cat("  3. Filtering: Stringent QC cutoffs removed more GWAT genes\n\n")
  cat("OPTIONS TO TRY:\n")
  cat("  A. Accept it - this may be real biology (GWAT ≠ iWAT)\n")
  cat("  B. Relax thresholds for GWAT only:\n")
  cat("     - Use FDR < 0.1 (instead of 0.05)\n")
  cat("     - Use |logFC| > 0.5 (instead of 1.0)\n")
  cat("  C. Use GSEA/fgsea ONLY for GWAT (doesn't require hard cutoffs)\n")
  cat("  D. Check if there are batch effects specific to GWAT samples\n\n")
} else {
  cat("✅ GWAT has sufficient down-regulated DEGs (", gwat_f_down_count, " F, ",
      gwat_m_down_count, " M)\n\n", sep = "")
  cat("The issue may be:\n")
  cat("  1. Down-regulated genes don't cluster into known pathways\n")
  cat("  2. Pathways are below significance threshold\n")
  cat("  3. Need to loosen pathway enrichment parameters\n\n")
  cat("OPTIONS TO TRY:\n")
  cat("  A. Lower pathway FDR cutoff (FDR < 0.1 instead of 0.05)\n")
  cat("  B. Use less stringent gene set sizes (minGSSize=3 instead of 5)\n")
  cat("  C. Try different pathway databases (Reactome, GO:MF, GO:CC)\n\n")
}

cat("==========================================================\n\n")

## Save summary
summary_file <- paste0("DIAGNOSTIC_pathway_analysis_", run_tag, ".txt")
sink(file.path(outdir, "logs", summary_file))
print(deg_summary)
sink()
cat("Saved summary to: ", file.path(outdir, "logs", summary_file), "\n\n")

cat("=== DIAGNOSTIC COMPLETE ===\n\n")
