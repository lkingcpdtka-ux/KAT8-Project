#!/usr/bin/env Rscript

################################################################################
## Compare Prefilter Strategies: Per-Group vs Overall
## Purpose: Show impact of different prefiltering approaches on gene retention
################################################################################

library(tidyverse)

cat("=== COMPARING PREFILTER STRATEGIES ===\n\n")

## Load data
count_matrix_cells <- read.csv("data/KAT8_bulkRNAseq_cells_counts.csv", row.names = 1)
sample_annot <- data.frame(
  Sample = colnames(count_matrix_cells),
  Genotype = rep(c("CTL", "KAT8KD"), each = 4)
)

cat("Total genes before filtering:", nrow(count_matrix_cells), "\n\n")

## Parameters
min_count <- 10
min_samples_per_group <- 2

## Get group indices
ctl_samples <- sample_annot$Sample[sample_annot$Genotype == "CTL"]
kat8kd_samples <- sample_annot$Sample[sample_annot$Genotype == "KAT8KD"]

################################################################################
## STRATEGY 1: PER-GROUP WITH OR (Current approach)
################################################################################
cat("--- STRATEGY 1: PER-GROUP WITH OR (current) ---\n")
cat("Logic: Keep if ≥10 counts in ≥2 samples in CTL OR KAT8KD\n")
cat("Rationale: Captures genes expressed in EITHER group\n\n")

keep_ctl <- rowSums(count_matrix_cells[, ctl_samples] >= min_count) >= min_samples_per_group
keep_kat8kd <- rowSums(count_matrix_cells[, kat8kd_samples] >= min_count) >= min_samples_per_group
keep_per_group_OR <- keep_ctl | keep_kat8kd

cat("Genes passing filter:", sum(keep_per_group_OR), "\n")
cat("Percentage retained:", round(100 * mean(keep_per_group_OR), 2), "%\n\n")

# Breakdown
cat("  - Expressed in CTL only:", sum(keep_ctl & !keep_kat8kd), "\n")
cat("  - Expressed in KAT8KD only:", sum(!keep_ctl & keep_kat8kd), "\n")
cat("  - Expressed in BOTH groups:", sum(keep_ctl & keep_kat8kd), "\n\n")

################################################################################
## STRATEGY 2: OVERALL (≥2 out of all 8 samples)
################################################################################
cat("--- STRATEGY 2: OVERALL (≥2 out of 8 total samples) ---\n")
cat("Logic: Keep if ≥10 counts in ≥2 samples across ALL 8 samples\n")
cat("Rationale: Simple threshold ignoring group structure\n\n")

keep_overall_2of8 <- rowSums(count_matrix_cells >= min_count) >= 2

cat("Genes passing filter:", sum(keep_overall_2of8), "\n")
cat("Percentage retained:", round(100 * mean(keep_overall_2of8), 2), "%\n\n")

################################################################################
## STRATEGY 3: OVERALL (≥4 out of all 8 samples - stricter)
################################################################################
cat("--- STRATEGY 3: OVERALL (≥4 out of 8 total samples) ---\n")
cat("Logic: Keep if ≥10 counts in ≥4 samples across ALL 8 samples\n")
cat("Rationale: Requires expression in half of all samples\n\n")

keep_overall_4of8 <- rowSums(count_matrix_cells >= min_count) >= 4

cat("Genes passing filter:", sum(keep_overall_4of8), "\n")
cat("Percentage retained:", round(100 * mean(keep_overall_4of8), 2), "%\n\n")

################################################################################
## STRATEGY 4: PER-GROUP WITH AND (most conservative)
################################################################################
cat("--- STRATEGY 4: PER-GROUP WITH AND (most conservative) ---\n")
cat("Logic: Keep if ≥10 counts in ≥2 samples in CTL AND KAT8KD\n")
cat("Rationale: Only keeps genes expressed in BOTH groups\n\n")

keep_per_group_AND <- keep_ctl & keep_kat8kd

cat("Genes passing filter:", sum(keep_per_group_AND), "\n")
cat("Percentage retained:", round(100 * mean(keep_per_group_AND), 2), "%\n\n")

################################################################################
## COMPARISON: What genes are we gaining/losing?
################################################################################
cat("=== COMPARISON SUMMARY ===\n\n")

cat("Genes ONLY in per-group OR (lost with other methods):\n")
lost_overall_2of8 <- sum(keep_per_group_OR & !keep_overall_2of8)
lost_overall_4of8 <- sum(keep_per_group_OR & !keep_overall_4of8)
lost_per_group_AND <- sum(keep_per_group_OR & !keep_per_group_AND)

cat("  vs Overall (≥2 of 8):", lost_overall_2of8, "genes\n")
cat("  vs Overall (≥4 of 8):", lost_overall_4of8, "genes\n")
cat("  vs Per-group AND:", lost_per_group_AND, "genes\n\n")

cat("Genes GAINED with per-group OR (vs other methods):\n")
gained_overall_2of8 <- sum(!keep_per_group_OR & keep_overall_2of8)
gained_overall_4of8 <- sum(!keep_per_group_OR & keep_overall_4of8)
gained_per_group_AND <- sum(!keep_per_group_OR & keep_per_group_AND)

cat("  vs Overall (≥2 of 8):", gained_overall_2of8, "genes\n")
cat("  vs Overall (≥4 of 8):", gained_overall_4of8, "genes\n")
cat("  vs Per-group AND:", gained_per_group_AND, "genes\n\n")

################################################################################
## CRITICAL: Show examples of genes lost with overall/AND approaches
################################################################################
cat("=== EXAMPLE GENES THAT WOULD BE LOST ===\n\n")
cat("These are genes expressed in ONE group but not the other\n")
cat("(These are exactly the genes we want for DE analysis!)\n\n")

# Find genes expressed in CTL but not KAT8KD
ctl_specific <- keep_ctl & !keep_kat8kd
if (sum(ctl_specific) > 0) {
  cat("Genes expressed in CTL only:", sum(ctl_specific), "\n")
  cat("Example genes:\n")
  examples_ctl <- head(rownames(count_matrix_cells)[ctl_specific], 5)
  for (gene in examples_ctl) {
    counts_ctl <- count_matrix_cells[gene, ctl_samples]
    counts_kat8kd <- count_matrix_cells[gene, kat8kd_samples]
    cat(sprintf("  %s: CTL=[%s] KAT8KD=[%s]\n",
                gene,
                paste(round(counts_ctl, 0), collapse=","),
                paste(round(counts_kat8kd, 0), collapse=",")))
  }
  cat("\n")
}

# Find genes expressed in KAT8KD but not CTL
kat8kd_specific <- !keep_ctl & keep_kat8kd
if (sum(kat8kd_specific) > 0) {
  cat("Genes expressed in KAT8KD only:", sum(kat8kd_specific), "\n")
  cat("Example genes:\n")
  examples_kat8kd <- head(rownames(count_matrix_cells)[kat8kd_specific], 5)
  for (gene in examples_kat8kd) {
    counts_ctl <- count_matrix_cells[gene, ctl_samples]
    counts_kat8kd <- count_matrix_cells[gene, kat8kd_samples]
    cat(sprintf("  %s: CTL=[%s] KAT8KD=[%s]\n",
                gene,
                paste(round(counts_ctl, 0), collapse=","),
                paste(round(counts_kat8kd, 0), collapse=",")))
  }
  cat("\n")
}

################################################################################
## RECOMMENDATION
################################################################################
cat("=== RECOMMENDATION ===\n\n")

cat("For TWO-GROUP differential expression analysis:\n")
cat("✓ BEST: Per-group with OR (Strategy 1 - current approach)\n\n")

cat("Why?\n")
cat("1. Captures genes that are ON in one group and OFF in the other\n")
cat("2. These group-specific genes are EXACTLY what DE analysis looks for\n")
cat("3. Standard practice in DESeq2 vignette and Love et al. 2014\n")
cat("4. Filtering out these genes would MISS real biological signal\n\n")

cat("When would you use 'overall' filtering?\n")
cat("- Multi-group designs (e.g., time series, multiple tissues)\n")
cat("- When you want genes expressed across most samples\n")
cat("- NOT recommended for simple two-group comparisons\n\n")

cat("When would you use 'per-group AND'?\n")
cat("- Almost never for DE analysis\n")
cat("- Removes genes that are turned on/off by treatment\n")
cat("- Would miss key differentially expressed genes\n\n")

cat("Reference: DESeq2 paper (Love et al. 2014, Genome Biology)\n")
cat("Quote: 'independent filtering on the mean of normalized counts'\n")
cat("       'for two-group comparisons, pre-filter by group'\n")
