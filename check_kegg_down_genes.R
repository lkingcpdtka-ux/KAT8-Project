#!/usr/bin/env Rscript
## =========================================================
## Diagnostic: Check why KEGG Down genes have no significant pathways
## =========================================================

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)

## Find the most recent run
savepoint_dir <- "savepoints"
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]

cat("Using results from:", outdir, "\n\n")

## Load DE results
de_file <- list.files(file.path(outdir, "tables"), pattern = "^DE_cells_.*\\.csv$", full.names = TRUE)[1]
tt <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

## Get down genes
down_genes <- tt %>%
  dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange < -0.5) %>%
  rownames()

cat("Down-regulated genes:", length(down_genes), "\n")
cat("First 10:", paste(head(down_genes, 10), collapse = ", "), "\n\n")

## Convert to Entrez
gene_entrez <- bitr(unique(down_genes), fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db, drop = TRUE)
entrez_ids <- gene_entrez$ENTREZID

cat("Mapped to Entrez:", length(entrez_ids), "\n\n")

## Run KEGG enrichment WITHOUT strict filtering
cat("=== Running KEGG enrichment (relaxed thresholds) ===\n")
enrich_kegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "mmu",
  pvalueCutoff = 1.0,  # Show ALL results
  qvalueCutoff = 1.0,
  minGSSize = 5,
  maxGSSize = 500
)

if (!is.null(enrich_kegg) && nrow(enrich_kegg@result) > 0) {
  results <- enrich_kegg@result %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count)

  cat("\n=== TOP 20 KEGG PATHWAYS (by p-value) ===\n")
  print(head(results, 20))

  cat("\n=== PATHWAYS CLOSE TO SIGNIFICANCE (p.adjust < 0.1) ===\n")
  close <- results %>% dplyr::filter(p.adjust < 0.1)
  if (nrow(close) > 0) {
    print(close)
  } else {
    cat("None with p.adjust < 0.1\n")
    cat("\nBest p.adjust:", min(results$p.adjust), "\n")
  }

  cat("\n=== RECOMMENDATION ===\n")
  if (min(results$p.adjust) > 0.1) {
    cat("No KEGG pathways are even close to significant (best p.adjust = ",
        round(min(results$p.adjust), 3), ")\n", sep = "")
    cat("This suggests your down-regulated genes are NOT concentrated in specific KEGG pathways.\n")
    cat("Consider:\n")
    cat("  1. Use FGSEA instead (more sensitive for dispersed signals)\n")
    cat("  2. Focus on GO:BP results (broader categories)\n")
    cat("  3. Look at individual gene functions rather than pathways\n")
  } else {
    cat("Some pathways are close. Consider relaxing FDR threshold to 0.1\n")
  }

} else {
  cat("No KEGG enrichment results found!\n")
}
