#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 2B: ORA WITH RELAXED PARAMETERS
## =========================================================
## This is a RELAXED version for GWAT contrasts only
## Use this if standard parameters give too few pathways
##
## RELAXED PARAMETERS:
## - FDR cutoff: 0.1 (instead of 0.05)
## - logFC cutoff: 0.5 (instead of 1.0)
## - Pathway p-value: 0.15 (instead of 0.1)
## - Pathway q-value: 0.25 (instead of 0.2)
## =========================================================

## Source the original script setup
source("part2_ora.R", local = TRUE)

cat("\n\n==========================================================\n")
cat("=== RELAXED ORA ANALYSIS (GWAT ONLY) ===\n")
cat("==========================================================\n\n")

## Override parameters
logFC_cut_relaxed <- 0.5   # More lenient
fdr_cut_relaxed   <- 0.1   # More lenient

cat("RELAXED PARAMETERS:\n")
cat("  DEG thresholds: FDR < ", fdr_cut_relaxed, ", |logFC| > ", logFC_cut_relaxed, "\n", sep = "")
cat("  Pathway enrichment: p < 0.15, q < 0.25\n\n")

## Re-run ORA for GWAT contrasts only with relaxed parameters
gwat_de_files <- de_files[grepl("gWAT", de_files)]

for (de_file in gwat_de_files) {

  de_filename <- basename(de_file)
  contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

  cat("\n=== RELAXED ORA: ", contrast_name, " ===\n", sep = "")

  de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  if (!"adj.P.Val" %in% colnames(de_table)) {
    if ("padj" %in% colnames(de_table)) {
      de_table$adj.P.Val <- de_table$padj
    } else {
      cat("[WARN] Cannot find FDR column\n")
      next
    }
  }

  if (!"logFC" %in% colnames(de_table)) {
    if ("log2FoldChange" %in% colnames(de_table)) {
      de_table$logFC <- de_table$log2FoldChange
    } else {
      cat("[WARN] Cannot find logFC column\n")
      next
    }
  }

  ## Get genes with RELAXED cutoffs
  up_genes <- de_table %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_relaxed, logFC > logFC_cut_relaxed) %>%
    rownames()

  down_genes <- de_table %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_relaxed, logFC < -logFC_cut_relaxed) %>%
    rownames()

  cat("[INFO] Up-regulated DEGs (relaxed): ", length(up_genes), "\n", sep = "")
  cat("[INFO] Down-regulated DEGs (relaxed): ", length(down_genes), "\n", sep = "")

  ## Run ORA with RELAXED pathway parameters
  if (length(down_genes) >= 5) {
    cat("\n--- Running RELAXED ORA for Down genes ---\n")
    tryCatch({
      gene_entrez <- bitr(
        unique(down_genes),
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db,
        drop     = TRUE
      )

      if (nrow(gene_entrez) >= 3) {
        ## GO:BP with relaxed parameters
        enrich_go_relaxed <- enrichGO(
          gene          = gene_entrez$ENTREZID,
          OrgDb         = org.Mm.eg.db,
          ont           = "BP",
          pvalueCutoff  = 0.15,   # RELAXED (was 0.1)
          qvalueCutoff  = 0.25,   # RELAXED (was 0.2)
          readable      = TRUE,
          minGSSize     = 3,      # RELAXED (was 5)
          maxGSSize     = 500,
          universe      = if (use_universe) universe_entrez else NULL
        )

        if (!is.null(enrich_go_relaxed) && nrow(enrich_go_relaxed@result) > 0) {
          enrich_go_simp <- simplify(enrich_go_relaxed, cutoff = 0.7, by = "p.adjust", select_fun = min)

          sig_results <- enrich_go_simp@result %>%
            dplyr::filter(p.adjust < 0.25) %>%  # RELAXED
            dplyr::arrange(p.adjust)

          cat("[OK] RELAXED GO:BP: ", nrow(sig_results), " significant pathways\n", sep = "")

          if (nrow(sig_results) > 0) {
            table_file <- paste0("ORA_RELAXED_gobp_", contrast_name, "_Down_", run_tag, ".csv")
            write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
            cat("[OK] Saved: ", table_file, "\n", sep = "")
          }
        }
      }
    }, error = function(e) {
      cat("[ERROR] ", conditionMessage(e), "\n", sep = "")
    })
  }
}

cat("\n=== RELAXED ORA COMPLETE ===\n\n")
