#!/usr/bin/env Rscript

## =========================================================
## DIAGNOSTIC: Gene Symbol to Entrez ID Mapping
## =========================================================
## This script analyzes which genes fail to map during
## pathway analysis and suggests ways to improve mapping rates
##
## Current typical mapping rate: ~80-85%
## Goal: Identify causes of failures and improve to 90%+
## =========================================================

library(dplyr)
library(org.Mm.eg.db)
library(AnnotationDbi)

## Find most recent Part 1 run
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))
tables_dir <- file.path(outdir, "tables")

cat("\n==========================================================\n")
cat("=== DIAGNOSTIC: GENE SYMBOL MAPPING ANALYSIS ===\n")
cat("==========================================================\n")
cat("Run directory: ", outdir, "\n", sep = "")
cat("Run tag: ", run_tag, "\n\n", sep = "")

## Load a representative DE table (use iWAT_F as example)
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("OVERALL", de_files)]

## Use iWAT_F (largest gene set)
iwat_f_file <- de_files[grepl("iWAT_F_KD_vs_CTL", de_files)]

if (length(iwat_f_file) == 0) {
  stop("Could not find iWAT_F DE table")
}

cat("[INFO] Analyzing gene mapping using: ", basename(iwat_f_file), "\n\n", sep = "")

de_table <- read.csv(iwat_f_file, row.names = 1, stringsAsFactors = FALSE)
all_gene_symbols <- rownames(de_table)

cat("==========================================================\n")
cat("=== STEP 1: STANDARD SYMBOL -> ENTREZID MAPPING ===\n")
cat("==========================================================\n\n")

cat("[INFO] Total genes in DE table: ", length(all_gene_symbols), "\n", sep = "")

## Attempt standard mapping (what clusterProfiler does)
mapped_genes <- tryCatch({
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys = all_gene_symbols,
    columns = c("ENTREZID", "GENENAME", "GENETYPE"),
    keytype = "SYMBOL"
  )
}, warning = function(w) {
  cat("[WARN] ", conditionMessage(w), "\n", sep = "")
  return(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = all_gene_symbols,
    columns = c("ENTREZID", "GENENAME", "GENETYPE"),
    keytype = "SYMBOL"
  ))
}, error = function(e) {
  cat("[ERROR] Mapping failed: ", conditionMessage(e), "\n", sep = "")
  return(NULL)
})

## Identify mapped vs unmapped
mapped_symbols <- unique(mapped_genes$SYMBOL[!is.na(mapped_genes$ENTREZID)])
unmapped_symbols <- setdiff(all_gene_symbols, mapped_symbols)

n_mapped <- length(mapped_symbols)
n_unmapped <- length(unmapped_symbols)
pct_mapped <- round(100 * n_mapped / length(all_gene_symbols), 1)

cat("\n--- MAPPING RESULTS ---\n")
cat("Mapped genes:     ", n_mapped, " (", pct_mapped, "%)\n", sep = "")
cat("Unmapped genes:   ", n_unmapped, " (", round(100 - pct_mapped, 1), "%)\n\n", sep = "")

## Check for 1:many mappings (one symbol -> multiple Entrez IDs)
duplicated_symbols <- mapped_genes$SYMBOL[duplicated(mapped_genes$SYMBOL) & !is.na(mapped_genes$ENTREZID)]
n_duplicated <- length(unique(duplicated_symbols))

if (n_duplicated > 0) {
  cat("[WARN] ", n_duplicated, " symbols map to multiple Entrez IDs (ambiguous)\n", sep = "")
  cat("Example: ", head(unique(duplicated_symbols), 5), "\n\n", sep = " ")
}

cat("==========================================================\n")
cat("=== STEP 2: ANALYZE UNMAPPED GENES ===\n")
cat("==========================================================\n\n")

if (n_unmapped > 0) {
  cat("[INFO] Analyzing ", n_unmapped, " unmapped genes...\n\n", sep = "")

  ## Show examples
  cat("--- First 30 unmapped gene symbols ---\n")
  print(head(unmapped_symbols, 30))
  cat("\n")

  ## Categorize by naming pattern
  cat("--- Gene Name Patterns ---\n")

  ## Predicted genes (Gm#####)
  gm_genes <- grep("^Gm[0-9]+$", unmapped_symbols, value = TRUE)
  cat("Predicted genes (Gm#####):        ", length(gm_genes), " (",
      round(100 * length(gm_genes) / n_unmapped, 1), "%)\n", sep = "")

  ## RIKEN genes (e.g., 1110008P14Rik)
  riken_genes <- grep("Rik$", unmapped_symbols, value = TRUE)
  cat("RIKEN cDNA genes (*Rik):          ", length(riken_genes), " (",
      round(100 * length(riken_genes) / n_unmapped, 1), "%)\n", sep = "")

  ## lncRNAs (often have -/. in names or specific patterns)
  lnc_genes <- grep("^[0-9]|^[A-Z][a-z]+[0-9]|\\-", unmapped_symbols, value = TRUE)
  cat("Potential lncRNAs/ncRNAs:         ", length(lnc_genes), " (",
      round(100 * length(lnc_genes) / n_unmapped, 1), "%)\n", sep = "")

  ## Pseudogenes (often have "ps" or "Gm" prefix)
  pseudo_genes <- grep("ps[0-9]*$|Gm", unmapped_symbols, value = TRUE)
  cat("Potential pseudogenes:            ", length(pseudo_genes), " (",
      round(100 * length(pseudo_genes) / n_unmapped, 1), "%)\n", sep = "")

  cat("\n")

  ## Check if unmapped genes are enriched in specific DE categories
  cat("--- Are unmapped genes differentially expressed? ---\n")

  logFC_cut <- 1
  fdr_cut <- 0.05

  unmapped_de <- de_table[unmapped_symbols, ]
  if ("adj.P.Val" %in% colnames(unmapped_de)) {
    fdr_col <- "adj.P.Val"
  } else if ("padj" %in% colnames(unmapped_de)) {
    fdr_col <- "padj"
  } else {
    fdr_col <- NULL
  }

  if (!is.null(fdr_col)) {
    if ("logFC" %in% colnames(unmapped_de)) {
      logfc_col <- "logFC"
    } else if ("log2FoldChange" %in% colnames(unmapped_de)) {
      logfc_col <- "log2FoldChange"
    } else {
      logfc_col <- NULL
    }

    if (!is.null(logfc_col)) {
      unmapped_sig <- sum(
        unmapped_de[[fdr_col]] < fdr_cut &
        abs(unmapped_de[[logfc_col]]) > logFC_cut,
        na.rm = TRUE
      )
      unmapped_total <- nrow(unmapped_de)
      pct_unmapped_sig <- round(100 * unmapped_sig / unmapped_total, 1)

      cat("Unmapped DEGs (FDR<0.05, |logFC|>1): ", unmapped_sig, " / ", unmapped_total,
          " (", pct_unmapped_sig, "%)\n", sep = "")

      ## Compare to mapped genes
      mapped_de <- de_table[mapped_symbols, ]
      mapped_sig <- sum(
        mapped_de[[fdr_col]] < fdr_cut &
        abs(mapped_de[[logfc_col]]) > logFC_cut,
        na.rm = TRUE
      )
      mapped_total <- nrow(mapped_de)
      pct_mapped_sig <- round(100 * mapped_sig / mapped_total, 1)

      cat("Mapped DEGs (FDR<0.05, |logFC|>1):   ", mapped_sig, " / ", mapped_total,
          " (", pct_mapped_sig, "%)\n\n", sep = "")

      if (pct_unmapped_sig > pct_mapped_sig) {
        cat("[WARN] Unmapped genes have HIGHER DE rate - losing significant signal!\n\n")
      } else {
        cat("[OK] Unmapped genes have similar/lower DE rate - less impact on pathways\n\n")
      }
    }
  }
}

cat("==========================================================\n")
cat("=== STEP 3: ALTERNATIVE MAPPING STRATEGIES ===\n")
cat("==========================================================\n\n")

## Strategy 1: Try ALIAS (alternative gene symbols)
cat("--- Strategy 1: Try ALIAS mapping ---\n")
cat("[INFO] Some genes may use outdated/alternative symbols\n")

alias_mapped <- tryCatch({
  AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unmapped_symbols,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ALIAS"
  )
}, warning = function(w) {
  return(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unmapped_symbols,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ALIAS"
  ))
}, error = function(e) {
  cat("[ERROR] ALIAS mapping failed: ", conditionMessage(e), "\n")
  return(NULL)
})

if (!is.null(alias_mapped)) {
  alias_rescued <- unique(alias_mapped$ALIAS[!is.na(alias_mapped$ENTREZID)])
  n_rescued_alias <- length(alias_rescued)
  cat("[OK] ALIAS mapping rescued ", n_rescued_alias, " additional genes (",
      round(100 * n_rescued_alias / n_unmapped, 1), "% of unmapped)\n", sep = "")

  if (n_rescued_alias > 0) {
    cat("Examples of rescued genes:\n")
    print(head(alias_rescued, 10))
    cat("\n")
  }
} else {
  n_rescued_alias <- 0
}

## Strategy 2: Check ENSEMBL IDs
cat("--- Strategy 2: Check if ENSEMBL IDs are available ---\n")

## Check if we have ENSEMBL IDs in the original data
if ("gene_id" %in% colnames(de_table)) {
  cat("[INFO] ENSEMBL gene_id column found in DE table\n")

  ## Try mapping ENSEMBL -> ENTREZ for unmapped genes
  unmapped_ensembl <- de_table[unmapped_symbols, "gene_id"]

  ensembl_mapped <- tryCatch({
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unmapped_ensembl,
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "ENSEMBL"
    )
  }, warning = function(w) {
    return(AnnotationDbi::select(
      org.Mm.eg.db,
      keys = unmapped_ensembl,
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "ENSEMBL"
    ))
  }, error = function(e) {
    cat("[ERROR] ENSEMBL mapping failed: ", conditionMessage(e), "\n")
    return(NULL)
  })

  if (!is.null(ensembl_mapped)) {
    ensembl_rescued <- unique(ensembl_mapped$ENSEMBL[!is.na(ensembl_mapped$ENTREZID)])
    n_rescued_ensembl <- length(ensembl_rescued)
    cat("[OK] ENSEMBL mapping rescued ", n_rescued_ensembl, " additional genes (",
        round(100 * n_rescued_ensembl / n_unmapped, 1), "% of unmapped)\n", sep = "")
  } else {
    n_rescued_ensembl <- 0
  }
} else {
  cat("[INFO] No ENSEMBL gene_id column found in DE table\n")
  n_rescued_ensembl <- 0
}
cat("\n")

## Strategy 3: Check gene types
cat("--- Strategy 3: Check gene biotypes ---\n")

## For mapped genes, check what types they are
if (!is.null(mapped_genes) && "GENETYPE" %in% colnames(mapped_genes)) {
  gene_types <- table(mapped_genes$GENETYPE[!is.na(mapped_genes$ENTREZID)])
  cat("Gene types in MAPPED genes:\n")
  print(head(sort(gene_types, decreasing = TRUE), 10))
  cat("\n")

  ## Count protein-coding
  n_protein_coding <- sum(grepl("protein.coding", names(gene_types)), na.rm = TRUE)
  if (n_protein_coding > 0) {
    cat("[INFO] Protein-coding genes: ", gene_types["protein-coding"], " (",
        round(100 * gene_types["protein-coding"] / sum(gene_types), 1), "%)\n", sep = "")
  }
}
cat("\n")

cat("==========================================================\n")
cat("=== STEP 4: COMBINED MAPPING IMPROVEMENT ===\n")
cat("==========================================================\n\n")

## Calculate potential improvement
total_rescued <- n_rescued_alias + n_rescued_ensembl
potential_mapped <- n_mapped + total_rescued
potential_pct <- round(100 * potential_mapped / length(all_gene_symbols), 1)

cat("--- SUMMARY ---\n")
cat("Original mapping:          ", n_mapped, " / ", length(all_gene_symbols),
    " (", pct_mapped, "%)\n", sep = "")
cat("Rescued by ALIAS:          +", n_rescued_alias, "\n", sep = "")
cat("Rescued by ENSEMBL:        +", n_rescued_ensembl, "\n", sep = "")
cat("POTENTIAL total mapping:   ", potential_mapped, " / ", length(all_gene_symbols),
    " (", potential_pct, "%)\n", sep = "")
cat("POTENTIAL improvement:     +", round(potential_pct - pct_mapped, 1), " percentage points\n\n", sep = "")

cat("==========================================================\n")
cat("=== RECOMMENDATIONS ===\n")
cat("==========================================================\n\n")

## Recommendation based on results
if (potential_pct > pct_mapped + 5) {
  cat("✅ RECOMMENDED: Implement improved mapping strategy\n\n")
  cat("OPTION 1: Use multi-strategy mapping in pathway scripts\n")
  cat("  1. Try SYMBOL -> ENTREZID (primary)\n")
  cat("  2. For failures, try ALIAS -> ENTREZID\n")
  cat("  3. For failures, try ENSEMBL -> ENTREZID\n")
  cat("  Expected improvement: +", round(potential_pct - pct_mapped, 1), " percentage points\n\n", sep = "")

  cat("OPTION 2: Pre-filter non-coding genes\n")
  cat("  - Remove predicted genes (Gm#####) before pathway analysis\n")
  cat("  - Remove RIKEN cDNAs (*Rik) unless they're DEGs\n")
  cat("  - Focus on protein-coding genes for cleaner pathway results\n\n")

  cat("OPTION 3: Update annotation database\n")
  cat("  - Check if org.Mm.eg.db package is up to date\n")
  cat("  - Consider using more recent annotation versions\n\n")

} else if (n_unmapped > 0) {
  cat("ℹ️  ANALYSIS: Unmapped genes are mostly non-coding\n\n")
  cat("Most unmapped genes appear to be:\n")
  cat("  - Predicted genes (Gm#####)\n")
  cat("  - Long non-coding RNAs\n")
  cat("  - Pseudogenes\n")
  cat("  - RIKEN cDNAs\n\n")
  cat("These genes typically DON'T have Entrez IDs because:\n")
  cat("  - They're not protein-coding\n")
  cat("  - They're not well-characterized\n")
  cat("  - They're not in standard pathway databases\n\n")
  cat("RECOMMENDATION: Current mapping rate (~80-85%) is EXPECTED and ACCEPTABLE\n")
  cat("  - Pathway databases focus on protein-coding genes\n")
  cat("  - Non-coding genes rarely appear in functional pathways\n")
  cat("  - Your current approach is standard practice\n\n")
}

## Check if unmapped genes are important DEGs
if (!is.null(fdr_col) && !is.null(logfc_col)) {
  unmapped_important <- unmapped_de %>%
    dplyr::filter(
      .data[[fdr_col]] < 0.01,
      abs(.data[[logfc_col]]) > 2
    )

  if (nrow(unmapped_important) > 10) {
    cat("⚠️  WARNING: You have ", nrow(unmapped_important), " highly significant unmapped DEGs\n", sep = "")
    cat("   (FDR < 0.01, |logFC| > 2)\n\n")
    cat("These genes may be biologically important even if they don't map to pathways.\n")
    cat("Consider:\n")
    cat("  1. Manual literature search for top unmapped DEGs\n")
    cat("  2. Alternative functional annotation (GO enrichment using gene names)\n")
    cat("  3. Reporting these as 'novel' or 'uncharacterized' genes of interest\n\n")

    cat("Top unmapped DEGs by significance:\n")
    unmapped_sorted <- unmapped_important %>%
      dplyr::arrange(.data[[fdr_col]]) %>%
      dplyr::select(all_of(c(logfc_col, fdr_col)))
    rownames(unmapped_sorted) <- rownames(unmapped_important)[order(unmapped_important[[fdr_col]])]
    print(head(unmapped_sorted, 15))
    cat("\n")
  }
}

cat("==========================================================\n")
cat("=== DIAGNOSTIC COMPLETE ===\n")
cat("==========================================================\n\n")

## Save detailed mapping table
mapping_summary <- data.frame(
  Gene_Symbol = all_gene_symbols,
  Mapped = all_gene_symbols %in% mapped_symbols,
  In_DE_Set = TRUE,
  stringsAsFactors = FALSE
)

if (!is.null(fdr_col) && !is.null(logfc_col)) {
  mapping_summary$logFC <- de_table[all_gene_symbols, logfc_col]
  mapping_summary$FDR <- de_table[all_gene_symbols, fdr_col]
  mapping_summary$Is_DEG <- mapping_summary$FDR < 0.05 & abs(mapping_summary$logFC) > 1
}

summary_file <- paste0("DIAGNOSTIC_gene_mapping_", run_tag, ".csv")
write.csv(mapping_summary, file = file.path(outdir, "tables", summary_file), row.names = FALSE)
cat("Saved detailed mapping table to: ", file.path(outdir, "tables", summary_file), "\n\n", sep = "")

cat("Key findings saved. Review recommendations above.\n")
