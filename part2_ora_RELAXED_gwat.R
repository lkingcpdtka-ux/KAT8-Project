#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - RELAXED ORA FOR GWAT ONLY
## =========================================================
## This script runs ORA with RELAXED parameters specifically
## for GWAT contrasts to find more down-regulated pathways
##
## RELAXED PARAMETERS:
## - DEG FDR cutoff: 0.1 (instead of 0.05)
## - DEG logFC cutoff: 0.5 (instead of 1.0)
## - Pathway p-value: 0.15 (instead of 0.1)
## - Pathway q-value: 0.25 (instead of 0.2)
## - minGSSize: 3 (instead of 5)
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "DOSE",
  "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(AnnotationDbi)
})

## 2) Find most recent Part 1 run --------------------------
cat("\n=== RELAXED ORA FOR GWAT CONTRASTS ===\n")
cat("Using relaxed parameters to find more down-regulated pathways\n\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run part1_main_analysis.R first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No Part 1 run directories found in savepoints/")
}

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using Part 1 results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n\n", sep = "")

tables_dir <- file.path(outdir, "tables")
if (!dir.exists(tables_dir)) {
  stop("tables directory not found in Part 1 run: ", outdir)
}
if (!dir.exists(file.path(outdir, "plots"))) {
  stop("plots directory not found in Part 1 run: ", outdir)
}

## 3) RELAXED Parameters ------------------------------------
cat("=== RELAXED PARAMETERS ===\n")
logFC_cut_relaxed <- 0.5   # More lenient (was 1.0)
fdr_cut_relaxed   <- 0.1   # More lenient (was 0.05)

cat("  DEG thresholds: FDR < ", fdr_cut_relaxed, ", |logFC| > ", logFC_cut_relaxed, "\n", sep = "")
cat("  Pathway p-cutoff: 0.15 (was 0.1)\n")
cat("  Pathway q-cutoff: 0.25 (was 0.2)\n")
cat("  minGSSize: 3 (was 5)\n\n")

## 4) Load DE tables ----------------------------------------
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("OVERALL", de_files)]

## Filter for GWAT only
gwat_de_files <- de_files[grepl("gWAT", de_files)]

if (length(gwat_de_files) == 0) {
  stop("No GWAT DE tables found")
}

cat("[INFO] Found ", length(gwat_de_files), " GWAT contrasts\n\n")

## 5) Build universe (optional) -----------------------------
use_universe <- TRUE
universe_entrez <- NULL

if (use_universe) {
  cat("[INFO] Building background gene set (universe) from all tested genes...\n")

  ## Use first DE file to get all tested genes
  de_ref <- read.csv(gwat_de_files[1], row.names = 1, stringsAsFactors = FALSE)
  all_genes_tested <- rownames(de_ref)

  universe_genes_mapped <- tryCatch({
    bitr(
      unique(all_genes_tested),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Mm.eg.db,
      drop     = TRUE
    )
  }, error = function(e) {
    cat("[WARN] Failed to map universe genes: ", conditionMessage(e), "\n")
    return(NULL)
  })

  if (!is.null(universe_genes_mapped) && nrow(universe_genes_mapped) > 0) {
    universe_entrez <- unique(universe_genes_mapped$ENTREZID)
    cat("[OK] Universe built: ", length(universe_entrez), " genes\n\n", sep = "")
  } else {
    cat("[WARN] Could not build universe; will use default\n\n")
    use_universe <- FALSE
  }
} else {
  cat("[INFO] NOT using background gene set (default/liberal approach)\n\n")
}

## 6) Run relaxed ORA for GWAT ------------------------------
for (de_file in gwat_de_files) {

  de_filename <- basename(de_file)
  contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

  cat("\n=== RELAXED ORA: ", contrast_name, " ===\n", sep = "")

  ## Load DE table
  de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Standardize column names
  if (!"adj.P.Val" %in% colnames(de_table)) {
    if ("padj" %in% colnames(de_table)) {
      de_table$adj.P.Val <- de_table$padj
    } else {
      cat("[WARN] Cannot find FDR column in ", de_file, "\n")
      next
    }
  }

  if (!"logFC" %in% colnames(de_table)) {
    if ("log2FoldChange" %in% colnames(de_table)) {
      de_table$logFC <- de_table$log2FoldChange
    } else {
      cat("[WARN] Cannot find logFC column in ", de_file, "\n")
      next
    }
  }

  ## Get up and down genes with RELAXED cutoffs
  up_genes_relaxed <- de_table %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_relaxed, logFC > logFC_cut_relaxed) %>%
    rownames()

  down_genes_relaxed <- de_table %>%
    dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_relaxed, logFC < -logFC_cut_relaxed) %>%
    rownames()

  cat("[INFO] Up DEGs (relaxed): ", length(up_genes_relaxed), "\n", sep = "")
  cat("[INFO] Down DEGs (relaxed): ", length(down_genes_relaxed), "\n", sep = "")

  ## Focus on DOWN genes (since that's the issue)
  if (length(down_genes_relaxed) >= 5) {
    cat("\n--- Running RELAXED ORA for DOWN genes ---\n")

    ## Convert to Entrez IDs
    gene_entrez <- tryCatch({
      bitr(
        unique(down_genes_relaxed),
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db,
        drop     = TRUE
      )
    }, error = function(e) {
      cat("[ERROR] Gene ID conversion failed: ", conditionMessage(e), "\n")
      return(data.frame(SYMBOL = character(), ENTREZID = character()))
    })

    if (nrow(gene_entrez) >= 3) {
      cat("[OK] Mapped ", nrow(gene_entrez), " of ", length(down_genes_relaxed), " genes to Entrez (",
          round(100 * nrow(gene_entrez) / length(down_genes_relaxed), 1), "%)\n", sep = "")

      ## GO:BP enrichment with RELAXED parameters
      cat("[INFO] Running GO:BP enrichment (relaxed)...\n")
      enrich_go_relaxed <- tryCatch({
        if (use_universe) {
          enrichGO(
            gene          = gene_entrez$ENTREZID,
            OrgDb         = org.Mm.eg.db,
            ont           = "BP",
            pvalueCutoff  = 0.15,   # RELAXED (was 0.1)
            qvalueCutoff  = 0.25,   # RELAXED (was 0.2)
            readable      = TRUE,
            minGSSize     = 3,      # RELAXED (was 5)
            maxGSSize     = 500,
            universe      = universe_entrez
          )
        } else {
          enrichGO(
            gene          = gene_entrez$ENTREZID,
            OrgDb         = org.Mm.eg.db,
            ont           = "BP",
            pvalueCutoff  = 0.15,
            qvalueCutoff  = 0.25,
            readable      = TRUE,
            minGSSize     = 3,
            maxGSSize     = 500
          )
        }
      }, error = function(e) {
        cat("[ERROR] GO:BP enrichment failed: ", conditionMessage(e), "\n")
        return(NULL)
      })

      ## Process GO:BP results
      if (!is.null(enrich_go_relaxed) && nrow(enrich_go_relaxed@result) > 0) {
        ## Simplify to remove redundant terms
        enrich_go_simp <- tryCatch({
          simplify(enrich_go_relaxed, cutoff = 0.7, by = "p.adjust", select_fun = min)
        }, error = function(e) {
          cat("[WARN] Simplification failed, using unsimplified results\n")
          return(enrich_go_relaxed)
        })

        sig_results <- enrich_go_simp@result %>%
          dplyr::filter(p.adjust < 0.25) %>%  # RELAXED
          dplyr::arrange(p.adjust)

        cat("[OK] RELAXED GO:BP: ", nrow(sig_results), " significant pathways\n", sep = "")

        if (nrow(sig_results) > 0) {
          ## Save table
          table_file <- paste0("ORA_RELAXED_gobp_", contrast_name, "_Down_", run_tag, ".csv")
          write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
          cat("[OK] Saved: ", table_file, "\n", sep = "")

          ## Create bar plot (top 20)
          if (nrow(sig_results) > 0) {
            plot_data <- sig_results %>%
              dplyr::arrange(p.adjust) %>%
              dplyr::slice_head(n = 20) %>%
              dplyr::mutate(
                Description = factor(Description, levels = rev(Description)),
                log10p = -log10(p.adjust)
              )

            p_bar <- ggplot(plot_data, aes(x = log10p, y = Description)) +
              geom_bar(stat = "identity", fill = "#0072B2", color = "black", linewidth = 0.3) +
              labs(
                title = paste0("RELAXED ORA GO:BP: ", contrast_name, " (Down)"),
                x = expression(-log[10]~FDR),
                y = NULL
              ) +
              theme_classic(base_size = 12) +
              theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                axis.text.y = element_text(size = 10, color = "black"),
                axis.text.x = element_text(size = 10, color = "black", face = "bold"),
                axis.title.x = element_text(size = 12, face = "bold"),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
              )

            plot_file <- paste0("ORA_RELAXED_barplot_gobp_", contrast_name, "_Down_", run_tag, ".png")
            ggsave(file.path(outdir, "plots", plot_file), plot = p_bar,
                   width = 12, height = max(6, nrow(plot_data) * 0.4), dpi = 300, bg = "white")
            cat("[OK] Saved plot: ", plot_file, "\n", sep = "")
          }
        }
      } else {
        cat("[INFO] No significant GO:BP pathways even with relaxed parameters\n")
      }

      ## KEGG enrichment with RELAXED parameters
      cat("[INFO] Running KEGG enrichment (relaxed)...\n")
      enrich_kegg_relaxed <- tryCatch({
        if (use_universe) {
          enrichKEGG(
            gene         = gene_entrez$ENTREZID,
            organism     = "mmu",
            pvalueCutoff = 0.15,
            qvalueCutoff = 0.25,
            minGSSize    = 3,
            maxGSSize    = 500,
            universe     = universe_entrez
          )
        } else {
          enrichKEGG(
            gene         = gene_entrez$ENTREZID,
            organism     = "mmu",
            pvalueCutoff = 0.15,
            qvalueCutoff = 0.25,
            minGSSize    = 3,
            maxGSSize    = 500
          )
        }
      }, error = function(e) {
        cat("[ERROR] KEGG enrichment failed: ", conditionMessage(e), "\n")
        return(NULL)
      })

      ## Process KEGG results
      if (!is.null(enrich_kegg_relaxed) && nrow(enrich_kegg_relaxed@result) > 0) {
        ## Make readable
        enrich_kegg_readable <- setReadable(enrich_kegg_relaxed, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

        sig_results_kegg <- enrich_kegg_readable@result %>%
          dplyr::filter(p.adjust < 0.25) %>%
          dplyr::arrange(p.adjust)

        cat("[OK] RELAXED KEGG: ", nrow(sig_results_kegg), " significant pathways\n", sep = "")

        if (nrow(sig_results_kegg) > 0) {
          ## Save table
          table_file_kegg <- paste0("ORA_RELAXED_kegg_", contrast_name, "_Down_", run_tag, ".csv")
          write.csv(sig_results_kegg, file = file.path(outdir, "tables", table_file_kegg), row.names = FALSE)
          cat("[OK] Saved: ", table_file_kegg, "\n", sep = "")

          ## Create bar plot (top 20)
          if (nrow(sig_results_kegg) > 0) {
            plot_data_kegg <- sig_results_kegg %>%
              dplyr::arrange(p.adjust) %>%
              dplyr::slice_head(n = 20) %>%
              dplyr::mutate(
                Description = factor(Description, levels = rev(Description)),
                log10p = -log10(p.adjust)
              )

            p_bar_kegg <- ggplot(plot_data_kegg, aes(x = log10p, y = Description)) +
              geom_bar(stat = "identity", fill = "#0072B2", color = "black", linewidth = 0.3) +
              labs(
                title = paste0("RELAXED ORA KEGG: ", contrast_name, " (Down)"),
                x = expression(-log[10]~FDR),
                y = NULL
              ) +
              theme_classic(base_size = 12) +
              theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                axis.text.y = element_text(size = 10, color = "black"),
                axis.text.x = element_text(size = 10, color = "black", face = "bold"),
                axis.title.x = element_text(size = 12, face = "bold"),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
              )

            plot_file_kegg <- paste0("ORA_RELAXED_barplot_kegg_", contrast_name, "_Down_", run_tag, ".png")
            ggsave(file.path(outdir, "plots", plot_file_kegg), plot = p_bar_kegg,
                   width = 12, height = max(6, nrow(plot_data_kegg) * 0.4), dpi = 300, bg = "white")
            cat("[OK] Saved plot: ", plot_file_kegg, "\n", sep = "")
          }
        }
      } else {
        cat("[INFO] No significant KEGG pathways even with relaxed parameters\n")
      }

    } else {
      cat("[WARN] Too few genes mapped to Entrez IDs\n")
    }
  } else {
    cat("[WARN] Too few down-regulated genes even with relaxed cutoffs\n")
  }
}

cat("\n======================================\n")
cat("=== RELAXED ORA COMPLETE ===\n")
cat("======================================\n\n")

cat("NOTE: These results use MORE LENIENT thresholds than standard analysis.\n")
cat("      - Higher false positive rate expected\n")
cat("      - Use these as exploratory/hypothesis-generating results\n")
cat("      - Compare with standard ORA results to assess robustness\n\n")

cat("=== Part 2 (RELAXED ORA) finished ===\n")
