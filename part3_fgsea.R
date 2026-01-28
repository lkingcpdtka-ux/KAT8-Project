#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 3: GSEA PATHWAY ANALYSIS
## =========================================================
## This script performs Gene Set Enrichment Analysis (gseGO/gseKEGG)
## for GO:BP and KEGG pathways using ranked gene lists from Part 1
##
## FIXES APPLIED:
## - Use DESeq2 Wald statistic for ranking
## - Align GO/KEGG GSEA via clusterProfiler
## - Better error handling
## - Direction-specific pathway bar plots
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "org.Mm.eg.db",
  "AnnotationDbi",
  "GO.db",
  "clusterProfiler",
  "GOSemSim"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(clusterProfiler)
})

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  source(file.path(utils_dir, "purge.R"))
  source(file.path(utils_dir, "dedupe.R"))
  use_save_core <- TRUE
} else {
  cat("[WARN] save_core not found; proceeding without it\n")
  use_save_core <- FALSE
}

## 2.5) Find most recent Part 1 run and use same directory ----
cat("\n=== FINDING PART 1 RESULTS ===\n")

## Look for run directories in savepoints
savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run part1_main_analysis.R first.")
}

## Find most recent run directory
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No Part 1 run directories found in savepoints/")
}

## Sort by modification time and get most recent
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]

## Extract run tag from directory name
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using Part 1 results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

## Verify directory structure exists
if (!dir.exists(file.path(outdir, "tables"))) {
  stop("tables directory not found in Part 1 run: ", outdir)
}
if (!dir.exists(file.path(outdir, "plots"))) {
  stop("plots directory not found in Part 1 run: ", outdir)
}

## 2.6) Load authoritative DEG lists from Part 1 -----------
## This ensures fgsea sanity checks can be cross-validated with ORA
## Note: GSEA uses ALL ranked genes for analysis, but sanity checks
## should report the authoritative DEG counts for cross-validation
cat("\n=== LOADING AUTHORITATIVE DEG LISTS FOR VALIDATION ===\n")
deg_lists_file <- list.files(file.path(outdir, "tables"), pattern = "^DEG_lists_authoritative_.*\\.rds$", full.names = TRUE)
if (length(deg_lists_file) == 0) {
  cat("[WARN] No authoritative DEG lists found\n")
  cat("[WARN] For full cross-validation, re-run part1_main_analysis.R first\n")
  authoritative_deg_lists <- NULL
} else {
  deg_lists_file <- deg_lists_file[order(file.info(deg_lists_file)$mtime, decreasing = TRUE)][1]
  authoritative_deg_lists <- readRDS(deg_lists_file)
  cat("[OK] Loaded authoritative DEG lists from: ", basename(deg_lists_file), "\n", sep = "")
  cat("[INFO] Contrasts available: ", paste(names(authoritative_deg_lists), collapse = ", "), "\n", sep = "")
  for (cn in names(authoritative_deg_lists)) {
    cat("  ", cn, ": Up=", length(authoritative_deg_lists[[cn]]$up),
        ", Down=", length(authoritative_deg_lists[[cn]]$down),
        " (FDR<", authoritative_deg_lists[[cn]]$fdr_cutoff,
        ", |logFC|>", authoritative_deg_lists[[cn]]$logFC_cutoff, ")\n", sep = "")
  }
}

## 3) Parameters (from central parameters.R) ----------------
## All cutoffs are now defined in parameters.R
fdr_cut   <- gsea_params$fdr_cutoff
simplify_go_bp <- gsea_params$simplify_go
simplify_go_cutoff <- gsea_params$simplify_cutoff
top_n_per_direction <- gsea_params$top_n_per_direction

## Local params that extend the central ones
local_gsea_params <- list(
  min_gs_size  = gsea_params$min_gs_size,
  max_gs_size  = gsea_params$max_gs_size,
  rank_metric  = gsea_params$rank_metric
)

gse_kegg_params <- list(
  min_gs_size     = gsea_params$min_gs_size,
  max_gs_size     = gsea_params$max_gs_size,
  pvalue_cutoff   = 0.1,  ## Initial filter (final uses fdr_cut)
  p_adjust_method = "BH",
  organism        = "mmu",
  seed            = TRUE
)

## ============================================================
## TROUBLESHOOTING: If GWAT has few down-regulated pathways
## ============================================================
## Try adjusting these parameters for GWAT contrasts:
##
## 1. RELAX FDR CUTOFF (for pathway significance):
##    fdr_cut <- 0.1  # More lenient (was 0.05)
##
## 2. RELAX GENE SET SIZES (allow smaller pathways):
##    In gseGO/gseKEGG calls below, change:
##    minGSSize = 3  # Was 5
##
## 3. USE DIFFERENT RANKING METRIC:
##    Instead of Wald stat, try: -log10(pvalue) * sign(logFC)
##
## Note: GSEA doesn't need hard DEG cutoffs - it uses the full
## ranked gene list, so it's MORE SENSITIVE than ORA for finding
## pathways with subtle shifts.
## ============================================================

## Library toggles
## UPDATED (2026-01-20): Re-enabled GO:BP and KEGG for comprehensive pathway analysis
## These were previously disabled, which limited pathway discovery
run_go_bp <- TRUE         # GO Biological Process (large, comprehensive)
run_kegg <- TRUE          # KEGG pathways (curated, metabolic focus)
run_wikipathways <- FALSE # WikiPathways (community-curated)
run_hallmark <- FALSE     # MSigDB Hallmark (50 well-defined signatures)

## ============================================================
## COMPLETE GSEA PARAMETER DOCUMENTATION
## ============================================================
## This section documents all GSEA parameters for reference
##
## METHOD: Gene Set Enrichment Analysis (clusterProfiler)
##   - Rank-based enrichment test
##   - Tests if pathway genes are enriched at top/bottom of ranked list
##   - Uses ALL genes (not just DEGs) - more sensitive than ORA
##   - Calculates Enrichment Score (ES) and Normalized ES (NES)
##
## GENE RANKING METRIC:
##   - DESeq2 Wald statistic (RECOMMENDED)
##   - Wald stat = log2FC / lfcSE (combines effect size + significance)
##   - Positive = up-regulated, Negative = down-regulated
##   - Alternative: -log10(pvalue) * sign(logFC) (mentioned in troubleshooting)
##
## GSEA CORE PARAMETERS:
##   - minGSSize: 5 (minimum genes in pathway, smaller pathways excluded)
##   - maxGSSize: 500 (maximum genes in pathway, larger pathways excluded)
##
## SIGNIFICANCE THRESHOLDS:
##   - FDR cutoff: 0.05 (Benjamini-Hochberg adjusted p-value)
##   - No logFC threshold (GSEA uses continuous ranking, not cutoffs)
##
## DATABASES ANALYZED:
##   - GO:BP: Gene Ontology Biological Process
##     * Source: org.Mm.eg.db (via AnnotationDbi)
##     * Size: ~17,000+ terms (large, comprehensive)
##     * Simplification: Yes (cutoff=0.7, semantic similarity via GOSemSim)
##   - KEGG: Kyoto Encyclopedia of Genes and Genomes
##     * Source: clusterProfiler::gseKEGG()
##     * Matches ORA's enrichKEGG() for consistency
##     * Requires ENTREZID mapping
##
## GO:BP SIMPLIFICATION (OPTIONAL):
##   - Enabled: simplify_go_bp = TRUE
##   - Cutoff: 0.7 (semantic similarity threshold)
##   - Method: GOSemSim::simplify() from clusterProfiler
##   - Purpose: Reduce redundancy in GO:BP results
##   - Keeps most significant term among similar ones
##
## VISUALIZATION:
##   - Bar plots: Top N pathways per direction (up/down-regulated)
##   - top_n_per_direction: 15 (shows top 15 up + top 15 down)
##   - X-axis: NES (Normalized Enrichment Score)
##     * Positive NES = pathway enriched in up-regulated genes
##     * Negative NES = pathway enriched in down-regulated genes
##   - Color: Direction-specific (orange=up, blue=down)
##
## OUTPUT FILES:
##   - fgsea_[database]_[contrast].csv: All significant pathways (FDR<0.05)
##   - fgsea_leadingEdge_[database]_[contrast].csv: Core enrichment genes
##   - fgsea_plot_[database]_[contrast].png: Bar plot visualization
##
## LEADING EDGE:
##   - Core subset of pathway genes driving enrichment signal
##   - Genes that appear before running ES peak in ranked list
##   - Saved as semicolon-separated strings in CSV
##
## SANITY TRACKING:
##   - fgsea_sanity_tracker logs: contrast, database,
##     N_Input_Genes, N_Pathways_Tested, N_Sig_Pathways,
##     N_Sig_Up, N_Sig_Down
##
## KEY DIFFERENCES FROM ORA:
##   - Uses continuous gene ranking (not discrete DEG list)
##   - Tests enrichment across entire ranked list
##   - More sensitive to subtle coordinated changes
##   - Doesn't require hard FDR/logFC cutoffs for input
##   - Can detect pathways with many small changes
## ============================================================

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## Sanity check tracker - includes DEG counts for cross-validation with ORA
## Note: N_Input_Genes = ranked genes for GSEA (all genes)
##       N_DEG_Up/N_DEG_Down = authoritative DEG counts (for cross-validation)
fgsea_sanity_tracker <- data.frame(
  Contrast = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_DEG_Up = integer(),
  N_DEG_Down = integer(),
  N_Pathways_Tested = integer(),
  N_Sig_Pathways = integer(),
  N_Sig_Up = integer(),
  N_Sig_Down = integer(),
  logFC_cutoff = numeric(),
  fdr_cutoff = numeric(),
  stringsAsFactors = FALSE
)

## 4) Main logic --------------------------------------------
tryCatch({
  
  ## 4.1) Load DE tables from Part 1 ------------------------
  cat("\n=== LOADING DE TABLES FROM PART 1 ===\n")
  
  ## Find all DE tables
  tables_dir <- file.path(outdir, "tables")
  if (!dir.exists(tables_dir)) {
    stop("tables directory not found in Part 1 run: ", outdir)
  }
  
  de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
  if (length(de_files) == 0) {
    stop("No DE tables found in: ", tables_dir)
  }
  
  cat("[INFO] Found ", length(de_files), " DE tables\n")
  
  ## 4.2) Build gene sets ------------------------------------
  cat("\n=== BUILDING GENE SETS ===\n")
  log_time("Starting gene set construction")
  
  ## GO:BP gene sets (disabled for now)
  go_bp_list <- NULL
  go_bp_term2gene <- NULL
  if (run_go_bp) {
    cat("[INFO] Building GO:BP gene sets from org.Mm.eg.db...\n")
    go_bp_genes <- tryCatch({
      AnnotationDbi::select(
        org.Mm.eg.db,
        keys = keys(org.Mm.eg.db, keytype = "GOALL"),
        columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
        keytype = "GOALL"
      )
    }, error = function(e) {
      cat("[ERROR] Failed to get GO:BP gene sets: ", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(go_bp_genes)) {
      ## Filter for BP ontology and remove NA symbols
      go_bp_genes <- go_bp_genes %>%
        dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
        dplyr::select(GOALL, SYMBOL)
      
      ## Convert to named list (pathway name -> gene vector)
      go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)
      go_bp_term2gene <- go_bp_genes %>% dplyr::select(GOALL, SYMBOL)
      
      ## Filter for gene set size
      go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= 5 & sapply(go_bp_list, length) <= 500]
      
      cat("[OK] GO:BP gene sets: ", length(go_bp_list), " pathways\n", sep = "")
    }
  } else {
    cat("[INFO] GO:BP gene sets disabled (run_go_bp = FALSE)\n")
  }
  
  ## ============================================================
  ## KEGG GSEA: Now uses clusterProfiler::gseKEGG() instead of fgsea+msigdbr
  ## This ensures consistency with ORA (which uses enrichKEGG)
  ## No pre-built gene sets needed - gseKEGG queries KEGG API directly
  ## ============================================================
  if (run_kegg) {
    cat("[INFO] KEGG GSEA will use clusterProfiler::gseKEGG() (not msigdbr)\n")
    cat("[INFO] This matches ORA's enrichKEGG() for consistency\n")
  } else {
    cat("[INFO] KEGG GSEA disabled (run_kegg = FALSE)\n")
  }
  
  ## WikiPathways/Hallmark gene sets are disabled in this workflow.
  log_time("Completed gene set construction")
  
  ## 4.3) GSEA analysis function ----------------------------
  ## NOTE: GO:BP uses gseGO and KEGG uses gseKEGG (clusterProfiler)
  run_fgsea_analysis <- function(ranked_genes, contrast_name,
                                 go_bp_list, go_bp_term2gene, wikipathways_list, hallmark_list,
                                 de_table,  ## Added: full DE table for ENTREZID mapping
                                 run_tag, outdir,
                                 fdr_cutoff,           ## From parameters.R (no default)
                                 simplify_go,          ## From parameters.R (no default)
                                 simplify_cutoff,      ## From parameters.R (no default)
                                 top_n_per_direction,  ## From parameters.R (no default)
                                 deg_up_count = NA,    ## Authoritative DEG count for cross-validation
                                 deg_down_count = NA,
                                 logfc_cutoff = NA,
                                 deg_fdr_cutoff = NA) {
    
    cat("\n--- GSEA: ", contrast_name, " ---\n", sep = "")
    cat("[INFO] Ranked gene list size: ", length(ranked_genes), " genes\n", sep = "")
    
    if (length(ranked_genes) < 10) {
      cat("[WARN] Too few genes for GSEA\n")
      return(NULL)
    }
    
    results <- list()
    
    ## GO:BP GSEA using clusterProfiler::gseGO()
    if (run_go_bp) {
      log_time(paste0("Starting GO:BP gseGO for ", contrast_name))
      tryCatch({
        cat("[INFO] Running clusterProfiler::gseGO() (ont = BP)...\n")

        ## Map SYMBOL to ENTREZID for gseGO
        go_symbols <- names(ranked_genes)
        go_symbol_to_entrez <- AnnotationDbi::select(
          org.Mm.eg.db,
          keys = go_symbols,
          columns = c("SYMBOL", "ENTREZID"),
          keytype = "SYMBOL"
        ) %>%
          dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL))

        go_symbol_to_entrez$rank_value <- ranked_genes[go_symbol_to_entrez$SYMBOL]
        go_symbol_to_entrez <- go_symbol_to_entrez %>%
          dplyr::group_by(ENTREZID) %>%
          dplyr::slice_max(order_by = abs(rank_value), n = 1, with_ties = FALSE) %>%
          dplyr::ungroup()

        go_genelist <- setNames(go_symbol_to_entrez$rank_value, go_symbol_to_entrez$ENTREZID)
        go_genelist <- sort(go_genelist, decreasing = TRUE)

        gsea_go <- gseGO(
          geneList     = go_genelist,
          OrgDb        = org.Mm.eg.db,
          ont          = "BP",
          keyType      = "ENTREZID",
          minGSSize    = gsea_params$min_gs_size,
          maxGSSize    = gsea_params$max_gs_size,
          pvalueCutoff = 1.0,
          pAdjustMethod = "BH",
          verbose      = FALSE,
          seed         = TRUE
        )

        if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
          ## SAVE gseaResult object for barcode plots (part5)
          gsea_obj_file <- paste0("gseaResult_gobp_", contrast_name, "_", run_tag, ".rds")
          saveRDS(gsea_go, file = file.path(outdir, "tables", gsea_obj_file))
          cat("[OK] Saved gseaResult object: ", gsea_obj_file, "\n", sep = "")

          go_results <- as.data.frame(gsea_go) %>%
            dplyr::rename(
              pathway = ID,
              padj = p.adjust,
              leadingEdge = core_enrichment
            ) %>%
            dplyr::arrange(padj)

          go_entrez_to_symbol <- go_symbol_to_entrez %>%
            dplyr::distinct(ENTREZID, SYMBOL)
          if ("leadingEdge" %in% colnames(go_results)) {
            go_results$leadingEdge_symbols <- vapply(
              go_results$leadingEdge,
              function(ids) {
                if (is.na(ids) || ids == "") {
                  return(NA_character_)
                }
                id_vec <- unlist(strsplit(ids, "/", fixed = TRUE))
                symbol_vec <- go_entrez_to_symbol$SYMBOL[match(id_vec, go_entrez_to_symbol$ENTREZID)]
                symbol_vec <- symbol_vec[!is.na(symbol_vec)]
                if (length(symbol_vec) == 0) {
                  return(NA_character_)
                }
                paste(symbol_vec, collapse = ",")
              },
              character(1)
            )
          }

          results$gobp <- go_results
          n_sig <- sum(go_results$padj < fdr_cutoff, na.rm = TRUE)
          n_sig_up <- sum(go_results$padj < fdr_cutoff & go_results$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(go_results$padj < fdr_cutoff & go_results$NES < 0, na.rm = TRUE)
          cat("[OK] gseGO GO:BP: ", nrow(go_results), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")

          ## Track in sanity checker (includes DEG counts for cross-validation)
          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "GO:BP",
              N_Input_Genes = length(go_genelist),
              N_DEG_Up = deg_up_count,
              N_DEG_Down = deg_down_count,
              N_Pathways_Tested = nrow(go_results),
              N_Sig_Pathways = n_sig,
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
              logFC_cutoff = logfc_cutoff,
              fdr_cutoff = deg_fdr_cutoff,
              stringsAsFactors = FALSE
            )
          )

          ## Optionally simplify GO:BP results to reduce redundancy
          if (simplify_go && nrow(go_results) > 0) {
            cat("[INFO] Simplifying GO:BP results (removing redundant terms)...\n")
            tryCatch({
              sig_go <- go_results %>% dplyr::filter(padj < fdr_cutoff)
              if (nrow(sig_go) > 0) {
                require(GOSemSim)
                godata <- godata('org.Mm.eg.db', ont = "BP", computeIC = TRUE)

                go_ids <- sig_go$pathway
                sim_mat <- mgoSim(go_ids, go_ids, semData = godata, measure = "Wang")

                keep_terms <- c()
                for (i in seq_along(go_ids)) {
                  if (i == 1 || all(sim_mat[i, keep_terms] < simplify_cutoff)) {
                    keep_terms <- c(keep_terms, i)
                  }
                }

                go_results_simplified <- sig_go[keep_terms, ]
                n_removed <- nrow(sig_go) - nrow(go_results_simplified)
                cat("[OK] Removed ", n_removed, " redundant GO:BP terms\n", sep = "")
                cat("[OK] Simplified GO:BP results: ", nrow(go_results_simplified), " terms\n", sep = "")

                results$gobp_simplified <- go_results_simplified
              }
            }, error = function(e) {
              cat("[WARN] GO:BP simplification failed: ", conditionMessage(e), "\n", sep = "")
            })
          }
        } else {
          cat("[INFO] gseGO returned no results\n")
        }
      }, error = function(e) {
        cat("[WARN] gseGO GO:BP failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "GO:BP",
            N_Input_Genes = length(ranked_genes),
            N_DEG_Up = deg_up_count,
            N_DEG_Down = deg_down_count,
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            logFC_cutoff = logfc_cutoff,
            fdr_cutoff = deg_fdr_cutoff,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed GO:BP gseGO for ", contrast_name))
    }
    
    ## ============================================================
    ## KEGG GSEA using clusterProfiler::gseKEGG()
    ## UPDATED: Replaced fgsea+msigdbr with gseKEGG for consistency with ORA
    ## ============================================================
    if (run_kegg) {
      log_time(paste0("Starting KEGG gseKEGG for ", contrast_name))
      tryCatch({
        cat("[INFO] Running clusterProfiler::gseKEGG() (organism = mmu)...\n")

        ## Step 1: Map SYMBOL to ENTREZID
        ## gseKEGG requires ENTREZID as gene identifiers
        gene_symbols <- names(ranked_genes)

        symbol_to_entrez <- AnnotationDbi::select(
          org.Mm.eg.db,
          keys = gene_symbols,
          columns = c("SYMBOL", "ENTREZID"),
          keytype = "SYMBOL"
        )

        ## Remove NA mappings
        symbol_to_entrez <- symbol_to_entrez %>%
          dplyr::filter(!is.na(ENTREZID), !is.na(SYMBOL))

        ## Create mapping: add ranking values
        symbol_to_entrez$rank_value <- ranked_genes[symbol_to_entrez$SYMBOL]

        ## Handle duplicates: keep gene with largest absolute rank value
        symbol_to_entrez <- symbol_to_entrez %>%
          dplyr::group_by(ENTREZID) %>%
          dplyr::slice_max(order_by = abs(rank_value), n = 1, with_ties = FALSE) %>%
          dplyr::ungroup()

        ## Create named vector: ENTREZID -> rank value (sorted decreasing)
        kegg_genelist <- setNames(symbol_to_entrez$rank_value, symbol_to_entrez$ENTREZID)
        kegg_genelist <- sort(kegg_genelist, decreasing = TRUE)

        cat("[INFO] Mapped ", length(kegg_genelist), " genes to ENTREZID for gseKEGG\n", sep = "")

        ## Step 2: Run gseKEGG
        gsea_kegg <- gseKEGG(
          geneList     = kegg_genelist,
          organism     = gse_kegg_params$organism,
          minGSSize    = gse_kegg_params$min_gs_size,
          maxGSSize    = gse_kegg_params$max_gs_size,
          pvalueCutoff = gse_kegg_params$pvalue_cutoff,        ## Initial filter (will use fdr_cutoff for final)
          pAdjustMethod = gse_kegg_params$p_adjust_method,
          verbose      = FALSE,
          seed         = gse_kegg_params$seed        ## For reproducibility
        )

        if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
          ## SAVE gseaResult object for barcode plots (part5)
          gsea_obj_file <- paste0("gseaResult_kegg_", contrast_name, "_", run_tag, ".rds")
          saveRDS(gsea_kegg, file = file.path(outdir, "tables", gsea_obj_file))
          cat("[OK] Saved gseaResult object: ", gsea_obj_file, "\n", sep = "")

          ## Convert to data frame
          kegg_results <- as.data.frame(gsea_kegg)

          ## Rename columns to match fgsea output format
          kegg_results <- kegg_results %>%
            dplyr::rename(
              pathway = ID,
              padj = p.adjust,
              leadingEdge = core_enrichment
            )

          ## Remove duplicate/legacy columns if present
          if ("leading_edge" %in% colnames(kegg_results)) {
            kegg_results <- kegg_results %>% dplyr::select(-leading_edge)
          }

          ## Map KEGG leading edge Entrez IDs back to gene symbols
          entrez_to_symbol <- symbol_to_entrez %>%
            dplyr::distinct(ENTREZID, SYMBOL)
          if ("leadingEdge" %in% colnames(kegg_results)) {
            kegg_results$leadingEdge_symbols <- vapply(
              kegg_results$leadingEdge,
              function(ids) {
                if (is.na(ids) || ids == "") {
                  return(NA_character_)
                }
                id_vec <- unlist(strsplit(ids, "/", fixed = TRUE))
                symbol_vec <- entrez_to_symbol$SYMBOL[match(id_vec, entrez_to_symbol$ENTREZID)]
                symbol_vec <- symbol_vec[!is.na(symbol_vec)]
                if (length(symbol_vec) == 0) {
                  return(NA_character_)
                }
                paste(symbol_vec, collapse = ",")
              },
              character(1)
            )
          }

          kegg_results <- kegg_results %>%
            dplyr::arrange(padj)

          results$kegg <- kegg_results
          n_sig <- sum(kegg_results$padj < fdr_cutoff, na.rm = TRUE)
          n_sig_up <- sum(kegg_results$padj < fdr_cutoff & kegg_results$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(kegg_results$padj < fdr_cutoff & kegg_results$NES < 0, na.rm = TRUE)
          cat("[OK] gseKEGG: ", nrow(kegg_results), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")

          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "KEGG",
              N_Input_Genes = length(kegg_genelist),
              N_DEG_Up = deg_up_count,
              N_DEG_Down = deg_down_count,
              N_Pathways_Tested = nrow(kegg_results),
              N_Sig_Pathways = n_sig,
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
              logFC_cutoff = logfc_cutoff,
              fdr_cutoff = deg_fdr_cutoff,
              stringsAsFactors = FALSE
            )
          )
        } else {
          cat("[INFO] gseKEGG returned no results\n")
          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "KEGG",
              N_Input_Genes = length(kegg_genelist),
              N_DEG_Up = deg_up_count,
              N_DEG_Down = deg_down_count,
              N_Pathways_Tested = 0,
              N_Sig_Pathways = 0,
              N_Sig_Up = 0,
              N_Sig_Down = 0,
              logFC_cutoff = logfc_cutoff,
              fdr_cutoff = deg_fdr_cutoff,
              stringsAsFactors = FALSE
            )
          )
        }
      }, error = function(e) {
        cat("[WARN] gseKEGG failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "KEGG",
            N_Input_Genes = length(ranked_genes),
            N_DEG_Up = deg_up_count,
            N_DEG_Down = deg_down_count,
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            logFC_cutoff = logfc_cutoff,
            fdr_cutoff = deg_fdr_cutoff,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed KEGG gseKEGG for ", contrast_name))
    }
    
    ## WikiPathways and Hallmark analyses are disabled in this workflow.
    
    ## Save results with direction-specific colors
    if (length(results) > 0) {
      for (db_name in names(results)) {
        fgsea_obj <- results[[db_name]]
        
        sig_results <- fgsea_obj %>%
          dplyr::filter(padj < fdr_cutoff) %>%
          dplyr::arrange(padj)
        
        if (nrow(sig_results) == 0) {
          cat("[INFO] No significant GSEA pathways in ", db_name, "\n", sep = "")
          next
        }
        ## Ensure pathway_name exists (for plotting)
        if (!"Description" %in% colnames(sig_results) && !"pathway_name" %in% colnames(sig_results)) {
          sig_results$Description <- sig_results$pathway
        }

        ## Convert list columns to character (leadingEdge is a list)
        for (col_name in colnames(sig_results)) {
          if (is.list(sig_results[[col_name]])) {
            sig_results[[col_name]] <- sapply(sig_results[[col_name]], function(x) {
              if (is.null(x) || length(x) == 0) {
                return(NA_character_)
              } else {
                return(paste(x, collapse = ","))
              }
            })
          }
        }

        ## Standardize column names for consistent output
        if ("pathway" %in% colnames(sig_results)) {
          sig_results <- sig_results %>% dplyr::rename(pathway_id = pathway)
        }
        if ("Description" %in% colnames(sig_results)) {
          sig_results <- sig_results %>% dplyr::rename(pathway_name = Description)
        }

        ## Remove qvalue column if present (keep padj only for clarity)
        if ("qvalue" %in% colnames(sig_results) && "padj" %in% colnames(sig_results)) {
          sig_results <- sig_results %>% dplyr::select(-qvalue)
        }
        if ("qvalues" %in% colnames(sig_results)) {
          sig_results <- sig_results %>% dplyr::select(-qvalues)
        }

        ## Add metadata columns for downstream analysis
        sig_results$contrast <- contrast_name
        sig_results$tissue <- ifelse(grepl("^iWAT", contrast_name), "iWAT",
                                     ifelse(grepl("^gWAT", contrast_name), "gWAT", "unknown"))
        sig_results$database <- ifelse(db_name == "gobp", "GO:BP",
                                       ifelse(db_name == "kegg", "KEGG", toupper(db_name)))
        sig_results$ranking_metric <- "DESeq2_Wald_stat"
        sig_results$analysis_type <- "fgsea"
        if (!is.na(logfc_cutoff)) sig_results$logFC_cutoff <- logfc_cutoff
        if (!is.na(deg_fdr_cutoff)) sig_results$fdr_cutoff <- deg_fdr_cutoff

        ## Save table
        table_file <- paste0("fgsea_", db_name, "_", contrast_name, "_", run_tag, ".csv")
        write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
        cat("[OK] Saved fgsea table: ", table_file, "\n", sep = "")

        if (db_name == "gobp_simplified") {
          next
        }
        
        ## Create combined plot (balanced up/down)
        top_up <- sig_results %>%
          dplyr::filter(NES > 0) %>%
          dplyr::arrange(padj, dplyr::desc(NES)) %>%
          dplyr::slice_head(n = top_n_per_direction)
        
        top_down <- sig_results %>%
          dplyr::filter(NES < 0) %>%
          dplyr::arrange(padj, NES) %>%
          dplyr::slice_head(n = top_n_per_direction)
        
        plot_data <- dplyr::bind_rows(top_down, top_up) %>%
          dplyr::arrange(NES) %>%
          dplyr::mutate(
            pathway_label = ifelse(!is.na(pathway_name) & pathway_name != "", pathway_name, pathway_id),
            pathway_label = factor(pathway_label, levels = pathway_label),
            Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")
          )
        
        if (nrow(plot_data) > 0) {
          ## Parameter caption
          if (db_name == "kegg") {
            param_caption <- paste(
              paste0(
                "FDR < ", fdr_cutoff,
                " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
              ),
              paste0("Method: gseKEGG (", gse_kegg_params$organism, ") | Color: adj. p-value"),
              sep = "\n"
            )
          } else if (db_name == "gobp") {
            simplify_label <- if (simplify_go) paste0("Simplified (cutoff=", simplify_cutoff, ")") else "Not simplified"
            param_caption <- paste(
              paste0(
                "FDR < ", fdr_cutoff,
                " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
              ),
              paste0("Method: gseGO (BP) | ", simplify_label, " | Color: adj. p-value"),
              sep = "\n"
            )
          } else {
            param_caption <- paste(
              paste0(
                "FDR < ", fdr_cutoff,
                " | Top ", top_n_per_direction, " per direction | Ranked by ", gsea_params$rank_metric
              ),
              "Method: GSEA | Color: adj. p-value",
              sep = "\n"
            )
          }

          ## Color by adjusted p-value with direction-specific gradients
          ## Orange gradient for up-regulated, blue gradient for down-regulated
          ## This matches the ORA bar plot style
          method_label <- "GSEA"

          ## Separate up and down pathways for color mapping
          ## Lower padj = more significant = darker color
          plot_data <- plot_data %>%
            dplyr::mutate(
              ## Use -log10(padj) for color intensity
              ## Higher value = more significant
              padj_safe = pmax(padj, .Machine$double.eps),
              neg_log10_padj = -log10(padj_safe)
            )

          ## Calculate nice breaks for the legend (show actual p-values)
          min_padj <- min(plot_data$padj, na.rm = TRUE)
          max_padj <- max(plot_data$padj, na.rm = TRUE)

          ## Create a custom fill aesthetic per row based on direction
          plot_data <- plot_data %>%
            dplyr::mutate(
              fill_value = ifelse(Direction == "Up-regulated", neg_log10_padj, -neg_log10_padj)
            )

          ## Get range for symmetric scale
          max_neg_log <- max(abs(plot_data$fill_value), na.rm = TRUE)

          p_fgsea <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = fill_value)) +
            geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
            scale_fill_gradientn(
              colors = c(
                "#0072B2",  ## Dark blue (most significant down)
                "#9ECAE1",  ## Light blue (least significant down)
                "grey95",   ## Near white (threshold)
                "#FDD49E",  ## Light orange (least significant up)
                "#E69F00"   ## Dark orange (most significant up)
              ),
              values = scales::rescale(c(-max_neg_log, -1.3, 0, 1.3, max_neg_log)),
              limits = c(-max_neg_log, max_neg_log),
              name = "Adj. P-value",
              breaks = c(-max_neg_log, -1.3, 0, 1.3, max_neg_log),
              labels = function(x) {
                sapply(x, function(val) {
                  if (abs(val) < 1e-6) return("1")
                  pval <- 10^(-abs(val))
                  format(pval, scientific = TRUE, digits = 1)
                })
              },
              guide = guide_colorbar(
                title = "Adj. P-value\n(blue=down, orange=up)",
                title.position = "top",
                barwidth = 1,
                barheight = 4
              )
            ) +
            labs(
              title = paste0(method_label, " ", toupper(db_name), ": ", contrast_name),
              x = "Normalized Enrichment Score (NES)",
              y = NULL,
              caption = param_caption
            ) +
            theme_classic(base_size = 12) +
            theme(
              plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
              axis.text.y = element_text(size = 10, color = "black"),
              axis.text.x = element_text(size = 10, color = "black", face = "bold"),
              axis.title.x = element_text(size = 12, face = "bold"),
              legend.title = element_text(size = 9, face = "bold"),
              legend.text = element_text(size = 8),
              legend.position = "right",
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8))
            ) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
            ## Add annotation for direction interpretation
            annotate("text", x = Inf, y = Inf, label = "Up", hjust = 1.1, vjust = -0.5,
                     size = 3, color = "#E69F00", fontface = "bold") +
            annotate("text", x = -Inf, y = Inf, label = "Down", hjust = -0.1, vjust = -0.5,
                     size = 3, color = "#0072B2", fontface = "bold")
          
          plot_file <- paste0("fgsea_plot_", db_name, "_", contrast_name, "_", run_tag, ".png")
          ggsave(file.path(outdir, "plots", plot_file), plot = p_fgsea, width = 12,
                 height = max(6, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
          cat("[OK] Saved fgsea plot: ", plot_file, "\n", sep = "")
        }
      }
    }
    
    return(results)
  }
  
  ## 4.4) Process each DE table ------------------------------
  cat("\n=== RUNNING GSEA FOR ALL CONTRASTS ===\n")
  
  for (de_file in de_files) {
    
    ## Extract contrast name from filename
    de_filename <- basename(de_file)
    contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)
    
    cat("\n=== Processing contrast: ", contrast_name, " ===\n", sep = "")
    log_time(paste0("Starting GSEA for contrast: ", contrast_name))
    
    ## Load DE table
    de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)
    
    ## Ensure we have the stat column (Wald statistic)
    if (!"stat" %in% colnames(de_table)) {
      cat("[WARN] Cannot find 'stat' column (Wald statistic) in ", de_file, "\n")
      next
    }
    
    ## Create ranked gene list using Wald statistic
    ranked_genes <- de_table %>%
      dplyr::filter(is.finite(stat)) %>%
      dplyr::arrange(dplyr::desc(stat))
    
    ## Create named vector (gene names -> Wald statistic)
    ranked_vec <- setNames(ranked_genes$stat, rownames(ranked_genes))
    n_pos <- sum(ranked_vec > 0, na.rm = TRUE)
    n_neg <- sum(ranked_vec < 0, na.rm = TRUE)
    n_zero <- sum(ranked_vec == 0, na.rm = TRUE)
    n_total <- length(ranked_vec)
    pct_pos <- if (n_total == 0) 0 else round(100 * n_pos / n_total, 1)
    pct_neg <- if (n_total == 0) 0 else round(100 * n_neg / n_total, 1)
    pct_zero <- if (n_total == 0) 0 else round(100 * n_zero / n_total, 1)
    cat("[INFO] Ranked stats distribution: +", n_pos, " (", pct_pos, "%), -",
        n_neg, " (", pct_neg, "%), 0=", n_zero, " (", pct_zero, "%)\n", sep = "")
    if (pct_neg < 5) {
      cat("[WARN] Low fraction of negative stats; down-regulated pathways may be scarce.\n")
    }
    
    cat("[INFO] Ranked gene list for GSEA: ", length(ranked_vec), " genes\n", sep = "")
    
    ## Get authoritative DEG counts for this contrast (for cross-validation)
    deg_up_count <- NA
    deg_down_count <- NA
    logfc_cutoff <- NA
    deg_fdr_cutoff <- NA
    if (!is.null(authoritative_deg_lists) && contrast_name %in% names(authoritative_deg_lists)) {
      deg_up_count <- length(authoritative_deg_lists[[contrast_name]]$up)
      deg_down_count <- length(authoritative_deg_lists[[contrast_name]]$down)
      logfc_cutoff <- authoritative_deg_lists[[contrast_name]]$logFC_cutoff
      deg_fdr_cutoff <- authoritative_deg_lists[[contrast_name]]$fdr_cutoff
      cat("[INFO] Authoritative DEG counts for cross-validation: Up=", deg_up_count,
          ", Down=", deg_down_count, "\n", sep = "")
    }

    ## Run GSEA (gseGO for GO:BP, gseKEGG for KEGG)
    if (length(ranked_vec) >= 10) {
      run_fgsea_analysis(
        ranked_genes  = ranked_vec,
        contrast_name = contrast_name,
        go_bp_list    = go_bp_list,
        go_bp_term2gene = go_bp_term2gene,
        wikipathways_list = wikipathways_list,
        hallmark_list = hallmark_list,
        de_table      = de_table,    ## Added for ENTREZID mapping in gseKEGG
        run_tag       = run_tag,
        outdir        = outdir,
        fdr_cutoff    = fdr_cut,
        simplify_go = simplify_go_bp,
        simplify_cutoff = simplify_go_cutoff,
        top_n_per_direction = top_n_per_direction,
        deg_up_count = deg_up_count,
        deg_down_count = deg_down_count,
        logfc_cutoff = logfc_cutoff,
        deg_fdr_cutoff = deg_fdr_cutoff
      )
      log_time(paste0("Completed GSEA for contrast: ", contrast_name))
    } else {
      cat("[WARN] Too few genes for GSEA\n")
    }
  }
  
  ## 4.5) Save sanity check table ----------------------------
  sanity_file <- paste0("fgsea_sanity_check_", run_tag, ".csv")
  write.csv(fgsea_sanity_tracker, file = file.path(outdir, "tables", sanity_file), row.names = FALSE)
  cat("\n[OK] Saved fgsea sanity check table: ", sanity_file, "\n", sep = "")
  
  ## 4.6) Print comprehensive sanity check summary -----------
  cat("\n==========================================================\n")
  cat("=== SANITY CHECK SUMMARY: GSEA PATHWAY ANALYSIS ===\n")
  cat("==========================================================\n")
  cat("Methods:\n")
  cat("  - GO:BP: clusterProfiler::gseGO() (rank-based GSEA)\n")
  cat("  - KEGG: clusterProfiler::gseKEGG() (matches ORA's enrichKEGG)\n")
  cat("  - WikiPathways/Hallmark: disabled in this run\n")
  cat("Note: GSEA uses ranked gene list as universe (no separate background needed)\n")
  cat("Ranking metric: DESeq2 Wald statistic\n\n")
  
  ## Summarize across all contrasts
  for (cn_name in unique(fgsea_sanity_tracker$Contrast)) {
    cat("--- ", cn_name, " ---\n", sep = "")
    subset_data <- fgsea_sanity_tracker[fgsea_sanity_tracker$Contrast == cn_name, ]
    
    for (i in seq_len(nrow(subset_data))) {
      row <- subset_data[i, ]
      cat("  ", row$Database, ":\n", sep = "")
      cat("    Ranked gene list size: ", row$N_Input_Genes, "\n", sep = "")
      cat("    Pathways tested: ", row$N_Pathways_Tested, "\n", sep = "")
      cat("    Significant pathways (padj<", fdr_cut, "): ", row$N_Sig_Pathways, "\n", sep = "")
      if (row$N_Sig_Pathways > 0) {
        cat("      Up-regulated (NES>0): ", row$N_Sig_Up, "\n", sep = "")
        cat("      Down-regulated (NES<0): ", row$N_Sig_Down, "\n", sep = "")
      }
    }
    cat("\n")
  }
  cat("==========================================================\n\n")

  ## Cross-validate with DEG_summary and ORA sanity check
  cat("\n==========================================================\n")
  cat("=== CROSS-VALIDATION: fgsea vs DEG_summary vs ORA ===\n")
  cat("==========================================================\n")
  cat("Note: GSEA uses ALL ranked genes for analysis.\n")
  cat("      N_DEG_Up/N_DEG_Down are recorded for cross-validation only.\n\n")

  all_pass <- TRUE

  ## Load DEG_summary
  deg_summary_file <- list.files(file.path(outdir, "tables"), pattern = "^DEG_summary_tissue_.*\\.csv$", full.names = TRUE)
  if (length(deg_summary_file) > 0) {
    deg_summary_file <- deg_summary_file[order(file.info(deg_summary_file)$mtime, decreasing = TRUE)][1]
    deg_summary <- read.csv(deg_summary_file, stringsAsFactors = FALSE)
    cat("[INFO] Loaded DEG_summary from: ", basename(deg_summary_file), "\n", sep = "")
  } else {
    deg_summary <- NULL
    cat("[WARN] DEG_summary not found\n")
  }

  ## Load ORA sanity check
  ora_sanity_file <- list.files(file.path(outdir, "tables"), pattern = "^ORA_sanity_check_.*\\.csv$", full.names = TRUE)
  if (length(ora_sanity_file) > 0) {
    ora_sanity_file <- ora_sanity_file[order(file.info(ora_sanity_file)$mtime, decreasing = TRUE)][1]
    ora_sanity <- read.csv(ora_sanity_file, stringsAsFactors = FALSE)
    cat("[INFO] Loaded ORA sanity check from: ", basename(ora_sanity_file), "\n\n", sep = "")
  } else {
    ora_sanity <- NULL
    cat("[WARN] ORA sanity check not found\n\n")
  }

  for (cn in unique(fgsea_sanity_tracker$Contrast)) {
    cat("--- ", cn, " ---\n", sep = "")

    ## Get fgsea DEG counts (from GO:BP row)
    fgsea_row <- fgsea_sanity_tracker %>%
      dplyr::filter(Contrast == cn, Database == "GO:BP")
    fgsea_up <- if (nrow(fgsea_row) > 0) fgsea_row$N_DEG_Up[1] else NA
    fgsea_down <- if (nrow(fgsea_row) > 0) fgsea_row$N_DEG_Down[1] else NA

    ## Get DEG_summary counts
    deg_up <- NA
    deg_down <- NA
    if (!is.null(deg_summary)) {
      deg_row <- deg_summary %>% dplyr::filter(Contrast == cn)
      if (nrow(deg_row) > 0) {
        deg_up <- deg_row$N_Up[1]
        deg_down <- deg_row$N_Down[1]
      }
    }

    ## Get ORA counts
    ora_up <- NA
    ora_down <- NA
    if (!is.null(ora_sanity)) {
      ora_up_row <- ora_sanity %>% dplyr::filter(Contrast == cn, Direction == "Up", Database == "GO:BP")
      ora_down_row <- ora_sanity %>% dplyr::filter(Contrast == cn, Direction == "Down", Database == "GO:BP")
      if (nrow(ora_up_row) > 0) ora_up <- ora_up_row$N_Input_Genes[1]
      if (nrow(ora_down_row) > 0) ora_down <- ora_down_row$N_Input_Genes[1]
    }

    ## Check consistency
    up_match <- "PASS"
    down_match <- "PASS"
    if (!is.na(fgsea_up) && !is.na(deg_up) && fgsea_up != deg_up) {
      up_match <- "FAIL"
      all_pass <- FALSE
    }
    if (!is.na(fgsea_up) && !is.na(ora_up) && fgsea_up != ora_up) {
      up_match <- "FAIL"
      all_pass <- FALSE
    }
    if (!is.na(fgsea_down) && !is.na(deg_down) && fgsea_down != deg_down) {
      down_match <- "FAIL"
      all_pass <- FALSE
    }
    if (!is.na(fgsea_down) && !is.na(ora_down) && fgsea_down != ora_down) {
      down_match <- "FAIL"
      all_pass <- FALSE
    }

    cat("  Up:   fgsea=", fgsea_up, " | DEG_summary=", deg_up, " | ORA=", ora_up, " [", up_match, "]\n", sep = "")
    cat("  Down: fgsea=", fgsea_down, " | DEG_summary=", deg_down, " | ORA=", ora_down, " [", down_match, "]\n", sep = "")
  }

  cat("\n")
  if (all_pass) {
    cat("VALIDATION RESULT: PASS - DEG counts match across all scripts\n")
    cat("  -> No drift in DEG membership between Part 1, Part 2, and Part 3\n")
    cat("  -> ORA input genes = fgsea sanity check genes = DEG_summary\n")
  } else {
    cat("VALIDATION RESULT: FAIL - Mismatch detected!\n")
    cat("  -> Re-run Part 1 to regenerate authoritative DEG lists\n")
    cat("  -> Then re-run Part 2 and Part 3\n")
  }
  cat("==========================================================\n\n")

  cat("\n======================================\n")
  cat("=== PART 3 (GSEA) COMPLETE ===\n")
  cat("======================================\n")
  cat("Methods used:\n")
  cat("  - GO:BP: clusterProfiler::gseGO()\n")
  cat("  - KEGG: clusterProfiler::gseKEGG()\n")
  cat("  - WikiPathways/Hallmark: disabled\n")
  cat("======================================\n\n")

  ## 4.7) Save session info for reproducibility -------------
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_fgsea_", run_tag, ".txt"))
  sink(session_file)
  cat("=== R SESSION INFORMATION ===\n\n")
  print(sessionInfo())
  cat("\n=== KEY PACKAGE VERSIONS ===\n\n")
  ## clusterProfiler used for: gseKEGG (KEGG GSEA)
  key_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "dplyr", "ggplot2", "GOSemSim")
  print(installed.packages()[intersect(key_pkgs, rownames(installed.packages())), c("Version", "Built")])
  sink()
  cat("[INFO] Session info saved to: ", session_file, "\n", sep = "")

}, error = function(e) {
  
  cat("\n======================================\n")
  cat("=== ERROR IN GSEA PIPELINE ===\n")
  cat("======================================\n")
  cat(conditionMessage(e), "\n\n")
  
  err_file <- paste0("ERROR_gsea_", run_tag, ".txt")
  writeLines(
    c(
      "gsea script error:",
      conditionMessage(e),
      "",
      "Traceback:",
      capture.output(traceback())
    ),
    con = file.path(outdir, "logs", err_file)
  )
  
  stop(e)
})

cat("\n=== Part 3 (GSEA) finished ===\n")
