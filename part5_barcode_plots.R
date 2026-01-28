#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 5: BARCODE PLOTS FOR TOP fGSEA PATHWAYS
## =========================================================
## Creates barcode (enrichment) plots for ALL top pathways
## that were plotted in the fGSEA bar graphs (part3)
##
## Uses enrichplot::gseaplot2() with the saved gseaResult objects
## Run AFTER parts 1-3 have completed
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "stringr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "clusterProfiler",
  "enrichplot",
  "DOSE"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
})

## 2.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

## 3) Find most recent run directory -----------------------
cat("\n=== FINDING ANALYSIS RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run parts 1-3 first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No run directories found in savepoints/")
}

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")

## Create barcode plots subdirectory
barcode_dir <- file.path(plots_dir, "barcode_plots")
if (!dir.exists(barcode_dir)) {
  dir.create(barcode_dir, recursive = TRUE)
}

## ============================================================
## PARAMETERS (from central parameters.R)
## ============================================================

## Use values from barcode_params (which references gsea_params for consistency)
top_n_per_direction <- barcode_params$top_n_per_direction
fdr_cutoff <- barcode_params$fdr_cutoff

cat("[INFO] Using parameters: top_n_per_direction=", top_n_per_direction,
    ", fdr_cutoff=", fdr_cutoff, "\n", sep = "")

## ============================================================
## SECTION 1: FIND AND LOAD gseaResult OBJECTS
## ============================================================

cat("\n=== LOADING gseaResult OBJECTS ===\n")

## Find gseaResult RDS files
gsea_rds_files <- list.files(tables_dir, pattern = "^gseaResult_(gobp|kegg)_.*\\.rds$", full.names = TRUE)

if (length(gsea_rds_files) == 0) {
  stop("No gseaResult RDS files found in: ", tables_dir, "\n",
       "Please re-run part3_fgsea.R to generate gseaResult objects.")
}

cat("[OK] Found ", length(gsea_rds_files), " gseaResult files\n", sep = "")

## ============================================================
## SECTION 2: CREATE BARCODE PLOTS FOR ALL TOP PATHWAYS
## ============================================================

cat("\n=== CREATING BARCODE PLOTS ===\n")

## Summary tracker
barcode_summary <- data.frame(
  Database = character(),
  Contrast = character(),
  Pathway_ID = character(),
  Pathway_Name = character(),
  NES = numeric(),
  padj = numeric(),
  Direction = character(),
  stringsAsFactors = FALSE
)

## Process each gseaResult file
for (rds_file in gsea_rds_files) {

  filename <- basename(rds_file)

  ## Parse filename: gseaResult_gobp_ContrastName_YYYYMMDD_HHMMSS.rds
  parts <- str_match(filename, "gseaResult_(gobp|kegg)_(.+)_[0-9]{8}_[0-9]{6}\\.rds")

  if (is.na(parts[1])) {
    cat("[WARN] Cannot parse filename: ", filename, "\n")
    next
  }

  db_name <- parts[2]
  contrast_name <- parts[3]
  database <- ifelse(db_name == "gobp", "GO:BP", "KEGG")

  cat("\n--- Processing: ", database, " | ", contrast_name, " ---\n", sep = "")

  ## Load gseaResult object
  gsea_result <- tryCatch({
    readRDS(rds_file)
  }, error = function(e) {
    cat("[WARN] Cannot load file: ", filename, "\n")
    return(NULL)
  })

  if (is.null(gsea_result)) {
    next
  }

  ## Diagnostic: Check gseaResult object structure
  cat("[DEBUG] Object class: ", class(gsea_result)[1], "\n", sep = "")
  cat("[DEBUG] Has @geneList: ", !is.null(gsea_result@geneList),
      " (length=", length(gsea_result@geneList), ")\n", sep = "")
  cat("[DEBUG] Has @geneSets: ", !is.null(gsea_result@geneSets),
      " (count=", length(gsea_result@geneSets), ")\n", sep = "")

  ## Check if geneList is valid (numeric and named)
  if (length(gsea_result@geneList) > 0) {
    gl <- gsea_result@geneList
    cat("[DEBUG] geneList type: ", typeof(gl), ", is.numeric: ", is.numeric(gl), "\n", sep = "")
    if (!is.numeric(gl)) {
      cat("[WARN] geneList is not numeric - this will cause plotting errors\n")
    }
  }

  ## Get results data frame
  result_df <- gsea_result@result

  if (nrow(result_df) == 0) {
    cat("[WARN] No pathways in gseaResult\n")
    next
  }

  ## Save full gseaResult results table for this contrast/database
  results_file <- paste0("gseaResult_full_", db_name, "_", contrast_name, "_", run_tag, ".csv")
  write.csv(result_df, file = file.path(tables_dir, results_file), row.names = FALSE)
  cat("[OK] Saved full gseaResult table: ", results_file, "\n", sep = "")

  ## Save ranked gene list used for barcode plots
  gene_list <- gsea_result@geneList
  if (length(gene_list) > 0) {
    gene_list_df <- data.frame(
      gene_id = names(gene_list),
      ranking_metric = as.numeric(gene_list),
      stringsAsFactors = FALSE
    )
    gene_list_file <- paste0("gseaResult_geneList_", db_name, "_", contrast_name, "_", run_tag, ".csv")
    write.csv(gene_list_df, file = file.path(tables_dir, gene_list_file), row.names = FALSE)
    cat("[OK] Saved gseaResult geneList: ", gene_list_file, "\n", sep = "")
  } else {
    cat("[WARN] gseaResult geneList is empty; skipping geneList export.\n")
  }

  ## Filter significant pathways
  sig_results <- result_df %>%
    dplyr::filter(p.adjust < fdr_cutoff)

  if (nrow(sig_results) == 0) {
    cat("[INFO] No significant pathways (FDR < ", fdr_cutoff, ")\n", sep = "")
    next
  }

  cat("[INFO] ", nrow(sig_results), " significant pathways\n", sep = "")

  ## Get top pathways per direction (same as bar plots in part3)
  top_up <- sig_results %>%
    dplyr::filter(NES > 0) %>%
    dplyr::arrange(p.adjust, dplyr::desc(NES)) %>%
    dplyr::slice_head(n = top_n_per_direction)

  top_down <- sig_results %>%
    dplyr::filter(NES < 0) %>%
    dplyr::arrange(p.adjust, NES) %>%
    dplyr::slice_head(n = top_n_per_direction)

  ## Combine top pathways
  top_pathways <- dplyr::bind_rows(top_up, top_down)

  if (nrow(top_pathways) == 0) {
    cat("[INFO] No top pathways to plot\n")
    next
  }

  cat("[INFO] Creating barcode plots for ", nrow(top_pathways), " pathways ",
      "(", nrow(top_up), " up, ", nrow(top_down), " down)\n", sep = "")

  ## Create barcode plot for each top pathway
  for (i in seq_len(nrow(top_pathways))) {
    pathway <- top_pathways[i, ]
    pathway_id <- pathway$ID
    pathway_name <- pathway$Description
    nes <- pathway$NES
    padj <- pathway$p.adjust
    direction <- ifelse(nes > 0, "Up", "Down")
    stats_label <- paste0("NES=", round(nes, 2), " | FDR=", signif(padj, 3))

    ## Create clean filename
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", pathway_name)
    safe_name <- gsub("_+", "_", safe_name)
    safe_name <- substr(safe_name, 1, 40)

    plot_file <- paste0("Barcode_", db_name, "_", contrast_name, "_", safe_name, "_", run_tag, ".png")

    ## Use enrichplot::gseaplot2 for the barcode plot
    tryCatch({
      ## Try gseaplot2 first - returns a patchwork object
      ## Note: Some versions have issues with saved/loaded objects
      p <- tryCatch({
        gseaplot2(
          gsea_result,
          geneSetID = pathway_id,
          title = paste0(str_wrap(pathway_name, width = 50), "\n",
                         contrast_name, " | ", database, " | ", direction, "-regulated", "\n",
                         stats_label),
          base_size = 10,
          pvalue_table = FALSE
        )
      }, error = function(e1) {
        ## Fallback to simpler gseaplot (single panel)
        tryCatch({
          enrichplot::gseaplot(
            gsea_result,
            geneSetID = pathway_id,
            title = paste0(pathway_name, "\n", contrast_name, " | ", direction, "\n", stats_label)
          )
        }, error = function(e2) {
          ## If both fail, return NULL
          NULL
        })
      })

      if (is.null(p)) {
        stop("Both gseaplot2 and gseaplot failed")
      }

      ggsave(
        file.path(barcode_dir, plot_file),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300,
        bg = "white"
      )

      cat("[OK] ", direction, ": ", substr(pathway_name, 1, 50), "...\n", sep = "")

      ## Add to summary
      barcode_summary <- rbind(
        barcode_summary,
        data.frame(
          Database = database,
          Contrast = contrast_name,
          Pathway_ID = pathway_id,
          Pathway_Name = pathway_name,
          NES = nes,
          padj = padj,
          Direction = direction,
          stringsAsFactors = FALSE
        )
      )

    }, error = function(e) {
      cat("[WARN] Failed to create plot for: ", pathway_name, "\n")
      cat("       Error: ", conditionMessage(e), "\n")
    })
  }
}

## ============================================================
## SECTION 3: SAVE SUMMARY AND FINISH
## ============================================================

cat("\n=== SAVING SUMMARY ===\n")

## Save summary table
summary_file <- paste0("barcode_plots_summary_", run_tag, ".csv")
write.csv(barcode_summary, file = file.path(tables_dir, summary_file), row.names = FALSE)
cat("[OK] Saved summary table: ", summary_file, "\n", sep = "")

## Print summary
cat("\n==========================================================\n")
cat("=== BARCODE PLOTS SUMMARY ===\n")
cat("==========================================================\n")
cat("Total plots created: ", nrow(barcode_summary), "\n\n", sep = "")

if (nrow(barcode_summary) > 0) {
  ## Summary by database and contrast
  summary_table <- barcode_summary %>%
    dplyr::group_by(Database, Contrast, Direction) %>%
    dplyr::summarise(N_Pathways = n(), .groups = "drop")

  for (i in seq_len(nrow(summary_table))) {
    row <- summary_table[i, ]
    cat("  ", row$Contrast, " (", row$Database, ", ", row$Direction, "): ",
        row$N_Pathways, " pathways\n", sep = "")
  }

  cat("\n--- Top Pathways by NES ---\n")

  ## Top 5 up and down overall
  top_up_all <- barcode_summary %>%
    dplyr::filter(Direction == "Up") %>%
    dplyr::arrange(desc(NES)) %>%
    dplyr::slice_head(n = 5)

  top_down_all <- barcode_summary %>%
    dplyr::filter(Direction == "Down") %>%
    dplyr::arrange(NES) %>%
    dplyr::slice_head(n = 5)

  cat("\nTop 5 UP-regulated pathways:\n")
  for (i in seq_len(nrow(top_up_all))) {
    row <- top_up_all[i, ]
    cat("  ", round(row$NES, 2), " | ", row$Contrast, " | ", substr(row$Pathway_Name, 1, 40), "\n", sep = "")
  }

  cat("\nTop 5 DOWN-regulated pathways:\n")
  for (i in seq_len(nrow(top_down_all))) {
    row <- top_down_all[i, ]
    cat("  ", round(row$NES, 2), " | ", row$Contrast, " | ", substr(row$Pathway_Name, 1, 40), "\n", sep = "")
  }

} else {
  cat("[WARN] No barcode plots were created.\n")
  cat("       Check that gseaResult RDS files exist in: ", tables_dir, "\n")
}
cat("\n==========================================================\n")

cat("\n======================================\n")
cat("=== PART 5 (BARCODE PLOTS) COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - Barcode plots for ALL top fGSEA pathways (", nrow(barcode_summary), " total)\n", sep = "")
cat("  - Summary table: ", summary_file, "\n", sep = "")
cat("\nAll plots saved to: ", barcode_dir, "\n\n", sep = "")
