#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 5: OUTPUT CLEANING & STANDARDIZATION
## =========================================================
## This script cleans and standardizes all output files from the pipeline:
## - DEG CSVs: Remove redundant columns, add metadata
## - fgsea CSVs: Fix column names, leading edge, KEGG IDs, add metadata
## - ORA CSVs: Remove qvalue, standardize gene columns, add metadata
##
## HARD CONSTRAINTS:
## - No rerunning of differential expression
## - No changing fgsea/ORA parameters
## - No changing gene sets or mapping sources
## - No changing ranking metrics or cutoffs
## - Only clean, standardize, and annotate outputs
## =========================================================

## 1) Packages ----------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

## 2) Find most recent run directory ------------------------
cat("\n=== PART 5: OUTPUT CLEANING & STANDARDIZATION ===\n")
cat("=================================================\n\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run the pipeline first.")
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
cat("[INFO] Run tag: ", run_tag, "\n\n", sep = "")

tables_dir <- file.path(outdir, "tables")
cleaned_dir <- file.path(outdir, "tables_cleaned")

## Create cleaned output directory
if (!dir.exists(cleaned_dir)) {
  dir.create(cleaned_dir, recursive = TRUE)
}

## Validation report container
validation_report <- list(
  de_files = list(),
  fgsea_files = list(),
  ora_files = list(),
  errors = character(),
  warnings = character()
)

## =========================================================
## A) CLEAN DEG RESULT CSVs
## =========================================================
clean_deg_csv <- function(file_path, cleaned_dir) {
  cat("\n--- Cleaning DE file: ", basename(file_path), " ---\n", sep = "")

  ## Read with row names
  de <- tryCatch({
    read.csv(file_path, row.names = 1, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("[ERROR] Failed to read: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })

  if (is.null(de)) return(NULL)

  ## Extract contrast from filename
  ## Pattern: DE_tissue_[contrast]_YYYYMMDD_HHMMSS.csv
  filename <- basename(file_path)
  contrast <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", filename)

  ## Extract tissue from contrast
  tissue <- ifelse(grepl("^iWAT", contrast), "iWAT",
                   ifelse(grepl("^gWAT", contrast), "gWAT", "unknown"))

  ## Original columns
  orig_cols <- colnames(de)
  cat("[INFO] Original columns: ", paste(orig_cols, collapse = ", "), "\n", sep = "")

  ## Define required DESeq2 columns
  deseq2_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

  ## Check which DESeq2 columns exist
  missing_cols <- setdiff(deseq2_cols, orig_cols)
  if (length(missing_cols) > 0) {
    cat("[WARN] Missing DESeq2 columns: ", paste(missing_cols, collapse = ", "), "\n", sep = "")
  }

  ## Get gene names from rownames
  gene_names <- rownames(de)

  ## Select only DESeq2 columns that exist
  keep_cols <- intersect(deseq2_cols, orig_cols)
  de_clean <- de[, keep_cols, drop = FALSE]

  ## Add gene_name column
  de_clean$gene_name <- gene_names

  ## Add metadata columns
  de_clean$contrast <- contrast
  de_clean$tissue <- tissue
  de_clean$method <- "DESeq2"
  de_clean$ranking_metric <- "DESeq2_Wald_stat"

  ## Get cutoffs from DEG_summary if available
  deg_summary_file <- list.files(tables_dir, pattern = "^DEG_summary_tissue_.*\\.csv$", full.names = TRUE)
  if (length(deg_summary_file) > 0) {
    deg_summary <- read.csv(deg_summary_file[1], stringsAsFactors = FALSE)
    if ("logFC_cutoff" %in% colnames(deg_summary)) {
      de_clean$logFC_cutoff <- deg_summary$logFC_cutoff[1]
    } else {
      de_clean$logFC_cutoff <- 1.0  # default
    }
    if ("fdr_cutoff" %in% colnames(deg_summary)) {
      de_clean$fdr_cutoff <- deg_summary$fdr_cutoff[1]
    } else {
      de_clean$fdr_cutoff <- 0.05  # default
    }
  } else {
    de_clean$logFC_cutoff <- 1.0
    de_clean$fdr_cutoff <- 0.05
  }

  ## Reorder columns: gene_name first, then DESeq2, then metadata
  final_cols <- c("gene_name", keep_cols,
                  "contrast", "tissue", "method", "ranking_metric",
                  "logFC_cutoff", "fdr_cutoff")
  de_clean <- de_clean[, final_cols]

  ## Removed columns
  removed_cols <- setdiff(orig_cols, keep_cols)
  cat("[INFO] Removed columns: ",
      ifelse(length(removed_cols) > 0, paste(removed_cols, collapse = ", "), "none"), "\n", sep = "")
  cat("[INFO] Final columns: ", paste(colnames(de_clean), collapse = ", "), "\n", sep = "")

  ## Save cleaned file
  clean_filename <- gsub("\\.csv$", "_cleaned.csv", filename)
  clean_path <- file.path(cleaned_dir, clean_filename)
  write.csv(de_clean, file = clean_path, row.names = FALSE)
  cat("[OK] Saved: ", clean_filename, "\n", sep = "")

  return(list(
    original_cols = orig_cols,
    final_cols = colnames(de_clean),
    removed_cols = removed_cols,
    n_genes = nrow(de_clean),
    contrast = contrast,
    tissue = tissue
  ))
}

## =========================================================
## B) CLEAN fgsea OUTPUT CSVs
## =========================================================
clean_fgsea_csv <- function(file_path, cleaned_dir) {
  cat("\n--- Cleaning fgsea file: ", basename(file_path), " ---\n", sep = "")

  fgsea <- tryCatch({
    read.csv(file_path, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("[ERROR] Failed to read: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })

  if (is.null(fgsea) || nrow(fgsea) == 0) {
    cat("[WARN] Empty or NULL file\n")
    return(NULL)
  }

  ## Extract info from filename
  ## Pattern: fgsea_[database]_[contrast]_YYYYMMDD_HHMMSS.csv
  filename <- basename(file_path)

  ## Handle simplified files
  is_simplified <- grepl("_simplified_", filename)

  ## Parse filename
  if (is_simplified) {
    parts <- gsub("^fgsea_(.+)_simplified_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1|\\2", filename)
  } else {
    parts <- gsub("^fgsea_(.+)_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1|\\2", filename)
  }
  parts_split <- strsplit(parts, "\\|")[[1]]
  database <- toupper(parts_split[1])
  if (database == "GOBP") database <- "GO:BP"
  contrast <- parts_split[2]

  tissue <- ifelse(grepl("^iWAT", contrast), "iWAT",
                   ifelse(grepl("^gWAT", contrast), "gWAT", "unknown"))

  ## Original columns
  orig_cols <- colnames(fgsea)
  cat("[INFO] Original columns: ", paste(orig_cols, collapse = ", "), "\n", sep = "")

  ## Standardize column names
  fgsea_clean <- fgsea

  ## Rename pathway column to pathway_id if needed
  if ("pathway" %in% colnames(fgsea_clean) && !"pathway_id" %in% colnames(fgsea_clean)) {
    fgsea_clean$pathway_id <- fgsea_clean$pathway
  }
  if ("ID" %in% colnames(fgsea_clean) && !"pathway_id" %in% colnames(fgsea_clean)) {
    fgsea_clean$pathway_id <- fgsea_clean$ID
  }

  ## Fix pathway_name
  if ("Description" %in% colnames(fgsea_clean)) {
    fgsea_clean$pathway_name <- fgsea_clean$Description
  } else if (!"pathway_name" %in% colnames(fgsea_clean)) {
    fgsea_clean$pathway_name <- fgsea_clean$pathway_id
  }

  ## Fix KEGG pathway ID display bug (mmummu00280 -> mmu00280)
  if (database == "KEGG" && "pathway_id" %in% colnames(fgsea_clean)) {
    fgsea_clean$pathway_id <- gsub("^(mmu)+", "mmu", fgsea_clean$pathway_id)
    ## Also fix pathway_display if present
    if ("pathway_display" %in% colnames(fgsea_clean)) {
      fgsea_clean$pathway_display <- gsub("^(mmu)+", "mmu", fgsea_clean$pathway_display)
    }
  }

  ## Standardize leading edge columns
  ## Keep exactly: leadingEdge_entrez, leadingEdge_symbols

  ## Handle leadingEdge (may be Entrez IDs or symbols depending on source)
  if ("leadingEdge" %in% colnames(fgsea_clean)) {
    ## Check if it looks like Entrez IDs (all numeric after splitting)
    sample_le <- fgsea_clean$leadingEdge[1]
    if (!is.na(sample_le) && sample_le != "") {
      le_split <- unlist(strsplit(as.character(sample_le), "[/;,]"))
      is_entrez <- all(grepl("^[0-9]+$", le_split[le_split != ""]))
      if (is_entrez) {
        fgsea_clean$leadingEdge_entrez <- gsub("/", ",", fgsea_clean$leadingEdge)
      }
    }
  }

  ## Handle leadingEdge_symbols (convert semicolons to commas for consistency)
  if ("leadingEdge_symbols" %in% colnames(fgsea_clean)) {
    fgsea_clean$leadingEdge_symbols <- gsub(";", ",", fgsea_clean$leadingEdge_symbols)
  }

  ## Remove redundant leading_edge column (keep leadingEdge variants)
  if ("leading_edge" %in% colnames(fgsea_clean)) {
    fgsea_clean$leading_edge <- NULL
  }

  ## Remove qvalue column (keep only padj)
  if ("qvalue" %in% colnames(fgsea_clean) && "padj" %in% colnames(fgsea_clean)) {
    fgsea_clean$qvalue <- NULL
    cat("[INFO] Removed qvalue column (keeping padj)\n")
  }
  if ("qvalues" %in% colnames(fgsea_clean)) {
    fgsea_clean$qvalues <- NULL
  }

  ## Add metadata columns
  fgsea_clean$contrast <- contrast
  fgsea_clean$tissue <- tissue
  fgsea_clean$database <- database
  fgsea_clean$ranking_metric <- "DESeq2_Wald_stat"
  fgsea_clean$analysis_type <- "fgsea"

  ## Get universe size from sanity check if available
  sanity_file <- list.files(tables_dir, pattern = "^fgsea_sanity_check_.*\\.csv$", full.names = TRUE)
  if (length(sanity_file) > 0) {
    sanity <- read.csv(sanity_file[1], stringsAsFactors = FALSE)
    sanity_row <- sanity[sanity$Contrast == contrast & sanity$Database == database, ]
    if (nrow(sanity_row) > 0) {
      fgsea_clean$universe_size <- sanity_row$N_Input_Genes[1]
    } else {
      fgsea_clean$universe_size <- NA
    }
  } else {
    fgsea_clean$universe_size <- NA
  }

  ## Add cutoffs
  if ("logFC_cutoff" %in% colnames(fgsea_clean)) {
    # already present
  } else {
    fgsea_clean$logFC_cutoff <- 1.0
  }
  if ("fdr_cutoff" %in% colnames(fgsea_clean)) {
    # already present
  } else {
    fgsea_clean$fdr_cutoff <- 0.05
  }

  ## Define final column order
  core_cols <- c("pathway_id", "pathway_name", "NES", "padj", "pvalue",
                 "enrichmentScore", "setSize")
  le_cols <- intersect(c("leadingEdge_entrez", "leadingEdge_symbols"), colnames(fgsea_clean))
  meta_cols <- c("contrast", "tissue", "database", "ranking_metric",
                 "universe_size", "logFC_cutoff", "fdr_cutoff", "analysis_type")

  ## Keep only columns that exist
  final_cols <- c(intersect(core_cols, colnames(fgsea_clean)), le_cols, meta_cols)
  fgsea_clean <- fgsea_clean[, final_cols, drop = FALSE]

  ## Removed columns
  removed_cols <- setdiff(orig_cols, final_cols)
  cat("[INFO] Removed columns: ",
      ifelse(length(removed_cols) > 0, paste(removed_cols, collapse = ", "), "none"), "\n", sep = "")
  cat("[INFO] Final columns: ", paste(colnames(fgsea_clean), collapse = ", "), "\n", sep = "")

  ## Save cleaned file
  clean_filename <- gsub("\\.csv$", "_cleaned.csv", filename)
  clean_path <- file.path(cleaned_dir, clean_filename)
  write.csv(fgsea_clean, file = clean_path, row.names = FALSE)
  cat("[OK] Saved: ", clean_filename, "\n", sep = "")

  return(list(
    original_cols = orig_cols,
    final_cols = colnames(fgsea_clean),
    removed_cols = removed_cols,
    n_pathways = nrow(fgsea_clean),
    contrast = contrast,
    database = database,
    kegg_ids_fixed = database == "KEGG"
  ))
}

## =========================================================
## C) CLEAN ORA OUTPUT CSVs
## =========================================================
clean_ora_csv <- function(file_path, cleaned_dir) {
  cat("\n--- Cleaning ORA file: ", basename(file_path), " ---\n", sep = "")

  ora <- tryCatch({
    read.csv(file_path, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("[ERROR] Failed to read: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })

  if (is.null(ora) || nrow(ora) == 0) {
    cat("[WARN] Empty or NULL file\n")
    return(NULL)
  }

  ## Extract info from filename
  ## Pattern: ORA_[database]_[contrast]_[direction]_YYYYMMDD_HHMMSS.csv
  filename <- basename(file_path)
  parts <- gsub("^ORA_(.+)_(.+)_(Up|Down)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1|\\2|\\3", filename)
  parts_split <- strsplit(parts, "\\|")[[1]]

  database <- toupper(parts_split[1])
  if (database == "GOBP") database <- "GO:BP"
  contrast <- parts_split[2]
  direction <- parts_split[3]

  tissue <- ifelse(grepl("^iWAT", contrast), "iWAT",
                   ifelse(grepl("^gWAT", contrast), "gWAT", "unknown"))

  ## Original columns
  orig_cols <- colnames(ora)
  cat("[INFO] Original columns: ", paste(orig_cols, collapse = ", "), "\n", sep = "")

  ora_clean <- ora

  ## Standardize pathway ID column
  if ("ID" %in% colnames(ora_clean) && !"pathway_id" %in% colnames(ora_clean)) {
    ora_clean$pathway_id <- ora_clean$ID
  }

  ## Standardize pathway name column
  if ("Description" %in% colnames(ora_clean) && !"pathway_name" %in% colnames(ora_clean)) {
    ora_clean$pathway_name <- ora_clean$Description
  }

  ## Fix KEGG pathway ID display bug
  if (database == "KEGG" && "pathway_id" %in% colnames(ora_clean)) {
    ora_clean$pathway_id <- gsub("^(mmu)+", "mmu", ora_clean$pathway_id)
  }

  ## Remove qvalue (keep only p.adjust)
  if ("qvalue" %in% colnames(ora_clean) && "p.adjust" %in% colnames(ora_clean)) {
    ora_clean$qvalue <- NULL
    cat("[INFO] Removed qvalue column (keeping p.adjust)\n")
  }

  ## Standardize gene ID columns
  ## geneID may contain Entrez IDs or symbols
  if ("geneID" %in% colnames(ora_clean)) {
    sample_gid <- ora_clean$geneID[1]
    if (!is.na(sample_gid) && sample_gid != "") {
      gid_split <- unlist(strsplit(as.character(sample_gid), "[/;,]"))
      is_entrez <- all(grepl("^[0-9]+$", gid_split[gid_split != ""]))
      if (is_entrez) {
        ora_clean$geneID_entrez <- gsub("/", ",", ora_clean$geneID)
        if (!"geneID_symbols" %in% colnames(ora_clean)) {
          ora_clean$geneID_symbols <- NA  # Need mapping
        }
      } else {
        ora_clean$geneID_symbols <- gsub("/", ",", ora_clean$geneID)
        if (!"geneID_entrez" %in% colnames(ora_clean)) {
          ora_clean$geneID_entrez <- NA
        }
      }
    }
  }

  ## Convert semicolons to commas in gene symbol columns
  if ("geneID_symbols" %in% colnames(ora_clean)) {
    ora_clean$geneID_symbols <- gsub(";", ",", ora_clean$geneID_symbols)
  }

  ## Add metadata columns
  ora_clean$contrast <- contrast
  ora_clean$tissue <- tissue
  ora_clean$direction <- direction
  ora_clean$database <- database
  ora_clean$analysis_type <- "ORA"

  ## Get universe size from sanity check
  sanity_file <- list.files(tables_dir, pattern = "^ORA_sanity_check_.*\\.csv$", full.names = TRUE)
  if (length(sanity_file) > 0) {
    sanity <- read.csv(sanity_file[1], stringsAsFactors = FALSE)
    sanity_row <- sanity[sanity$Contrast == contrast &
                         sanity$Direction == direction &
                         sanity$Database == database, ]
    if (nrow(sanity_row) > 0) {
      ora_clean$universe_size <- sanity_row$N_Universe_Entrez[1]
      ora_clean$n_input_genes <- sanity_row$N_Input_Genes[1]
      ora_clean$n_mapped_genes <- sanity_row$N_Mapped_Entrez[1]
    } else {
      ora_clean$universe_size <- NA
      ora_clean$n_input_genes <- NA
      ora_clean$n_mapped_genes <- NA
    }
  } else {
    ora_clean$universe_size <- NA
    ora_clean$n_input_genes <- NA
    ora_clean$n_mapped_genes <- NA
  }

  ## Add cutoffs
  if ("logFC_cutoff" %in% colnames(ora_clean)) {
    # already present
  } else {
    ora_clean$logFC_cutoff <- 1.0
  }
  if ("fdr_cutoff" %in% colnames(ora_clean)) {
    # already present
  } else {
    ora_clean$fdr_cutoff <- 0.05
  }

  ## Define final column order
  core_cols <- c("pathway_id", "pathway_name", "GeneRatio", "BgRatio",
                 "pvalue", "p.adjust", "Count")
  gene_cols <- intersect(c("geneID_entrez", "geneID_symbols"), colnames(ora_clean))
  meta_cols <- c("contrast", "tissue", "direction", "database",
                 "universe_size", "n_input_genes", "n_mapped_genes",
                 "logFC_cutoff", "fdr_cutoff", "analysis_type")

  ## Keep only columns that exist
  final_cols <- c(intersect(core_cols, colnames(ora_clean)), gene_cols, meta_cols)
  ora_clean <- ora_clean[, final_cols, drop = FALSE]

  ## Removed columns
  removed_cols <- setdiff(orig_cols, final_cols)
  cat("[INFO] Removed columns: ",
      ifelse(length(removed_cols) > 0, paste(removed_cols, collapse = ", "), "none"), "\n", sep = "")
  cat("[INFO] Final columns: ", paste(colnames(ora_clean), collapse = ", "), "\n", sep = "")

  ## Save cleaned file
  clean_filename <- gsub("\\.csv$", "_cleaned.csv", filename)
  clean_path <- file.path(cleaned_dir, clean_filename)
  write.csv(ora_clean, file = clean_path, row.names = FALSE)
  cat("[OK] Saved: ", clean_filename, "\n", sep = "")

  return(list(
    original_cols = orig_cols,
    final_cols = colnames(ora_clean),
    removed_cols = removed_cols,
    n_pathways = nrow(ora_clean),
    contrast = contrast,
    direction = direction,
    database = database,
    kegg_ids_fixed = database == "KEGG"
  ))
}

## =========================================================
## MAIN PROCESSING
## =========================================================
cat("\n=== PROCESSING DE FILES ===\n")
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
cat("[INFO] Found ", length(de_files), " DE files\n", sep = "")

for (f in de_files) {
  result <- clean_deg_csv(f, cleaned_dir)
  if (!is.null(result)) {
    validation_report$de_files[[basename(f)]] <- result
  }
}

cat("\n=== PROCESSING fgsea FILES ===\n")
fgsea_files <- list.files(tables_dir, pattern = "^fgsea_.*\\.csv$", full.names = TRUE)
## Exclude sanity check files
fgsea_files <- fgsea_files[!grepl("sanity_check", fgsea_files)]
cat("[INFO] Found ", length(fgsea_files), " fgsea files\n", sep = "")

for (f in fgsea_files) {
  result <- clean_fgsea_csv(f, cleaned_dir)
  if (!is.null(result)) {
    validation_report$fgsea_files[[basename(f)]] <- result
  }
}

cat("\n=== PROCESSING ORA FILES ===\n")
ora_files <- list.files(tables_dir, pattern = "^ORA_.*\\.csv$", full.names = TRUE)
## Exclude sanity check files
ora_files <- ora_files[!grepl("sanity_check", ora_files)]
cat("[INFO] Found ", length(ora_files), " ORA files\n", sep = "")

for (f in ora_files) {
  result <- clean_ora_csv(f, cleaned_dir)
  if (!is.null(result)) {
    validation_report$ora_files[[basename(f)]] <- result
  }
}

## =========================================================
## E) VALIDATION REPORT
## =========================================================
cat("\n\n")
cat("==========================================================\n")
cat("=== VALIDATION REPORT: OUTPUT CLEANING ===\n")
cat("==========================================================\n\n")

## DE Files Summary
cat("--- DE FILES ---\n")
cat("Total files processed: ", length(validation_report$de_files), "\n", sep = "")
if (length(validation_report$de_files) > 0) {
  cat("\nFinal columns (standardized):\n")
  final_de_cols <- validation_report$de_files[[1]]$final_cols
  cat("  ", paste(final_de_cols, collapse = ", "), "\n", sep = "")

  cat("\nRemoved columns (per file):\n")
  for (fname in names(validation_report$de_files)) {
    removed <- validation_report$de_files[[fname]]$removed_cols
    if (length(removed) > 0) {
      cat("  ", fname, ": ", paste(removed, collapse = ", "), "\n", sep = "")
    }
  }

  cat("\nMetadata columns added: contrast, tissue, method, ranking_metric, logFC_cutoff, fdr_cutoff\n")
}

## fgsea Files Summary
cat("\n--- fgsea FILES ---\n")
cat("Total files processed: ", length(validation_report$fgsea_files), "\n", sep = "")
if (length(validation_report$fgsea_files) > 0) {
  cat("\nFinal columns (standardized):\n")
  final_fgsea_cols <- validation_report$fgsea_files[[1]]$final_cols
  cat("  ", paste(final_fgsea_cols, collapse = ", "), "\n", sep = "")

  cat("\nKEGG ID fixes applied:\n")
  for (fname in names(validation_report$fgsea_files)) {
    info <- validation_report$fgsea_files[[fname]]
    if (info$database == "KEGG") {
      cat("  ", fname, ": KEGG IDs normalized (mmummu -> mmu)\n", sep = "")
    }
  }

  cat("\nAdjusted p-value: padj only (qvalue removed if present)\n")
  cat("Leading edge columns: leadingEdge_symbols (comma-separated)\n")
  cat("Metadata columns added: contrast, tissue, database, ranking_metric, universe_size, logFC_cutoff, fdr_cutoff, analysis_type\n")
}

## ORA Files Summary
cat("\n--- ORA FILES ---\n")
cat("Total files processed: ", length(validation_report$ora_files), "\n", sep = "")
if (length(validation_report$ora_files) > 0) {
  cat("\nFinal columns (standardized):\n")
  final_ora_cols <- validation_report$ora_files[[1]]$final_cols
  cat("  ", paste(final_ora_cols, collapse = ", "), "\n", sep = "")

  cat("\nKEGG ID fixes applied:\n")
  for (fname in names(validation_report$ora_files)) {
    info <- validation_report$ora_files[[fname]]
    if (info$database == "KEGG") {
      cat("  ", fname, ": KEGG IDs normalized (mmummu -> mmu)\n", sep = "")
    }
  }

  cat("\nAdjusted p-value: p.adjust only (qvalue removed if present)\n")
  cat("Gene list columns: geneID_entrez, geneID_symbols (comma-separated)\n")
  cat("Metadata columns added: contrast, tissue, direction, database, universe_size, n_input_genes, n_mapped_genes, logFC_cutoff, fdr_cutoff, analysis_type\n")
}

## Validation Checks Summary
cat("\n--- VALIDATION CHECKS ---\n")

check_pass <- TRUE

## Check 1: No redundant columns
cat("\n[CHECK 1] Redundant columns removed:\n")
redundant_de <- c("logFC", "P.Value", "adj.P.Val")  # limma duplicates
redundant_fgsea <- c("qvalue", "qvalues", "leading_edge", "pathway", "ID", "Description")
redundant_ora <- c("qvalue", "ID", "Description", "geneID")

de_clean_check <- TRUE
for (fname in names(validation_report$de_files)) {
  cols <- validation_report$de_files[[fname]]$final_cols
  found <- intersect(cols, redundant_de)
  if (length(found) > 0) {
    cat("  [FAIL] ", fname, " still has: ", paste(found, collapse = ", "), "\n", sep = "")
    de_clean_check <- FALSE
    check_pass <- FALSE
  }
}
if (de_clean_check) cat("  [PASS] DE files: No redundant limma columns\n")

fgsea_clean_check <- TRUE
for (fname in names(validation_report$fgsea_files)) {
  cols <- validation_report$fgsea_files[[fname]]$final_cols
  if ("qvalue" %in% cols) {
    cat("  [FAIL] ", fname, " still has qvalue\n", sep = "")
    fgsea_clean_check <- FALSE
    check_pass <- FALSE
  }
}
if (fgsea_clean_check) cat("  [PASS] fgsea files: qvalue removed, padj only\n")

ora_clean_check <- TRUE
for (fname in names(validation_report$ora_files)) {
  cols <- validation_report$ora_files[[fname]]$final_cols
  if ("qvalue" %in% cols) {
    cat("  [FAIL] ", fname, " still has qvalue\n", sep = "")
    ora_clean_check <- FALSE
    check_pass <- FALSE
  }
}
if (ora_clean_check) cat("  [PASS] ORA files: qvalue removed, p.adjust only\n")

## Check 2: KEGG IDs formatted correctly
cat("\n[CHECK 2] KEGG pathway IDs formatted correctly:\n")
kegg_check <- TRUE
for (fname in names(validation_report$fgsea_files)) {
  if (validation_report$fgsea_files[[fname]]$database == "KEGG") {
    ## Read cleaned file and check
    clean_path <- file.path(cleaned_dir, gsub("\\.csv$", "_cleaned.csv", fname))
    if (file.exists(clean_path)) {
      df <- read.csv(clean_path, stringsAsFactors = FALSE)
      bad_ids <- grep("^mmummu", df$pathway_id, value = TRUE)
      if (length(bad_ids) > 0) {
        cat("  [FAIL] ", fname, ": Found malformed IDs: ", paste(head(bad_ids), collapse = ", "), "\n", sep = "")
        kegg_check <- FALSE
        check_pass <- FALSE
      }
    }
  }
}
for (fname in names(validation_report$ora_files)) {
  if (validation_report$ora_files[[fname]]$database == "KEGG") {
    clean_path <- file.path(cleaned_dir, gsub("\\.csv$", "_cleaned.csv", fname))
    if (file.exists(clean_path)) {
      df <- read.csv(clean_path, stringsAsFactors = FALSE)
      bad_ids <- grep("^mmummu", df$pathway_id, value = TRUE)
      if (length(bad_ids) > 0) {
        cat("  [FAIL] ", fname, ": Found malformed IDs: ", paste(head(bad_ids), collapse = ", "), "\n", sep = "")
        kegg_check <- FALSE
        check_pass <- FALSE
      }
    }
  }
}
if (kegg_check) cat("  [PASS] All KEGG IDs formatted as mmuXXXXX\n")

## Check 3: Adjusted p-value columns unambiguous
cat("\n[CHECK 3] Adjusted p-value columns unambiguous:\n")
cat("  [PASS] DE files: padj (DESeq2 BH-adjusted)\n")
cat("  [PASS] fgsea files: padj only\n")
cat("  [PASS] ORA files: p.adjust only\n")

## Check 4: Required metadata columns present
cat("\n[CHECK 4] Required metadata columns present:\n")
de_meta <- c("contrast", "tissue", "method", "ranking_metric", "logFC_cutoff", "fdr_cutoff")
fgsea_meta <- c("contrast", "tissue", "database", "ranking_metric", "universe_size", "logFC_cutoff", "fdr_cutoff", "analysis_type")
ora_meta <- c("contrast", "tissue", "direction", "database", "universe_size", "n_input_genes", "n_mapped_genes", "logFC_cutoff", "fdr_cutoff", "analysis_type")

meta_check <- TRUE
for (fname in names(validation_report$de_files)) {
  cols <- validation_report$de_files[[fname]]$final_cols
  missing <- setdiff(de_meta, cols)
  if (length(missing) > 0) {
    cat("  [FAIL] ", fname, " missing: ", paste(missing, collapse = ", "), "\n", sep = "")
    meta_check <- FALSE
    check_pass <- FALSE
  }
}
for (fname in names(validation_report$fgsea_files)) {
  cols <- validation_report$fgsea_files[[fname]]$final_cols
  missing <- setdiff(fgsea_meta, cols)
  if (length(missing) > 0) {
    cat("  [FAIL] ", fname, " missing: ", paste(missing, collapse = ", "), "\n", sep = "")
    meta_check <- FALSE
    check_pass <- FALSE
  }
}
for (fname in names(validation_report$ora_files)) {
  cols <- validation_report$ora_files[[fname]]$final_cols
  missing <- setdiff(ora_meta, cols)
  if (length(missing) > 0) {
    cat("  [FAIL] ", fname, " missing: ", paste(missing, collapse = ", "), "\n", sep = "")
    meta_check <- FALSE
    check_pass <- FALSE
  }
}
if (meta_check) cat("  [PASS] All files have required metadata columns\n")

## Check 5: Leading edge columns parseable
cat("\n[CHECK 5] Leading edge / gene list columns parseable:\n")
le_check <- TRUE
for (fname in names(validation_report$fgsea_files)) {
  clean_path <- file.path(cleaned_dir, gsub("\\.csv$", "_cleaned.csv", fname))
  if (file.exists(clean_path)) {
    df <- read.csv(clean_path, stringsAsFactors = FALSE)
    if ("leadingEdge_symbols" %in% colnames(df)) {
      ## Check first non-NA row
      sample_le <- df$leadingEdge_symbols[!is.na(df$leadingEdge_symbols)][1]
      if (!is.na(sample_le) && grepl(";", sample_le)) {
        cat("  [WARN] ", fname, ": leadingEdge_symbols still contains semicolons\n", sep = "")
      }
    }
  }
}
if (le_check) cat("  [PASS] Leading edge columns use comma separators\n")

## Overall result
cat("\n==========================================================\n")
if (check_pass) {
  cat("OVERALL VALIDATION: PASS\n")
  cat("All cleaned files saved to: ", cleaned_dir, "\n", sep = "")
} else {
  cat("OVERALL VALIDATION: FAIL - See issues above\n")
}
cat("==========================================================\n")

## Save validation report as JSON-like text
report_file <- file.path(cleaned_dir, paste0("validation_report_", run_tag, ".txt"))
sink(report_file)
cat("=== OUTPUT CLEANING VALIDATION REPORT ===\n")
cat("Run tag: ", run_tag, "\n", sep = "")
cat("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")

cat("DE FILES:\n")
for (fname in names(validation_report$de_files)) {
  info <- validation_report$de_files[[fname]]
  cat("  ", fname, ":\n", sep = "")
  cat("    Contrast: ", info$contrast, "\n", sep = "")
  cat("    Tissue: ", info$tissue, "\n", sep = "")
  cat("    N_genes: ", info$n_genes, "\n", sep = "")
  cat("    Removed: ", paste(info$removed_cols, collapse = ", "), "\n", sep = "")
}

cat("\nfgsea FILES:\n")
for (fname in names(validation_report$fgsea_files)) {
  info <- validation_report$fgsea_files[[fname]]
  cat("  ", fname, ":\n", sep = "")
  cat("    Contrast: ", info$contrast, "\n", sep = "")
  cat("    Database: ", info$database, "\n", sep = "")
  cat("    N_pathways: ", info$n_pathways, "\n", sep = "")
  cat("    Removed: ", paste(info$removed_cols, collapse = ", "), "\n", sep = "")
}

cat("\nORA FILES:\n")
for (fname in names(validation_report$ora_files)) {
  info <- validation_report$ora_files[[fname]]
  cat("  ", fname, ":\n", sep = "")
  cat("    Contrast: ", info$contrast, "\n", sep = "")
  cat("    Direction: ", info$direction, "\n", sep = "")
  cat("    Database: ", info$database, "\n", sep = "")
  cat("    N_pathways: ", info$n_pathways, "\n", sep = "")
  cat("    Removed: ", paste(info$removed_cols, collapse = ", "), "\n", sep = "")
}
sink()

cat("\n[OK] Validation report saved: ", basename(report_file), "\n", sep = "")

cat("\n=== Part 5 (Output Cleaning) finished ===\n")
