#!/usr/bin/env Rscript

## =========================================================
## KAT8 RNA-seq - PIPELINE DRIVER
## =========================================================
## Runs the pipeline in DEPENDENCY ORDER and stops at the first failure.
##
## WHY THIS EXISTS
##   The parts must run in order: Part 1 writes the DE tables, the authoritative
##   DEG lists and the run-folder pointer that Parts 2/3/4 consume. Running them
##   out of order, or continuing after a mid-chain failure, leaves a partially
##   populated run folder that the next script will happily read -- producing
##   results that look fine but are stale or incomplete.
##
## HOW TO RUN  (both work; sourcing is SAFE and will not close your session)
##
##   From RStudio / the R console:
##       setwd("path/to/KAT8KD_RNAseq")     # must be the project root
##       source("run_all.R")                 # runs everything
##
##       # to run only part of it, set this BEFORE sourcing:
##       RUN_TARGET <- "tissue"; source("run_all.R")
##
##   From a terminal:
##       Rscript run_all.R                   # everything
##       Rscript run_all.R tissue            # parts 1-4 only
##       Rscript run_all.R cells             # cells DE/ORA/GSEA only
##       Rscript run_all.R downstream        # the four mechanism scripts only
##
## Targets: all (default) | tissue | cells | downstream
##
## Before running anything the driver PRE-FLIGHTS: it checks every script it is
## about to run actually exists, so a missing file fails in one second with a
## clear list, instead of after Part 1 has spent ten minutes on DESeq2.
##
## Between the cells step and the downstream steps it AUTOMATICALLY copies the
## newest DE tables from savepoints/RUN_*/tables into downstream_analysis/data/
## under the canonical names the downstream scripts expect.
## =========================================================

run_pipeline <- function(target = NULL, project_root = getwd()) {

  ## ---- resolve target -------------------------------------------------
  if (is.null(target)) {
    args <- commandArgs(trailingOnly = TRUE)
    target <- if (length(args) > 0) tolower(args[1]) else
              if (exists("RUN_TARGET", envir = globalenv()))
                tolower(get("RUN_TARGET", envir = globalenv())) else "all"
  }
  target <- tolower(target)

  cat("Project root : ", project_root, "\n", sep = "")
  cat("Target       : ", target, "\n", sep = "")

  TISSUE <- list(
    list(label = "Part 1  DE / QC / DEG lists", path = "part1_main_analysis.R", required = TRUE),
    list(label = "Part 2  ORA",                 path = "part2_ora.R",           required = TRUE),
    list(label = "Part 3  GSEA",                path = "part3_fgsea.R",         required = TRUE),
    list(label = "Part 4  Visualisation",       path = "part4_visualization.R", required = TRUE))
  CELLS <- list(
    list(label = "Cells   DE / ORA / GSEA", path = "KAT8_bulk_cells_comprehensive.R", required = TRUE))
  DOWNSTREAM <- list(
    list(label = "Part 4b  Is it cell-autonomous?  (concordance + programs)",
         path = "downstream_analysis/part4b_cell_tissue_concordance.R", required = TRUE),
    list(label = "Part 5a  What does the cell broadcast?  (secretome)",
         path = "downstream_analysis/part5a_secretome.R", required = TRUE),
    list(label = "Part 5b  Which TF drives it?  (3 layers + PROGENy)",
         path = "downstream_analysis/part5b_tf_analysis.R", required = TRUE),
    list(label = "Part 5c  Is the KAT8/NSL program engaged?  (complex output)",
         path = "downstream_analysis/part5c_kat8_complex.R", required = TRUE))

  steps <- switch(target,
    "tissue"     = TISSUE,
    "cells"      = CELLS,
    "downstream" = DOWNSTREAM,
    "all"        = c(TISSUE, CELLS, DOWNSTREAM),
    NULL)
  if (is.null(steps)) {
    cat("\n[ERROR] Unknown target '", target, "'. Use: all | tissue | cells | downstream\n", sep = "")
    return(invisible(NULL))
  }

  ## ---- PRE-FLIGHT: fail in one second, not ten minutes ----------------
  cat("\n--- pre-flight ---\n")
  problems <- character(0)
  if (!file.exists(file.path(project_root, "parameters.R")))
    problems <- c(problems, "parameters.R  (are you in the project root?)")
  for (s in steps)
    if (!file.exists(file.path(project_root, s$path))) problems <- c(problems, s$path)

  if (length(problems) > 0) {
    cat("[ERROR] These required files are missing from the project root:\n")
    for (p in problems) cat("   - ", p, "\n", sep = "")
    cat("\n  Working directory is:\n    ", project_root, "\n", sep = "")
    present <- list.files(project_root, pattern = "\\.R$")
    cat("  R scripts actually here (", length(present), "):\n", sep = "")
    if (length(present)) cat("    ", paste(present, collapse = "\n    "), "\n", sep = "") else
      cat("    (none)\n")
    dsa <- file.path(project_root, "downstream_analysis")
    cat("  downstream_analysis/ exists: ", dir.exists(dsa), "\n", sep = "")
    if (dir.exists(dsa)) cat("    contains: ",
        paste(list.files(dsa, pattern = "\\.R$"), collapse = ", "), "\n", sep = "")
    cat("\n  Fix: setwd() to the folder holding parameters.R and the part*.R files,\n")
    cat("  and make sure you downloaded the whole repo (including downstream_analysis/).\n")
    return(invisible(NULL))
  }
  cat("[OK] all ", length(steps), " scripts found\n", sep = "")

  ## ---- staging helper --------------------------------------------------
  stage_downstream_inputs <- function() {
    sp   <- file.path(project_root, "savepoints")
    dest <- file.path(project_root, "downstream_analysis", "data")
    if (!dir.exists(sp)) { cat("[WARN] no savepoints/ yet - nothing to stage\n"); return(FALSE) }
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    newest <- function(pat) {
      hits <- list.files(sp, pattern = pat, recursive = TRUE, full.names = TRUE)
      if (length(hits) == 0) return(NA_character_)
      hits[order(basename(hits), decreasing = TRUE)][1]   ## RUN tag sorts chronologically
    }
    wanted <- list(
      c("^DE_cells_KAT8KD_vs_CTL_.*\\.csv$",   "DE_cells_KAT8KD_vs_CTL.csv"),
      c("^DE_tissue_iWAT_KD_vs_CTL_.*\\.csv$", "DE_tissue_iWAT_KD_vs_CTL.csv"),
      c("^DE_tissue_gWAT_KD_vs_CTL_.*\\.csv$", "DE_tissue_gWAT_KD_vs_CTL.csv"))
    cat("\n--- staging downstream inputs ---\n")
    ok <- TRUE
    for (w in wanted) {
      src <- newest(w[1])
      if (is.na(src)) { cat("[MISS] ", w[2], "\n", sep = ""); ok <- FALSE; next }
      file.copy(src, file.path(dest, w[2]), overwrite = TRUE)
      cat("[OK]   ", w[2], "  <- ", basename(src), "\n", sep = "")
    }
    if (!ok) cat("[WARN] some DE tables not found; downstream steps will be skipped.\n")
    ok
  }

  ## ---- run -------------------------------------------------------------
  cat("\n==========================================================\n")
  cat("RUNNING: ", toupper(target), "  (", length(steps), " steps)\n", sep = "")
  cat("==========================================================\n")
  for (s in steps) cat("   - ", s$label, "\n", sep = "")

  is_down <- function(p) grepl("^downstream_analysis/", p)
  staged  <- if (identical(target, "downstream")) stage_downstream_inputs() else NA
  results <- data.frame(step = character(), status = character(),
                        seconds = numeric(), stringsAsFactors = FALSE)
  t_all <- Sys.time()

  for (s in steps) {
    if (is_down(s$path) && is.na(staged)) staged <- stage_downstream_inputs()
    if (is_down(s$path) && isFALSE(staged)) {
      cat("\n[SKIP] ", s$label, " - inputs not staged\n", sep = "")
      results <- rbind(results, data.frame(step = s$label, status = "SKIPPED", seconds = 0))
      next
    }

    cat("\n----------------------------------------------------------\n")
    cat(">>> ", s$label, "   [", s$path, "]\n", sep = "")
    cat("----------------------------------------------------------\n")

    t0 <- Sys.time()
    ## Each part runs in its OWN environment so variables cannot leak between
    ## scripts and silently change behaviour.
    env <- new.env(parent = globalenv())
    ok <- tryCatch({ source(file.path(project_root, s$path), local = env, echo = FALSE); TRUE },
                   error = function(e) {
                     cat("\n[FAIL] ", s$label, "\n  ", conditionMessage(e), "\n", sep = ""); FALSE })
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    results <- rbind(results, data.frame(step = s$label,
                                         status = ifelse(ok, "OK", "FAILED"), seconds = el))
    if (!ok) {
      cat("\n==========================================================\n")
      cat("PIPELINE HALTED at: ", s$label, "\n", sep = "")
      cat("  Later steps were NOT run, so no stale/partial results were produced.\n")
      cat("  Fix the error above and re-run.\n")
      cat("==========================================================\n")
      print(results, row.names = FALSE)
      ## NOTE: return, never quit(). quit() would terminate an interactive
      ## RStudio session and lose the user's workspace.
      return(invisible(results))
    }
    cat("[OK] ", s$label, "  (", el, "s)\n", sep = "")
  }

  cat("\n==========================================================\n")
  cat("PIPELINE COMPLETE  (",
      round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min)\n", sep = "")
  cat("==========================================================\n")
  print(results, row.names = FALSE)

  ptr <- file.path(project_root, "savepoints", "LATEST_RUN.txt")
  if (file.exists(ptr))
    cat("\nCurrent run folder: ", trimws(readLines(ptr, warn = FALSE))[1], "\n", sep = "")
  cat("\nBefore interpreting anything, check the SANITY / VALIDATION blocks in the\n")
  cat("console output above (sample annotation, DEG counts, Part 1-3 cross-validation).\n")
  invisible(results)
}

## ---------------------------------------------------------
## Entry point
## ---------------------------------------------------------
## Safe to source() interactively: on failure it RETURNS, it does not quit().
## Only a non-interactive `Rscript run_all.R` sets a failing exit status.
.res <- run_pipeline()
if (!interactive() && !is.null(.res) && any(.res$status %in% c("FAILED"))) quit(status = 1)
