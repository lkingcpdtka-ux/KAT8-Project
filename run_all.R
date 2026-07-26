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
##   results that look fine but are stale or incomplete. This driver makes the
##   order explicit and fails loudly.
##
## USAGE
##   Rscript run_all.R                 # tissue pipeline (parts 1-4)
##   Rscript run_all.R all             # tissue + cells + downstream mechanism
##   Rscript run_all.R tissue          # parts 1-4 only
##   Rscript run_all.R cells           # cells DE/ORA/GSEA only
##   Rscript run_all.R downstream      # TF activity, concordance, mechanism
##                                     #   (requires tissue + cells to have run)
## =========================================================

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args) == 0) "tissue" else tolower(args[1])

PROJECT_ROOT <- getwd()
cat("Project root: ", PROJECT_ROOT, "\n", sep = "")
if (!file.exists(file.path(PROJECT_ROOT, "parameters.R")))
  stop("parameters.R not found in ", PROJECT_ROOT,
       " -- run this from the project root.", call. = FALSE)

## Each step: label, path, and whether a missing file is fatal.
TISSUE <- list(
  list(label = "Part 1  DE / QC / DEG lists",      path = "part1_main_analysis.R",  required = TRUE),
  list(label = "Part 2  ORA",                      path = "part2_ora.R",            required = TRUE),
  list(label = "Part 3  GSEA",                     path = "part3_fgsea.R",          required = TRUE),
  list(label = "Part 4  Visualisation",            path = "part4_visualization.R",  required = TRUE)
)
CELLS <- list(
  list(label = "Cells   DE / ORA / GSEA",          path = "KAT8_bulk_cells_comprehensive.R", required = TRUE)
)
DOWNSTREAM <- list(
  list(label = "Part 4  TF activity (all contrasts)",
       path = file.path("downstream_analysis", "part4_tf_activity.R"), required = FALSE),
  list(label = "Part 4b Cell<->tissue concordance",
       path = file.path("downstream_analysis", "part4b_cell_tissue_concordance.R"), required = FALSE),
  list(label = "Part 5  Cell-autonomous mechanism",
       path = file.path("downstream_analysis", "part5_cells_mechanism.R"), required = FALSE)
)

steps <- switch(what,
  "tissue"     = TISSUE,
  "cells"      = CELLS,
  "downstream" = DOWNSTREAM,
  "all"        = c(TISSUE, CELLS, DOWNSTREAM),
  stop("Unknown target '", what, "'. Use: tissue | cells | downstream | all", call. = FALSE))

cat("\n==========================================================\n")
cat("RUNNING: ", toupper(what), "  (", length(steps), " steps)\n", sep = "")
cat("==========================================================\n")

results <- data.frame(step = character(), status = character(),
                      seconds = numeric(), stringsAsFactors = FALSE)
t_all <- Sys.time()

for (s in steps) {
  cat("\n----------------------------------------------------------\n")
  cat(">>> ", s$label, "   [", s$path, "]\n", sep = "")
  cat("----------------------------------------------------------\n")

  if (!file.exists(s$path)) {
    msg <- paste0("script not found: ", s$path)
    if (isTRUE(s$required)) stop("FATAL - ", msg, call. = FALSE)
    cat("[SKIP] ", msg, "\n", sep = "")
    results <- rbind(results, data.frame(step = s$label, status = "SKIPPED", seconds = 0))
    next
  }

  t0 <- Sys.time()
  ## Run each part in its OWN environment so leftover variables from one script
  ## cannot leak into the next and silently change its behaviour.
  env <- new.env(parent = globalenv())
  ok <- tryCatch({
    source(s$path, local = env, echo = FALSE)
    TRUE
  }, error = function(e) {
    cat("\n[FAIL] ", s$label, "\n  ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  results <- rbind(results, data.frame(step = s$label,
                                       status = ifelse(ok, "OK", "FAILED"),
                                       seconds = round(el, 1)))
  if (!ok) {
    cat("\n==========================================================\n")
    cat("PIPELINE HALTED at: ", s$label, "\n", sep = "")
    cat("  Downstream steps were NOT run, so no stale/partial results\n")
    cat("  are produced. Fix the error above and re-run.\n")
    cat("==========================================================\n")
    print(results, row.names = FALSE)
    quit(status = 1)
  }
  cat("[OK] ", s$label, "  (", round(el, 1), "s)\n", sep = "")
}

cat("\n==========================================================\n")
cat("PIPELINE COMPLETE  (", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1),
    " min total)\n", sep = "")
cat("==========================================================\n")
print(results, row.names = FALSE)

## Point at the run folder that was just produced, so it is obvious which
## results are current.
ptr <- file.path(PROJECT_ROOT, "savepoints", "LATEST_RUN.txt")
if (file.exists(ptr)) {
  cat("\nCurrent run folder: ", trimws(readLines(ptr, warn = FALSE))[1], "\n", sep = "")
}
cat("\nReminder: check the SANITY / VALIDATION blocks in each part's console\n")
cat("output before interpreting results (sample annotation, DEG counts,\n")
cat("cross-validation between Parts 1-3).\n")
