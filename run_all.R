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
##   Rscript run_all.R                 # EVERYTHING (tissue + cells + downstream)
##   Rscript run_all.R tissue          # parts 1-4 only
##   Rscript run_all.R cells           # cells DE/ORA/GSEA only
##   Rscript run_all.R downstream      # the four mechanism scripts only
##                                     #   (needs tissue + cells to have run first)
##
## Between the cells step and the downstream steps the driver AUTOMATICALLY
## copies the newest DE tables from savepoints/RUN_*/tables into
## downstream_analysis/data/ under the canonical names the downstream scripts
## expect. Without that staging step the downstream scripts cannot find their
## input, which is the one thing that stops this running end to end.
## =========================================================

args <- commandArgs(trailingOnly = TRUE)
what <- if (length(args) == 0) "all" else tolower(args[1])

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
## One script per question, run in this order.
DOWNSTREAM <- list(
  list(label = "Part 4b  Is it cell-autonomous?  (concordance + programs)",
       path = file.path("downstream_analysis", "part4b_cell_tissue_concordance.R"), required = FALSE),
  list(label = "Part 5a  What does the cell broadcast?  (secretome)",
       path = file.path("downstream_analysis", "part5a_secretome.R"), required = FALSE),
  list(label = "Part 5b  Which TF drives it?  (3 layers + PROGENy)",
       path = file.path("downstream_analysis", "part5b_tf_analysis.R"), required = FALSE),
  list(label = "Part 5c  Is the KAT8/NSL program engaged?  (complex output)",
       path = file.path("downstream_analysis", "part5c_kat8_complex.R"), required = FALSE)
)

## ---------------------------------------------------------
## Stage downstream inputs
## ---------------------------------------------------------
## The tissue and cells scripts write DE tables into their own timestamped
## savepoints/RUN_*/tables folders (and the cells run is a SEPARATE folder from
## the tissue run). The downstream scripts read fixed paths under
## downstream_analysis/data/. This copies the newest of each across, so the
## whole pipeline runs end to end without a manual file-shuffling step.
stage_downstream_inputs <- function(root = PROJECT_ROOT) {
  sp <- file.path(root, "savepoints")
  dest <- file.path(root, "downstream_analysis", "data")
  if (!dir.exists(sp)) { cat("[WARN] no savepoints/ yet - nothing to stage\n"); return(FALSE) }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  ## newest file matching a pattern anywhere under savepoints/
  newest <- function(pat) {
    hits <- list.files(sp, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(hits) == 0) return(NA_character_)
    hits[order(basename(hits), decreasing = TRUE)][1]   # RUN tag sorts chronologically
  }
  wanted <- list(
    c("^DE_cells_KAT8KD_vs_CTL_.*\\.csv$",    "DE_cells_KAT8KD_vs_CTL.csv"),
    c("^DE_tissue_iWAT_KD_vs_CTL_.*\\.csv$",  "DE_tissue_iWAT_KD_vs_CTL.csv"),
    c("^DE_tissue_gWAT_KD_vs_CTL_.*\\.csv$",  "DE_tissue_gWAT_KD_vs_CTL.csv"))

  cat("\n--- staging downstream inputs ---\n")
  ok <- TRUE
  for (w in wanted) {
    src <- newest(w[1])
    if (is.na(src)) {
      cat("[MISS] ", w[2], "  (no match for ", w[1], ")\n", sep = ""); ok <- FALSE; next
    }
    file.copy(src, file.path(dest, w[2]), overwrite = TRUE)
    cat("[OK]   ", w[2], "  <- ", basename(src), "\n", sep = "")
  }
  if (!ok) cat("[WARN] some DE tables were not found; downstream steps will be skipped.\n")
  ok
}

steps <- switch(what,
  "tissue"     = TISSUE,
  "cells"      = CELLS,
  "downstream" = DOWNSTREAM,
  "all"        = c(TISSUE, CELLS, DOWNSTREAM),
  stop("Unknown target '", what, "'. Use: tissue | cells | downstream | all", call. = FALSE))

cat("\n==========================================================\n")
cat("RUNNING: ", toupper(what), "  (", length(steps), " steps)\n", sep = "")
cat("==========================================================\n")
for (s in steps) cat("   - ", s$label, "\n", sep = "")

## If downstream steps are part of this run, stage their inputs first (for
## "downstream" alone) or after the upstream steps finish (for "all").
STAGE_BEFORE <- identical(what, "downstream")
if (STAGE_BEFORE) staged_ok <- stage_downstream_inputs() else staged_ok <- NA

results <- data.frame(step = character(), status = character(),
                      seconds = numeric(), stringsAsFactors = FALSE)
t_all <- Sys.time()

is_downstream <- function(p) grepl("^downstream_analysis/", p)

for (s in steps) {
  ## First downstream step in a combined run: stage inputs now that the
  ## upstream scripts have produced fresh DE tables.
  if (is_downstream(s$path) && is.na(staged_ok)) staged_ok <- stage_downstream_inputs()
  if (is_downstream(s$path) && isFALSE(staged_ok)) {
    cat("\n[SKIP] ", s$label, " - required DE tables were not staged\n", sep = "")
    results <- rbind(results, data.frame(step = s$label, status = "SKIPPED", seconds = 0))
    next
  }

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
