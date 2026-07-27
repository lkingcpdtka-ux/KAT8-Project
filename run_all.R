#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   pipeline driver; one-shot RUN_* flags; run-scoped staging
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## KAT8 RNA-seq - PIPELINE DRIVER
## =========================================================
## Runs the pipeline in dependency order, RESUMING where it left off.
##
## RESUME BEHAVIOUR (the default)
##   The driver reuses the most recent savepoints/RUN_* folder and SKIPS any
##   step that already completed, so a failure late in the pipeline no longer
##   costs you another full ORA + GSEA run.
##
##   A step counts as complete ONLY if this driver recorded it finishing, in
##   RUN_<tag>/.done/<key>.done. Existing output files are deliberately NOT
##   treated as proof: reusing an old completed run folder would then skip
##   steps whose figures are stale. The plan printed at the start shows the
##   completion time of every step it intends to skip.
##
## HOW TO RUN  (sourcing is SAFE - it will not close your R session)
##
##   setwd("path/to/KAT8KD_RNAseq")   # the folder containing parameters.R
##   source("run_all.R")              # resume: run only what is missing
##
##   RUN_FRESH <- TRUE;  source("run_all.R")   # ignore old runs, new folder
##   RUN_FORCE <- "all"; source("run_all.R")   # same folder, redo everything
##   RUN_FORCE <- c("part2","part3"); source("run_all.R")  # redo just those
##   RUN_TARGET <- "tissue"; source("run_all.R")           # subset of steps
##
##   Terminal equivalents:
##     Rscript run_all.R                 # resume everything
##     Rscript run_all.R all fresh       # force a brand-new run folder
##     Rscript run_all.R tissue          # parts 1-4 only
##
## STEP KEYS for RUN_FORCE:
##   part1  qc  part2  part3  part4  cells  p4b  p5a  p5b  p5c   (or "all")
##
## TARGETS: all (default) | tissue | cells | downstream
##
## OUTPUT: one timestamped folder holds the whole run --
##   savepoints/RUN_<tag>/{tissue,cells,downstream,data}
## =========================================================

run_pipeline <- function(target = NULL, fresh = NULL, force = NULL,
                         project_root = getwd()) {

  ## ---- resolve options -------------------------------------------------
  gv <- function(nm) if (exists(nm, envir = globalenv())) get(nm, envir = globalenv()) else NULL
  args <- commandArgs(trailingOnly = TRUE)
  if (is.null(target)) target <- if (length(args) > 0) tolower(args[1]) else
                                 if (!is.null(gv("RUN_TARGET"))) tolower(gv("RUN_TARGET")) else "all"
  if (is.null(fresh))  fresh  <- ("fresh" %in% tolower(args)) || isTRUE(gv("RUN_FRESH"))
  if (is.null(force))  force  <- gv("RUN_FORCE")
  force <- if (is.null(force)) character(0) else tolower(as.character(force))
  target <- tolower(target)

  ## RUN_FRESH / RUN_FORCE / RUN_TARGET are ONE-SHOT.
  ##
  ## They live in the global environment and used to persist there for the
  ## whole R session. Setting RUN_FRESH <- TRUE for one run and then asking
  ## for RUN_FORCE <- "part4" later did NOT do what it looks like: RUN_FRESH
  ## was still set, so fresh mode won, a new folder was created, no step had
  ## a completion marker, and the ENTIRE pipeline re-ran instead of just
  ## Part 4. Reading them here and clearing them immediately makes each flag
  ## apply to exactly the run you set it for.
  used <- intersect(c("RUN_FRESH", "RUN_FORCE", "RUN_TARGET"), ls(globalenv()))
  if (length(used)) {
    rm(list = used, envir = globalenv())
    cat("[INFO] applied and cleared: ", paste(used, collapse = ", "),
        "  (these are one-shot; set again for the next run)\n", sep = "")
  }

  cat("Project root : ", project_root, "\n", sep = "")
  cat("Target       : ", target, "\n", sep = "")
  cat("Mode         : ", if (fresh) "FRESH (new run folder)" else "RESUME (reuse latest)", "\n", sep = "")
  if (length(force)) cat("Force re-run : ", paste(force, collapse = ", "), "\n", sep = "")

  ## A fresh run costs ~30-40 minutes. If that was not the intention it should
  ## be obvious BEFORE the first slow step, not thirty minutes in.
  if (fresh) {
    cat("\n**********************************************************\n")
    cat("FRESH RUN: every step will be recomputed (~30-40 min).\n")
    cat("Completed work in previous run folders is IGNORED, not reused.\n")
    if (length(force))
      cat("NOTE: RUN_FORCE (", paste(force, collapse = ", "),
          ") has no effect here -- a fresh folder has nothing to skip.\n", sep = "")
    cat("To run only some steps instead: press Esc now, then use\n")
    cat("  RUN_FORCE <- \"part4\"   (without RUN_FRESH)\n")
    cat("**********************************************************\n")
  }

  ## ---- step definitions -------------------------------------------------
  ## Completion is recorded ONLY by a marker the driver writes after a step
  ## actually succeeds (RUN_<tag>/.done/<key>.done).
  ##
  ## An earlier version also treated "a matching output file exists" as done.
  ## That silently mis-fired: reusing a months-old COMPLETED run folder made
  ## Part 4 look finished because old .png files were sitting in plots/, so the
  ## PCA / MDS / density / volcano / heatmap figures were never regenerated.
  ## Markers are unambiguous -- they mean "this driver ran this step here".
  ## fatal = FALSE marks an ADVISORY step (diagnostics): if it fails the
  ## pipeline records it and carries on, rather than blocking the analysis.
  st <- function(key, label, path, fatal = TRUE)
    list(key = key, label = label, path = path, fatal = fatal)

  TISSUE <- list(
    st("part1", "Part 1  DE / QC / DEG lists", "part1_main_analysis.R"),
    st("qc",    "Part 1b QC diagnostics",       "part1b_qc_diagnostics.R", fatal = FALSE),
    st("part2", "Part 2  ORA  (slow)",         "part2_ora.R"),
    st("part3", "Part 3  GSEA (slow)",         "part3_fgsea.R"),
    st("part4", "Part 4  Visualisation",       "part4_visualization.R"))
  CELLS <- list(
    st("cells", "Cells   DE / ORA / GSEA (slow)", "KAT8_bulk_cells_comprehensive.R"))
  DOWNSTREAM <- list(
    st("p4b", "Part 4b  Is it cell-autonomous?", "downstream_analysis/part4b_cell_tissue_concordance.R"),
    st("p5a", "Part 5a  What does the cell broadcast?", "downstream_analysis/part5a_secretome.R"),
    st("p5b", "Part 5b  Which TF drives it?", "downstream_analysis/part5b_tf_analysis.R"),
    st("p5c", "Part 5c  Is the KAT8/NSL program engaged?", "downstream_analysis/part5c_kat8_complex.R"))

  steps <- switch(target, "tissue" = TISSUE, "cells" = CELLS, "downstream" = DOWNSTREAM,
                  "all" = c(TISSUE, CELLS, DOWNSTREAM), NULL)
  if (is.null(steps)) {
    cat("\n[ERROR] Unknown target '", target, "'. Use: all | tissue | cells | downstream\n", sep = "")
    return(invisible(NULL))
  }

  ## ---- pre-flight -------------------------------------------------------
  cat("\n--- pre-flight ---\n")
  problems <- character(0)
  if (!file.exists(file.path(project_root, "parameters.R")))
    problems <- c(problems, "parameters.R  (are you in the project root?)")
  for (s in steps)
    if (!file.exists(file.path(project_root, s$path))) problems <- c(problems, s$path)
  if (length(problems) > 0) {
    cat("[ERROR] Missing required files:\n"); for (p in problems) cat("   - ", p, "\n", sep = "")
    cat("\n  Working directory:\n    ", project_root, "\n", sep = "")
    present <- list.files(project_root, pattern = "\\.R$")
    cat("  R scripts here (", length(present), "): ",
        if (length(present)) paste(present, collapse = ", ") else "(none)", "\n", sep = "")
    dsa <- file.path(project_root, "downstream_analysis")
    cat("  downstream_analysis/ exists: ", dir.exists(dsa), "\n", sep = "")
    cat("\n  Fix: setwd() to the folder holding parameters.R, and make sure the\n")
    cat("  whole repo was downloaded (including downstream_analysis/).\n")
    return(invisible(NULL))
  }
  cat("[OK] all ", length(steps), " scripts found\n", sep = "")

  ## ---- version check ----------------------------------------------------
  ## Every script carries a "## SCRIPT VERSION: YYYY-MM-DD" line in its header.
  ## They are all updated together, so they should all read the same date. A
  ## file with an older date is a leftover copy -- the exact situation that
  ## produced duplicate contrasts and stale figures before. This reports the
  ## dates BEFORE any analysis runs, so a bad copy is caught while it still
  ## costs nothing to fix.
  .script_version <- function(path) {
    ln <- tryCatch(readLines(path, n = 25, warn = FALSE), error = function(e) character(0))
    m  <- grep("^##\\s*SCRIPT VERSION:", ln, value = TRUE)
    if (!length(m)) return(NA_character_)
    trimws(sub("^##\\s*SCRIPT VERSION:", "", m[1]))
  }
  vfiles <- c("parameters.R", vapply(steps, function(s) s$path, character(1)))
  vers   <- vapply(file.path(project_root, vfiles), .script_version, character(1))
  names(vers) <- vfiles
  newest <- if (all(is.na(vers))) NA_character_ else max(vers, na.rm = TRUE)
  odd    <- names(vers)[is.na(vers) | vers != newest]
  if (length(odd) == 0) {
    cat("[OK] all scripts are version ", newest, "\n", sep = "")
  } else {
    cat("[WARN] script versions DISAGREE - newest is ", newest, "\n", sep = "")
    for (f in names(vers))
      cat(sprintf("   %-52s %s%s\n", f,
                  ifelse(is.na(vers[[f]]), "(no version stamp)", vers[[f]]),
                  ifelse(f %in% odd, "   <-- update this file", "")))
    cat("  Re-download the flagged file(s) before trusting the output.\n")
  }

  ## ---- pick the run folder: reuse newest, or make a fresh one ----------
  sp <- file.path(project_root, "savepoints")
  dir.create(sp, recursive = TRUE, showWarnings = FALSE)
  existing <- list.dirs(sp, recursive = FALSE, full.names = TRUE)
  existing <- existing[grepl("^RUN_", basename(existing))]
  existing <- existing[order(basename(existing), decreasing = TRUE)]  ## name = timestamp

  if (!fresh && length(existing) > 0) {
    KAT8_RUN_DIR <- existing[1]
    cat("Run folder   : ", KAT8_RUN_DIR, "  (REUSING)\n", sep = "")
  } else {
    KAT8_RUN_DIR <- file.path(sp, paste0("RUN_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(KAT8_RUN_DIR, recursive = TRUE, showWarnings = FALSE)
    cat("Run folder   : ", KAT8_RUN_DIR, "  (NEW)\n", sep = "")
  }
  assign("KAT8_RUN_DIR", KAT8_RUN_DIR, envir = globalenv())
  on.exit(suppressWarnings(rm("KAT8_RUN_DIR", envir = globalenv())), add = TRUE)
  done_dir <- file.path(KAT8_RUN_DIR, ".done")
  dir.create(done_dir, recursive = TRUE, showWarnings = FALSE)

  ## When Part 1 is SKIPPED, Parts 2/3/4 still need a pointer to the tables.
  ## Handle both layouts: the nested one (RUN_<tag>/tissue/tables) and the
  ## legacy flat one (RUN_<tag>/tables) used by older runs.
  ptr_target <- NULL
  if (dir.exists(file.path(KAT8_RUN_DIR, "tissue", "tables")))      ptr_target <- file.path(KAT8_RUN_DIR, "tissue")
  else if (dir.exists(file.path(KAT8_RUN_DIR, "tables")))           ptr_target <- KAT8_RUN_DIR
  if (!is.null(ptr_target)) {
    writeLines(normalizePath(ptr_target, winslash = "/", mustWork = FALSE),
               file.path(sp, "LATEST_RUN.txt"))
    cat("Tables from  : ", ptr_target, "\n", sep = "")
  }

  ## ---- completion test --------------------------------------------------
  is_done <- function(s) {
    if ("all" %in% force || s$key %in% force) return(FALSE)          ## forced re-run
    file.exists(file.path(done_dir, paste0(s$key, ".done")))
  }
  done_when <- function(s) {
    f <- file.path(done_dir, paste0(s$key, ".done"))
    if (!file.exists(f)) return("")
    paste0(" (", trimws(readLines(f, warn = FALSE))[1], ")")
  }
  mark_done <- function(s)
    writeLines(format(Sys.time()), file.path(done_dir, paste0(s$key, ".done")))

  ## ---- plan -------------------------------------------------------------
  todo <- vapply(steps, function(s) !is_done(s), logical(1))
  cat("\n==========================================================\n")
  cat("PLAN  (", sum(todo), " to run, ", sum(!todo), " already complete)\n", sep = "")
  cat("==========================================================\n")
  for (i in seq_along(steps))
    cat(sprintf("  %-6s %-46s %s\n", steps[[i]]$key, steps[[i]]$label,
                if (todo[i]) "RUN" else paste0("skip - done", done_when(steps[[i]]))))
  if (!any(todo)) {
    cat("\nNothing to do. To redo work:\n")
    cat("   RUN_FORCE <- \"all\";  source(\"run_all.R\")   # same folder\n")
    cat("   RUN_FRESH <- TRUE;    source(\"run_all.R\")   # brand-new folder\n")
    return(invisible(NULL))
  }

  ## ---- staging helper ---------------------------------------------------
  stage_downstream_inputs <- function() {
    dest <- file.path(project_root, "downstream_analysis", "data")
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    ## Search THIS RUN'S folder only.
    ##
    ## This used to search all of savepoints/ and take the "newest" file by
    ## sorting basenames in reverse. That is not a date sort. A file named
    ## DE_tissue_iWAT_KD_vs_CTL_TISSUE_PANELS_20260515_163857.csv sorts ABOVE
    ## DE_tissue_iWAT_KD_vs_CTL_20260727_152755.csv, because "T" beats "2"
    ## before either timestamp is reached -- so a brand-new run silently fed
    ## MONTHS-OLD DE tables to every downstream step.
    ##
    ## A run must be self-contained: the only correct source is the run folder
    ## the rest of this pipeline just wrote into. If a file is not there, that
    ## is a real failure and downstream is skipped rather than quietly reaching
    ## into another run.
    newest <- function(pat) {
      hits <- list.files(KAT8_RUN_DIR, pattern = pat, recursive = TRUE, full.names = TRUE)
      hits <- hits[!grepl("[/\\\\](data|_superseded_[0-9]+)[/\\\\]", hits)]
      if (length(hits) == 0) return(NA_character_)
      if (length(hits) > 1) {
        ## Several candidates within one run: take the most recently written.
        hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
        cat("[INFO] ", length(hits), " matches for ", pat, "; using ",
            basename(hits[1]), "\n", sep = "")
      }
      hits[1]
    }
    wanted <- list(
      c("^DE_cells_KAT8KD_vs_CTL_.*\\.csv$",   "DE_cells_KAT8KD_vs_CTL.csv"),
      c("^DE_tissue_iWAT_KD_vs_CTL_.*\\.csv$", "DE_tissue_iWAT_KD_vs_CTL.csv"),
      c("^DE_tissue_gWAT_KD_vs_CTL_.*\\.csv$", "DE_tissue_gWAT_KD_vs_CTL.csv"))
    cat("\n--- staging downstream inputs ---\n")
    ok <- TRUE
    for (w in wanted) {
      src <- newest(w[1])
      if (is.na(src)) {
        cat("[MISS] ", w[2], "  (not found in this run folder)\n", sep = "")
        ok <- FALSE; next
      }
      file.copy(src, file.path(dest, w[2]), overwrite = TRUE)
      rd <- file.path(KAT8_RUN_DIR, "data"); dir.create(rd, showWarnings = FALSE)
      file.copy(src, file.path(rd, w[2]), overwrite = TRUE)
      cat("[OK]   ", w[2], "  <- ", basename(src), "\n", sep = "")
    }
    if (!ok) cat("[WARN] some DE tables not found; downstream steps will be skipped.\n")
    ok
  }

  ## ---- run --------------------------------------------------------------
  is_down <- function(p) grepl("^downstream_analysis/", p)
  staged  <- if (identical(target, "downstream")) stage_downstream_inputs() else NA
  results <- data.frame(step = character(), status = character(),
                        seconds = numeric(), stringsAsFactors = FALSE)
  t_all <- Sys.time()

  for (i in seq_along(steps)) {
    s <- steps[[i]]
    if (!todo[i]) {
      results <- rbind(results, data.frame(step = s$label, status = "skipped (done)", seconds = 0))
      next
    }
    if (is_down(s$path) && is.na(staged)) staged <- stage_downstream_inputs()
    if (is_down(s$path) && isFALSE(staged)) {
      cat("\n[SKIP] ", s$label, " - inputs not staged\n", sep = "")
      results <- rbind(results, data.frame(step = s$label, status = "SKIPPED", seconds = 0)); next
    }

    cat("\n----------------------------------------------------------\n")
    cat(">>> ", s$label, "   [", s$path, "]\n", sep = "")
    cat("----------------------------------------------------------\n")
    t0 <- Sys.time()
    env <- new.env(parent = globalenv())   ## isolate each script
    ok <- tryCatch({ source(file.path(project_root, s$path), local = env, echo = FALSE); TRUE },
                   error = function(e) {
                     cat("\n[FAIL] ", s$label, "\n  ", conditionMessage(e), "\n", sep = ""); FALSE })
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    results <- rbind(results, data.frame(step = s$label,
                                         status = ifelse(ok, "OK", "FAILED"), seconds = el))
    if (!ok && isFALSE(s$fatal)) {
      cat("[WARN] ", s$label, " failed but is ADVISORY - continuing.\n", sep = "")
      next
    }
    if (!ok) {
      cat("\n==========================================================\n")
      cat("PIPELINE HALTED at: ", s$label, "\n", sep = "")
      cat("  Completed steps are recorded, so re-running resumes from HERE --\n")
      cat("  fix the error, then simply: source(\"run_all.R\")\n")
      cat("==========================================================\n")
      print(results, row.names = FALSE)
      return(invisible(results))   ## never quit(): that would kill an RStudio session
    }
    mark_done(s)
    cat("[OK] ", s$label, "  (", el, "s)\n", sep = "")
  }

  cat("\n==========================================================\n")
  cat("PIPELINE COMPLETE  (",
      round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min)\n", sep = "")
  cat("==========================================================\n")
  print(results, row.names = FALSE)
  cat("\nRun folder: ", KAT8_RUN_DIR, "\n", sep = "")
  cat("To redo work:  RUN_FORCE <- \"all\" (same folder)  |  RUN_FRESH <- TRUE (new folder)\n")
  cat("               one step only:  RUN_FORCE <- \"part4\"  then run_pipeline()\n")
  cat("               (these flags are one-shot -- they do not carry over)\n")
  cat("\nCheck the SANITY / VALIDATION blocks above before interpreting anything.\n")
  invisible(results)
}

## ---------------------------------------------------------
## Entry point -- safe to source(); returns instead of quitting.
## ---------------------------------------------------------
.res <- run_pipeline()
if (!interactive() && !is.null(.res) && any(.res$status == "FAILED")) quit(status = 1)
