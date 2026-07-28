#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   QC diagnostics A-F; contamination vs genotype within depot
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## KAT8 bulk RNA-seq - PART 1b: QC DIAGNOSTICS
## =========================================================
## Runs after Part 1, before the pathway analyses. Everything here is a check
## you should look at BEFORE trusting a result -- these are the diagnostics
## that catch a broken model or a bad sample, and they were the main gap in
## this pipeline's QC.
##
##   A. p-value histograms          is the DE model calibrated?
##   B. independent filtering       how many genes got padj = NA, and why
##   C. unstable estimates          genes whose log2FC is driven by ~1 sample
##   D. contamination screen        foreign tissue in a fat sample, PER SAMPLE
##   E. PCA scree + PC-covariate    what actually drives the main components
##   F. sample-sample correlation   outlier detection independent of PCA
##
## Reads Part 1's QC_data_*.rds and DE_tissue_*.csv from the current run folder.
## Writes plots + a QC_FLAGS_*.csv summarising everything that needs attention.
## =========================================================

required_pkgs <- c("dplyr", "ggplot2", "tidyr", "pheatmap", "RColorBrewer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
  library(pheatmap); library(RColorBrewer)
})

params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) source(params_file)

## ---- locate the run folder (same logic as Parts 2/3/4) ----
savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) stop("savepoints/ not found - run Part 1 first.", call. = FALSE)
outdir <- NULL
ptr <- file.path(savepoint_dir, "LATEST_RUN.txt")
if (file.exists(ptr)) {
  cand <- trimws(readLines(ptr, warn = FALSE)); cand <- cand[nzchar(cand)]
  if (length(cand) && dir.exists(cand[1])) outdir <- cand[1]
}
if (is.null(outdir)) {
  rd <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  rd <- rd[grepl("^RUN_", basename(rd))]
  if (!length(rd)) stop("No RUN_ folders found.", call. = FALSE)
  outdir <- rd[order(basename(rd), decreasing = TRUE)][1]
}
.clean_run_tag <- function(dir) {
  b <- basename(dir); if (b %in% c("tissue", "cells")) b <- basename(dirname(dir))
  b <- sub("^RUN_", "", b)
  m <- regmatches(b, regexpr("[0-9]{8}_[0-9]{6}$", b))
  if (length(m) == 1 && nzchar(m)) m else b
}
run_tag   <- .clean_run_tag(outdir)
tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
cat("[INFO] QC diagnostics for: ", outdir, "\n", sep = "")

FLAGS <- list()
## A diagnostic script must never take the pipeline down with it. Each section
## runs inside this wrapper: a failure is recorded as a WARN flag and the rest
## of the QC still runs.
safely <- function(label, expr) {
  tryCatch(force(expr), error = function(e) {
    FLAGS[[length(FLAGS) + 1]] <<- data.frame(check = paste0(label, " [section failed]"),
      value = conditionMessage(e), status = "WARN",
      action = "this QC section did not run; others still did",
      stringsAsFactors = FALSE)
    cat(sprintf("  [WARN] %-44s %s\n", paste0(label, " failed"), conditionMessage(e)))
    NULL })
}
flag <- function(check, value, status = "OK", action = "") {
  FLAGS[[length(FLAGS) + 1]] <<- data.frame(check = check, value = as.character(value),
                                            status = status, action = action,
                                            stringsAsFactors = FALSE)
  cat(sprintf("  [%-4s] %-44s %s\n", status, check, value))
  ## The REASONING is the product of this script, not the pass/fail letter.
  ## `action` used to be stored and then surfaced only in the end-of-run table,
  ## and that table lists only checks that need attention -- so a verdict like
  ## "INCONCLUSIVE: Alb is the only restricted marker above background, so
  ## 'do several move together' cannot be asked" was computed in full and then
  ## thrown away, precisely because the check passed. The passing checks are
  ## where the interpretation lives. Print it wherever there is one.
  if (nzchar(action))
    cat(strwrap(action, width = 96, prefix = "         "), sep = "\n")
}

## =========================================================
## A. p-VALUE HISTOGRAMS  -- the single best model diagnostic
## =========================================================
## A well-behaved contrast gives a FLAT histogram with a spike near 0.
##   flat + spike at 0      = healthy
##   hill in the middle     = test is mis-calibrated / model misspecified
##   spike at 1             = over-conservative filtering
## Anything other than the first pattern means the DE results are suspect,
## regardless of how good the biology looks.
cat("\n=== A. p-value histograms ===\n")
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
if (!length(de_files)) stop("No DE_tissue_*.csv found in ", tables_dir, call. = FALSE)

safely("A. p-value histograms", {
pv_list <- list()
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  pv <- d$pvalue[is.finite(d$pvalue)]
  pv_list[[cn]] <- data.frame(contrast = cn, pvalue = pv)

  ## Fraction in the flat region: for a calibrated test the density of p in
  ## [0.5,1] should be ~= the null density. A big hump in the middle shows up
  ## as the mid bins exceeding the high bins.
  h   <- hist(pv, breaks = seq(0, 1, 0.05), plot = FALSE)$counts
  mid <- mean(h[9:12]); high <- mean(h[17:20])
  ratio <- mid / max(high, 1)
  flag(paste0(cn, ": p-hist mid/high ratio"), sprintf("%.2f", ratio),
       ifelse(ratio > 1.5, "WARN", "OK"),
       ifelse(ratio > 1.5, "hump in the middle -> check model/covariates", ""))
  flag(paste0(cn, ": fraction p < 0.01"),
       sprintf("%.1f%%", 100 * mean(pv < 0.01)), "INFO")
}
pv_all <- dplyr::bind_rows(pv_list)
p <- ggplot(pv_all, aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.02), fill = "grey70", colour = "black", linewidth = .2) +
  facet_wrap(~ contrast, scales = "free_y") +
  labs(title = "p-value histograms (raw p, all tested genes)",
       subtitle = "healthy = flat with a spike at 0; a hump in the middle means the model is mis-specified",
       x = "raw p-value", y = "genes") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave(file.path(plots_dir, paste0("QC_pvalue_histograms_", run_tag, ".png")),
       p, width = 9, height = 4.2, dpi = 300, bg = "white")
})

## =========================================================
## B. INDEPENDENT FILTERING
## =========================================================
cat("\n=== B. independent filtering ===\n")
safely("B. independent filtering", {
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  n_na <- sum(is.na(d$padj))
  flag(paste0(cn, ": genes with padj = NA"),
       sprintf("%d / %d (%.1f%%)", n_na, nrow(d), 100 * n_na / nrow(d)),
       ifelse(n_na / nrow(d) > 0.5, "WARN", "OK"),
       ifelse(n_na / nrow(d) > 0.5, "most genes filtered out - check prefilter", ""))
}
})

## =========================================================
## C. UNSTABLE ESTIMATES  -- log2FC driven by very few samples
## =========================================================
## A huge |log2FC| paired with a huge lfcSE is not a big biological effect; it
## means the gene is ~zero in one group and nonzero in a couple of samples of
## the other. These dominate volcano plots and look impressive but are usually
## contamination or a single outlier library.
cat("\n=== C. unstable log2FC estimates ===\n")
safely("C. unstable log2FC estimates", {
unstable_all <- list()
for (f in de_files) {
  cn <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  d  <- read.csv(f, stringsAsFactors = FALSE)
  d$gene <- if ("gene_name" %in% names(d)) d$gene_name else d[[1]]
  u <- d %>% dplyr::filter(!is.na(padj), padj < 0.05, lfcSE > 2) %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, padj) %>%
    dplyr::mutate(contrast = cn)
  unstable_all[[cn]] <- u
  flag(paste0(cn, ": significant genes with lfcSE > 2"), nrow(u),
       ifelse(nrow(u) > 0, "WARN", "OK"),
       ifelse(nrow(u) > 0, "inspect per-sample counts before reporting these", ""))
  if (nrow(u)) cat("        ", paste(utils::head(u$gene, 10), collapse = ", "), "\n")
}
unstable <- dplyr::bind_rows(unstable_all)
if (nrow(unstable))
  write.csv(unstable, file.path(tables_dir, paste0("QC_unstable_estimates_", run_tag, ".csv")),
            row.names = FALSE)
})

## =========================================================
## D. CONTAMINATION SCREEN  -- per sample
## =========================================================
## Fat pads are dissected next to other tissues. A single sample carrying
## epididymis / testis / muscle can create "significant" DEGs on its own at
## n = 5 per group. This scores every sample against marker panels and, more
## importantly, tests whether contamination is CONFOUNDED WITH GENOTYPE --
## which is the case that actually biases results.
cat("\n=== D. contamination screen (per sample) ===\n")
qc_files <- list.files(tables_dir, pattern = "^QC_data_.*\\.rds$", full.names = TRUE)
if (!length(qc_files)) {
  flag("QC_data_*.rds", "NOT FOUND - contamination screen skipped", "WARN",
       "re-run Part 1 to regenerate")
} else {
  qc <- readRDS(qc_files[order(basename(qc_files), decreasing = TRUE)][1])
  ## TWO different matrices, used for different jobs:
  ##   vst_screen   pre-filter (~78k genes) -- best for CONTAMINATION, because a
  ##                marker for a foreign tissue may be removed by the prefilter.
  ##   vst_analysis post-filter -- the correct set for PCA/correlation. It also
  ##                excludes all-zero genes, which have ZERO VARIANCE and make
  ##                prcomp(scale. = TRUE) fail outright.
  vst_screen   <- qc$vst_mat_before; if (is.null(vst_screen))   vst_screen   <- qc$vst_mat_qc
  vst_analysis <- qc$vst_mat_qc;     if (is.null(vst_analysis)) vst_analysis <- vst_screen
  vst <- vst_analysis
  si  <- qc$sample_info
  si  <- si[match(colnames(vst), si$Sample), , drop = FALSE]
  flag("QC matrices", sprintf("screen %d genes / analysis %d genes",
                              nrow(vst_screen), nrow(vst_analysis)), "OK")

  ## The Liver panel is deliberately widened. With only Alb, Apoa1, Ttr and
  ## Serpina1a it was impossible to tell contamination from regulation: three of
  ## those four are also expressed in adipose, so a shift in them says nothing,
  ## and Alb moving alone is not enough to convict. Hepatocyte-restricted genes
  ## are added so the panel can answer the question it is asked.
  PANELS <- list(
    Sperm_testis = c("Prm1","Prm2","Tnp1","Tnp2","Smcp","Oaz3","Odf1","Akap4","Spata19","Izumo1"),
    Blood        = c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Alas2"),
    Muscle       = c("Acta1","Myh1","Myh2","Ckm","Tnnt3","Mylpf"),
    Liver        = c("Alb","Apoa1","Ttr","Serpina1a",
                     "Apob","Tf","Ahsg","F2","Cyp2e1","Fga","Fgb","Apoa2"),
    LymphNode    = c("Cd19","Ms4a1","Cd3e","Ccl21a","Ccl19"))

  ## Which members are RESTRICTED to the foreign tissue, and which are also
  ## expressed in fat. Real carry-over drags the restricted genes with it -- a
  ## piece of liver in the tube brings Alb AND Apob AND Fga. A shift confined to
  ## the shared genes is regulation happening in the adipocyte, and the two
  ## cannot be told apart from a panel mean. Anything not listed here is
  ## treated as shared, i.e. as weak evidence.
  EXCLUSIVE <- list(
    Sperm_testis = c("Prm1","Prm2","Tnp1","Tnp2","Smcp","Oaz3","Odf1","Akap4","Spata19","Izumo1"),
    Blood        = c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt","Alas2"),
    Muscle       = c("Acta1","Myh1","Myh2","Ckm","Tnnt3","Mylpf"),
    ## Cyp2e1 is deliberately NOT here. It sits at VST ~10.2 in these fat pads,
    ## the highest of the whole panel and above Alb -- it is a genuine adipose
    ## gene, not a hepatocyte-restricted one, and treating it as restricted
    ## produced a false CARRY-OVER call.
    Liver        = c("Alb","Apob","Tf","Ahsg","F2","Fga","Fgb","Apoa2"),
    LymphNode    = c("Cd19","Ms4a1","Cd3e"))

  ## ONE background threshold for the whole section. A marker sitting at the
  ## noise floor cannot testify either way, and both the panel-score test and
  ## the per-gene verdict need to agree on where that floor is -- two copies of
  ## "7.5" in two places is how the z > 3 / z > 5 contradiction happened.
  ## Most of the liver panel reads VST ~6 in fat, which is background.
  BACKGROUND_VST <- 7.5

  score_rows <- list()
  per_gene   <- list()
  panel_has_outlier <- list()
  safely("D. contamination screen", {
  for (nm in names(PANELS)) {
    g <- intersect(PANELS[[nm]], rownames(vst_screen))
    if (length(g) < 2) { flag(paste0(nm, " panel"), "too few genes measured", "INFO"); next }
    sc <- colMeans(vst_screen[g, colnames(vst), drop = FALSE])
    z  <- (sc - median(sc)) / (mad(sc) + 1e-9)          ## robust z, outlier-resistant
    score_rows[[nm]] <- data.frame(panel = nm, Sample = names(sc), score = as.numeric(sc),
                                   robust_z = as.numeric(z),
                                   Genotype = si$Genotype, Depot = si$Depot, Sex = si$Sex,
                                   stringsAsFactors = FALSE)
    hot <- names(sc)[z > 5]
    panel_has_outlier[[nm]] <- length(hot) > 0
    flag(paste0(nm, ": samples with robust z > 5"),
         ifelse(length(hot), paste(hot, collapse = ", "), "none"),
         ifelse(length(hot), "WARN", "OK"),
         ifelse(length(hot), "inspect these libraries", ""))

    ## The dangerous case: contamination tracking with genotype.
    ##
    ## This is tested TWICE, and the second test is the one that matters.
    ## Pooling iWAT and gWAT together confounds the comparison with depot: if
    ## one depot simply carries more of a marker, a pooled test can report a
    ## genotype effect that is really a depot effect (Simpson's paradox). DE is
    ## run within depot, so the within-depot test is what actually predicts
    ## whether contamination can generate false DEGs.
    wtest <- function(x, grp) {
      grp <- factor(grp)
      if (nlevels(grp) != 2 || any(table(grp) < 2)) return(NA_real_)
      tryCatch(suppressWarnings(
        wilcox.test(x[grp == levels(grp)[1]], x[grp == levels(grp)[2]])$p.value),
        error = function(e) NA_real_)
    }
    if (length(unique(si$Genotype)) == 2) {
      pv <- wtest(sc, si$Genotype)
      flag(paste0(nm, ": contamination vs GENOTYPE (pooled)"), sprintf("p = %.3g", pv),
           ifelse(!is.na(pv) && pv < 0.05, "INFO", "OK"),
           ifelse(!is.na(pv) && pv < 0.05,
                  "pooled across depots - check the within-depot test below", ""))

      ## Within each depot, which is how the DE model is actually fitted.
      ##
      ## A genotype association ALONE does not mean contamination. Dissection
      ## carry-over is sporadic: it lands in a few libraries and shows up as
      ## outlier samples. A panel that shifts smoothly with genotype while NO
      ## sample stands out is far more likely to be genuine regulation of those
      ## genes in adipose -- none of these panels is truly tissue-exclusive
      ## (Apoa1, Ttr and Serpina1a are all expressed in fat). So the WARN
      ## requires BOTH signals; without an outlier sample it is reported as a
      ## biological observation instead of a QC failure.
      ## ONE outlier threshold, stated. This used z > 3 while the panel-level
      ## flag above reported z > 5, so the same panel could print "samples with
      ## robust z > 5: none" and then "outlier samples present" three lines
      ## later. z > 3 is the right sensitivity for "is anything carrying extra
      ## signal", but the count has to be shown rather than asserted.
      n_warm <- sum(z > 3)
      any_outlier <- n_warm > 0
      for (dp in unique(si$Depot)) {
        k   <- si$Depot == dp
        pvd <- wtest(sc[k], si$Genotype[k])
        if (is.na(pvd)) next
        hit <- pvd < 0.05
        flag(sprintf("%s: marker score vs GENOTYPE within %s", nm, dp),
             sprintf("p = %.3g", pvd),
             if (hit && any_outlier) "WARN" else if (hit) "INFO" else "OK",
             if (hit && any_outlier)
               sprintf("CONFOUNDED with genotype AND %d sample(s) at z > 3 -> can create false DEGs", n_warm)
             else if (hit)
               "shifts with genotype but NO sample above z > 3 -> reads as regulation of these genes in fat, not carry-over"
             else "")

        ## WHICH GENE IS CARRYING IT?
        ##
        ## The test above is CIRCULAR for the question actually being asked.
        ## `sc` is a plain mean of raw VST across the panel, so it is dominated
        ## by whichever members are expressed highest -- for Liver in iWAT that
        ## is Alb (8.2) and Cyp2e1 (9.8), while the other nine sit at ~6 and
        ## contribute noise. Asking "is the liver score confounded with
        ## genotype?" when Alb is the gene under suspicion returns "yes,
        ## p = 0.002", which is just "Alb went down in the knockdown" restated.
        ## It is not independent evidence that liver tissue is in the tube.
        ##
        ## So ask the question that can actually separate the two. Drop the
        ## single strongest mover and re-test. Carry-over transfers a whole
        ## tissue, so the association survives losing any one member; the
        ## regulation of one gene does not survive losing that gene.
        if (hit) {
          smp    <- colnames(vst)[k]
          gt     <- factor(si$Genotype[k])
          mu     <- rowMeans(vst_screen[g, smp, drop = FALSE])
          g_test <- g[is.finite(mu) & mu >= BACKGROUND_VST]
          if (length(g_test) < 2) {
            flag(sprintf("%s in %s: is it one gene?", nm, dp),
                 sprintf("%d of %d above VST %.1f", length(g_test), length(g), BACKGROUND_VST),
                 "INFO",
                 sprintf(paste0("the panel mean here rests on %s -- every other member is at ",
                                "background, so the p-value above is that one gene restated ",
                                "rather than independent evidence of a whole tissue"),
                         if (length(g_test)) paste(g_test, collapse = ", ") else "nothing"))
          } else {
            d_gene <- vapply(g_test, function(gg) {
              m <- tapply(vst_screen[gg, smp], gt, mean); m[[2]] - m[[1]]
            }, numeric(1))
            drop_g <- g_test[which.max(abs(d_gene))]
            keep   <- setdiff(g_test, drop_g)
            p2     <- wtest(colMeans(vst_screen[keep, smp, drop = FALSE]), gt)
            surv   <- !is.na(p2) && p2 < 0.05
            ## Surviving is necessary but NOT sufficient. What survives has to
            ## be tissue-RESTRICTED to mean anything: the shared members are
            ## expressed in fat by definition, so a panel that keeps tracking
            ## genotype only through Cyp2e1 or Apoa1 is still describing
            ## adipose biology. Without this split the test would report
            ## "keep carry-over on the table" for Liver/iWAT on the strength of
            ## Cyp2e1 -- the highest-expressed gene of the panel in these pads,
            ## and the one already excluded from EXCLUSIVE for that reason.
            keep_restricted <- intersect(keep, EXCLUSIVE[[nm]])
            real <- surv && length(keep_restricted) > 0
            flag(sprintf("%s in %s: survives dropping %s?", nm, dp, drop_g),
                 sprintf("p = %.3g -> %.3g", pvd, p2),
                 if (real) "WARN" else "OK",
                 if (real)
                   sprintf(paste0("the association is NOT one gene: the tissue-restricted ",
                                  "member(s) %s still track genotype after %s is removed. ",
                                  "That is what a whole tissue looks like -- treat carry-over ",
                                  "as live for this panel."),
                           paste(keep_restricted, collapse = ", "), drop_g)
                 else if (surv)
                   sprintf(paste0("survives dropping %s, but only through %s -- shared ",
                                  "member(s) that fat expresses in its own right. Nothing ",
                                  "tissue-restricted is left carrying the signal, so this is ",
                                  "not carry-over evidence."),
                           drop_g, paste(keep, collapse = ", "))
                 else
                   sprintf(paste0("the whole association is %s. Remove that one gene and the ",
                                  "panel stops tracking genotype (p = %.3g). Carry-over cannot ",
                                  "show up in a single marker, so this is regulation of %s in ",
                                  "fat, not foreign tissue in the tube."), drop_g, p2, drop_g))
          }
        }
      }

      ## And the benign explanation, reported so the pooled result can be read
      ## correctly rather than guessed at.
      pvdep <- tryCatch(suppressWarnings(
        kruskal.test(sc ~ factor(si$Depot))$p.value), error = function(e) NA_real_)
      flag(paste0(nm, ": marker score vs DEPOT"), sprintf("p = %.3g", pvdep), "INFO",
           "a strong depot effect explains a pooled genotype signal")

      ## Per-gene breakdown for any panel that moved with genotype. A panel
      ## score is a mean, and a mean hides which gene produced it. True
      ## carry-over moves the WHOLE panel together, and above all it moves the
      ## tissue-EXCLUSIVE member (Alb for liver; Prm1/Prm2 for sperm). One or
      ## two shared genes moving while the exclusive marker sits still is
      ## regulation, not contamination.
      if (!is.na(pv) && pv < 0.05) {
        per_gene[[nm]] <- do.call(rbind, lapply(g, function(gg) {
          x <- vst_screen[gg, colnames(vst)]
          do.call(rbind, lapply(unique(si$Depot), function(dp) {
            k  <- si$Depot == dp
            gt <- factor(si$Genotype[k])
            if (nlevels(gt) != 2) return(NULL)
            m  <- tapply(x[k], gt, mean)
            data.frame(panel = nm, gene = gg, Depot = dp,
                       exclusive = gg %in% EXCLUSIVE[[nm]],
                       mean_vst = mean(x[k]),
                       mean_ctl = m[[1]], mean_kd = m[[2]],
                       delta_vst = m[[2]] - m[[1]],
                       p_genotype = wtest(x[k], si$Genotype[k]),
                       stringsAsFactors = FALSE)
          }))
        }))
      }
    }
  }
  pg <- dplyr::bind_rows(per_gene)
  if (nrow(pg)) {
    write.csv(pg, file.path(tables_dir, paste0("QC_contamination_per_gene_", run_tag, ".csv")),
              row.names = FALSE)
    cat("  per-gene breakdown for flagged panels -> QC_contamination_per_gene_",
        run_tag, ".csv\n", sep = "")
    print(pg[order(pg$panel, pg$Depot, pg$p_genotype), ], row.names = FALSE, digits = 3)

    ## THE VERDICT THE PANEL MEAN CANNOT GIVE.
    ## Carry-over moves the tissue-RESTRICTED members together. If only the
    ## shared members move, the genes are being regulated in the fat and there
    ## is nothing to correct for.
    for (nm in unique(pg$panel)) {
      for (dp in unique(pg$Depot[pg$panel == nm])) {
        q  <- pg[pg$panel == nm & pg$Depot == dp, , drop = FALSE]
        ## A marker sitting at the noise floor cannot testify either way.
        ## Most of the liver panel reads VST ~6 in fat, which is background;
        ## counting those as "did not move" is as wrong as counting them as
        ## evidence, so they are excluded from the tally and reported apart.
        floor_vst <- BACKGROUND_VST
        ex  <- q[q$exclusive & is.finite(q$mean_vst) & q$mean_vst >= floor_vst, , drop = FALSE]
        n_quiet <- sum(q$exclusive & (!is.finite(q$mean_vst) | q$mean_vst < floor_vst))
        if (!nrow(ex)) {
          flag(sprintf("%s in %s: restricted markers testable", nm, dp),
               sprintf("0 / %d above VST %.1f", n_quiet, floor_vst), "OK",
               "every tissue-restricted marker is at background - no carry-over detectable")
          next
        }
        sig <- ex[!is.na(ex$p_genotype) & ex$p_genotype < 0.05, , drop = FALSE]

        ## THREE conditions, because any one alone gives false positives:
        ##  (a) at least two restricted markers move
        ##  (b) they move in the SAME direction -- carry-over adds a whole
        ##      tissue, it cannot raise one marker and lower another
        ##  (c) some sample is actually contaminated. Carry-over is ADDITIVE and
        ##      sporadic: it shows up as outlier libraries. A panel that shifts
        ##      smoothly with genotype and has no outlier sample is regulation.
        same_dir <- nrow(sig) >= 2 && length(unique(sign(sig$delta_vst))) == 1
        hot_any  <- isTRUE(panel_has_outlier[[nm]])
        carry    <- nrow(sig) >= 2 && same_dir && hot_any
        verdict <- if (carry)
          "CARRY-OVER: restricted markers move together, in one direction, with outlier samples present"
        else if (nrow(sig) >= 2 && !same_dir)
          paste0("NOT carry-over: markers move in OPPOSITE directions (",
                 paste(sig$gene, sprintf("%+.2f", sig$delta_vst), collapse = ", "),
                 ") - contamination cannot do that")
        else if (nrow(sig) >= 2 && !hot_any)
          "NOT carry-over: markers move but NO sample is an outlier - carry-over is additive and sporadic, this is regulation"
        else if (nrow(sig) == 1 && nrow(ex) == 1)
          paste0("INCONCLUSIVE by this test: ", sig$gene[1],
                 " is the ONLY restricted marker expressed above background, so",
                 " 'do several move together' cannot be asked. Judge on the",
                 " outlier and direction evidence instead.")
        else if (nrow(sig) == 1)
          paste0("NOT carry-over: only ", sig$gene[1], " moves of ", nrow(ex),
                 " testable - a single gene, not a tissue")
        else
          "no restricted marker moves - regulation, not contamination"
        flag(sprintf("%s in %s: restricted markers moving", nm, dp),
             sprintf("%d / %d testable (%d at background, excluded)",
                     nrow(sig), nrow(ex), n_quiet),
             ifelse(carry, "WARN", "OK"), verdict)
      }
    }
  }

  cont <- dplyr::bind_rows(score_rows)
  if (nrow(cont)) {
    write.csv(cont, file.path(tables_dir, paste0("QC_contamination_scores_", run_tag, ".csv")),
              row.names = FALSE)
    p <- ggplot(cont, aes(x = Sample, y = score, fill = Genotype)) +
      geom_col(colour = "black", linewidth = .2) +
      facet_wrap(~ panel, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c(CTL = "#0072B2", KAT8KD = "#D55E00")) +
      labs(title = "Contamination screen: foreign-tissue marker score per sample",
           subtitle = "one sample standing far above the rest = likely dissection carry-over",
           x = NULL, y = "mean VST of panel") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 6),
            panel.grid.minor = element_blank())
    ggsave(file.path(plots_dir, paste0("QC_contamination_", run_tag, ".png")),
           p, width = 11, height = 12, dpi = 300, bg = "white")
  }
  })

  ## =========================================================
  ## E. PCA SCREE + WHAT DRIVES EACH COMPONENT
  ## =========================================================
  ## A PCA plot shows structure but not its CAUSE. This tests every PC against
  ## the known covariates, so "PC1 is depot" or "PC2 is library size" becomes a
  ## number rather than an impression.
  cat("\n=== E. PCA scree + PC-covariate association ===\n")
  ## idx (the top-500-variable genes) drives the sample-correlation heatmap in
  ## section F. It is built out here so a failure in E cannot starve F of it.
  rv  <- apply(vst, 1, stats::var)
  idx <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

  safely("E. PCA scree / PC-covariate", {
  ## The scree and the PC-covariate tests MUST describe the same PCA that
  ## Part 4 plots, otherwise "PC3 tracks genotype" refers to a component the
  ## figure does not show. So the convention comes from parameters.R rather
  ## than being hard-coded here.
  .ntop  <- if (!is.null(qc_plot_params$pca_ntop))  qc_plot_params$pca_ntop  else Inf
  .scale <- if (!is.null(qc_plot_params$pca_scale)) qc_plot_params$pca_scale else TRUE
  pmat <- vst
  rv_p <- rv
  if (isTRUE(.scale)) {                       ## scale. = TRUE cannot rescale a
    kv   <- is.finite(rv_p) & rv_p > 0        ## constant column
    pmat <- pmat[kv, , drop = FALSE]; rv_p <- rv_p[kv]
  }
  keep_n <- if (is.infinite(.ntop)) length(rv_p) else min(.ntop, length(rv_p))
  pidx   <- order(rv_p, decreasing = TRUE)[seq_len(keep_n)]
  cat("  PCA on ", keep_n, " genes, scale = ", .scale, "\n", sep = "")
  flag("PCA convention", sprintf("%d genes, scale = %s", keep_n, .scale), "OK",
       "matches the PCA drawn by Part 4")
  pca <- prcomp(t(pmat[pidx, , drop = FALSE]), scale. = .scale)
  ve  <- summary(pca)$importance[2, ] * 100
  scree <- data.frame(PC = factor(paste0("PC", 1:min(10, length(ve))),
                                  levels = paste0("PC", 1:min(10, length(ve)))),
                      variance = as.numeric(ve[1:min(10, length(ve))]))
  p <- ggplot(scree, aes(PC, variance)) +
    geom_col(fill = "grey70", colour = "black", linewidth = .2) +
    geom_text(aes(label = sprintf("%.1f%%", variance)), vjust = -0.4, size = 3) +
    labs(title = "PCA scree", x = NULL, y = "% variance explained") +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(plots_dir, paste0("QC_PCA_scree_", run_tag, ".png")),
         p, width = 7, height = 4, dpi = 300, bg = "white")

  libsize <- if (!is.null(qc$lib_sizes)) qc$lib_sizes[colnames(vst)] else NULL
  covs <- list(Genotype = si$Genotype, Depot = si$Depot, Sex = si$Sex)
  if (!is.null(libsize)) covs$LibrarySize <- as.numeric(libsize)
  for (pc in 1:min(4, ncol(pca$x))) {
    for (cvn in names(covs)) {
      cv <- covs[[cvn]]
      pv <- tryCatch({
        if (is.numeric(cv)) cor.test(pca$x[, pc], cv)$p.value
        else summary(aov(pca$x[, pc] ~ factor(cv)))[[1]][["Pr(>F)"]][1]
      }, error = function(e) NA_real_)
      if (!is.na(pv) && pv < 0.05)
        flag(sprintf("PC%d ~ %s", pc, cvn), sprintf("p = %.2g", pv), "INFO",
             "this covariate explains part of this component")
    }
  }
  if (is.null(libsize))
    flag("PC ~ LibrarySize", "lib_sizes not saved by Part 1", "INFO",
         "re-run Part 1 to enable this test")
  })

  ## =========================================================
  ## F. SAMPLE-SAMPLE CORRELATION  -- outliers PCA can miss
  ## =========================================================
  cat("\n=== F. sample-sample correlation ===\n")
  safely("F. sample-sample correlation", {
  cm <- cor(vst[idx, , drop = FALSE], method = "spearman")
  ann <- data.frame(Genotype = si$Genotype, Depot = si$Depot, row.names = si$Sample)
  png(file.path(plots_dir, paste0("QC_sample_correlation_", run_tag, ".png")),
      width = 1600, height = 1500, res = 180)
  pheatmap(cm, annotation_col = ann, annotation_row = ann,
           color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
           main = "Sample-sample Spearman correlation (top 500 variable genes)",
           fontsize = 7)
  dev.off()
  mean_cor <- (rowSums(cm) - 1) / (ncol(cm) - 1)
  z <- (mean_cor - median(mean_cor)) / (mad(mean_cor) + 1e-9)
  low <- names(mean_cor)[z < -5]

  ## A NAME ALONE IS NOT ACTIONABLE. "JS_06, JS_07, JS_09, JS_10" reads as four
  ## bad libraries until you look up what they are -- and if they all belong to
  ## ONE group, the far likelier explanation is that the group is genuinely
  ## different from everything else. Correlation here is against every other
  ## sample, most of which are controls or the other depot, so the arm with the
  ## largest treatment effect is EXPECTED to correlate least. That is biology,
  ## not a failed library. The annotation is printed so the difference is
  ## obvious without a manual lookup.
  if (length(low)) {
    li <- si[match(low, si$Sample), , drop = FALSE]
    grp <- paste(li$Depot, li$Genotype, li$Sex, sep = "_")
    tab <- table(grp)
    all_one <- length(tab) == 1
    flag("samples poorly correlated with the rest (z < -5)",
         paste(low, collapse = ", "),
         if (all_one) "INFO" else "WARN",
         if (all_one)
           paste0("all ", length(low), " are ", names(tab)[1],
                  " -- one group differing from the rest is a treatment effect, not a bad library")
         else "candidate outlier libraries - check against submission records")
    print(data.frame(Sample = low, Depot = li$Depot, Genotype = li$Genotype,
                     Sex = li$Sex, mean_cor = round(mean_cor[low], 3),
                     robust_z = round(z[low], 1)), row.names = FALSE)
  } else {
    flag("samples poorly correlated with the rest (z < -5)", "none", "OK")
  }
  })

  ## The PCA convention comparison that used to sit here has been REMOVED.
  ## The convention is now settled (all filtered genes, unit-scaled -- set in
  ## parameters.R) and drawing the rejected alternative beside it every run
  ## only invited second-guessing a decision that has been made.
}

## ---- summary ----
fl <- dplyr::bind_rows(FLAGS)
write.csv(fl, file.path(tables_dir, paste0("QC_FLAGS_", run_tag, ".csv")), row.names = FALSE)
nw <- sum(fl$status == "WARN")
cat("\n==========================================================\n")
cat("QC DIAGNOSTICS: ", nrow(fl), " checks | ", nw, " need attention\n", sep = "")
if (nw > 0) print(fl[fl$status == "WARN", c("check", "value", "action")], row.names = FALSE)
cat("==========================================================\n")
cat("[DONE] Part 1b. Plots in ", plots_dir, "\n", sep = "")
