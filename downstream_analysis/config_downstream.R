## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   shared config, sanity helpers, FDR_CUT
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## KAT8 downstream analysis - SHARED CONFIG
## =========================================================
## Sourced by:
##   part4b_cell_tissue_concordance.R  is the effect cell-autonomous?
##   part5a_secretome.R                what does the KAT8-null cell broadcast?
##   part5b_tf_analysis.R              which TF drives it?
##   part5c_kat8_complex.R             is the KAT8/NSL program engaged?
##
## All four scripts consume the per-gene DESeq2 result tables that
## your existing pipeline already writes (they contain gene_name,
## log2FoldChange, stat, pvalue, padj). No raw counts required.
## =========================================================

## ---------------------------------------------------------
## 0) INHERIT THE PROJECT-WIDE PARAMETERS
## ---------------------------------------------------------
## Pull in parameters.R so these scripts honour the same GENES_OF_INTEREST list
## as the main pipeline. Sourced only if not already loaded, and silently
## skipped when these scripts are run somewhere without it (the marker list
## then falls back to Kat8 alone).
if (!exists("GENES_OF_INTEREST")) {
  for (.p in c("parameters.R", file.path("..", "parameters.R"))) {
    if (file.exists(.p)) { suppressWarnings(try(source(.p), silent = TRUE)); break }
  }
}

## ---------------------------------------------------------
## 1) INPUT DE TABLES  (<<< EDIT THESE PATHS >>>)
## ---------------------------------------------------------
## Point each to the DESeq2 CSV your pipeline produced. Columns
## expected (any order): gene_name, log2FoldChange, stat, pvalue, padj.
## The first (unnamed) column of your CSVs is the gene name and is
## read as the row name automatically.
## Paths are relative to the PROJECT ROOT (the folder holding parameters.R),
## because that is the working directory the scripts are sourced from.
## BUGFIX: these previously read "data/..." , which resolves to <root>/data/ --
## a folder that does not exist. The files live in downstream_analysis/data/.
DE_PATHS <- list(
  cells = "downstream_analysis/data/DE_cells_KAT8KD_vs_CTL.csv",
  iWAT  = "downstream_analysis/data/DE_tissue_iWAT_KD_vs_CTL.csv",
  gWAT  = "downstream_analysis/data/DE_tissue_gWAT_KD_vs_CTL.csv"
)

## Which contrast is the clean, cell-autonomous reference?
CELL_KEY   <- "cells"
TISSUE_KEYS <- c("iWAT", "gWAT")

## ---------------------------------------------------------
## 2) OUTPUT
## ---------------------------------------------------------
## When run_all.R drives the pipeline everything lands in ONE run folder
## (savepoints/RUN_<tag>/downstream). Run standalone, results go to
## downstream_analysis/results as before.
OUTDIR <- if (exists("KAT8_RUN_DIR", envir = globalenv())) {
  file.path(get("KAT8_RUN_DIR", envir = globalenv()), "downstream")
} else {
  file.path("downstream_analysis", "results")
}
RUN_TAG <- format(Sys.time(), "%Y%m%d_%H%M%S")

## ---------------------------------------------------------
## 3) ORGANISM / REPRODUCIBILITY
## ---------------------------------------------------------
ORGANISM   <- "mouse"
MASTER_SEED <- 12345

## ---------------------------------------------------------
## 4) STATISTIC USED FOR RANKING / INFERENCE
## ---------------------------------------------------------
## The DESeq2 Wald statistic ("stat" = log2FC / lfcSE) is the
## variance-aware signed statistic. It is the correct input for
## both TF activity inference and rank-based concordance, and it
## needs NO fold-change threshold. This is why a low DEG count
## does not limit these analyses.
RANK_STAT <- "stat"

## ---------------------------------------------------------
## 5) TF ACTIVITY PARAMETERS (part5b)
## ---------------------------------------------------------
tf_params <- list(
  ## Minimum number of a TF's target genes that must be present
  ## (and expressed) in a contrast for that TF to be scored.
  ## 5 is the decoupleR default; do not go below 5.
  min_regulon_size = 5,

  ## Primary method. "ulm" (univariate linear model) is the
  ## benchmarked-robust default in decoupleR. Set run_consensus = TRUE
  ## to ALSO compute the multi-method consensus (ulm + mlm + wsum) as a
  ## robustness cross-check (recommended for the low-signal cells).
  method        = "ulm",
  run_consensus = TRUE,

  ## FDR threshold for calling a TF's activity significant.
  fdr_cutoff = 0.05,

  ## How many top TFs (by |activity|) to show in the cells barplot
  ## and in the cross-contrast heatmap.
  top_n_cells   = 25,
  top_n_heatmap = 30,

  ## Cache the CollecTRI network here so runs are reproducible and
  ## do not depend on OmniPath being reachable every time.
  network_cache = file.path("downstream_analysis", "collectri_mouse.rds"),

  ## Signed-regulon source. Defaults to the BUNDLED DoRothEA (mouse,
  ## confidence A-C) CSV committed alongside this config -- fully offline, so
  ## it sidesteps the OmnipathR static-table bug entirely. 273 TFs, all with
  ## >= 5 targets, incl. Pparg/Cebpa/Srebf1/Nr1h3/Ppara and the inflammatory
  ## set. Columns: source, target, mor.
  ## Set to NULL to instead fetch live CollecTRI from OmniPath (once your
  ## OmnipathR is fixed and OmniPath's server is up).
  network_csv = file.path("downstream_analysis", "regulon_dorothea_mm_ABC.csv"),

  ## Empirical null: permute gene labels to confirm the analytic p-values are
  ## calibrated (worth having at n = 4/group). Set to 0 to SKIP it entirely --
  ## everything else in part5b is unaffected, only the emp_p column is dropped.
  ## Permutations are run in batched columns, so 1000 takes ~1 min rather than
  ## the 10-30 min a per-shuffle loop cost. 200-500 is plenty for a check.
  n_permutations = 1000,
  perm_chunk     = 100,   ## columns per decoupleR call (memory vs speed)

  ## --- Biological highlight sets (for the cells read-out) ---
  ## Adipogenic / lipogenic master regulators. The concordance result
  ## predicts these should be FLAT in the cells (their collapse in iWAT
  ## tissue is non-cell-autonomous), and are the key TFs to watch.
  tf_adipogenic = c("Pparg","Cebpa","Cebpd","Srebf1","Srebf2","Nr1h3","Nr1h2",
                    "Rxra","Foxo1","Klf15","Klf5","Ebf1","Mlxipl","Ppara",
                    "Esrra","Nr3c1","Gata2","Gata3","Stat5a","Creb1"),
  ## Inflammatory / stress regulators (expected active in TISSUE, driven
  ## partly by infiltrate; the cells test whether any are cell-autonomous).
  tf_inflammatory = c("Nfkb1","Rela","Rel","Relb","Stat1","Stat3","Irf1","Irf3",
                      "Irf7","Irf8","Jun","Junb","Fos","Fosb","Atf3","Cebpb",
                      "Spi1","Egr1","Egr2","Nfe2l2","Hif1a")
)

## ---------------------------------------------------------
## 6) CONCORDANCE PARAMETERS (part4b)
## ---------------------------------------------------------
concordance_params <- list(
  ## Significance used only to DEFINE the well-powered tissue gene
  ## sets that we then test for direction in the cells. Cells are
  ## never thresholded, so their low power does not bite here.
  tissue_fdr = 0.05,

  ## Optional |log2FC| floor for defining tissue signature gene sets
  ## (0 = FDR only, recommended).
  tissue_lfc = 0,

  ## Correlation method for the rank-concordance of Wald statistics.
  cor_method = "spearman",

  ## Genes always labelled on the log2FC quadrant plot (from your slides).
  ## Metabolic / adipocyte markers = direct-candidate program.
  marker_metabolic = c("Ucp1", "Cidea", "Ppara", "Pparg", "Adipoq",
                       "Acox1", "Pck1", "Acacb", "Cd36", "Plin1"),
  ## Inflammatory / immune = candidate non-cell-autonomous program.
  marker_immune = c("Ccl2", "Ccl8", "Cxcl2", "Cxcl10", "Il1b",
                    "Il12b", "Itgax", "Ctss", "Cd44", "Gdf15"),

  ## Always-labelled anchors. Inherits GENES_OF_INTEREST from parameters.R when
  ## that has been sourced (the usual case via run_all.R); falls back to Kat8
  ## alone if these scripts are run without the main pipeline's parameters.
  marker_anchor = if (exists("GENES_OF_INTEREST")) GENES_OF_INTEREST else c("Kat8"),

  plot_width  = 7.5,
  plot_height = 7.0
)

## ---------------------------------------------------------
## 7) SMALL SHARED HELPERS
## ---------------------------------------------------------

## Read one DESeq2 CSV into a standard tibble:
##   gene, stat, log2FC, pvalue, padj
load_de_table <- function(path) {
  ## Be forgiving about the working directory: try the configured path, then
  ## the same file under downstream_analysis/data/ and ./data/. This means the
  ## scripts work whether they are sourced from the project root or from
  ## inside downstream_analysis/.
  if (!file.exists(path)) {
    alt <- c(file.path("downstream_analysis", "data", basename(path)),
             file.path("data", basename(path)),
             file.path("..", "downstream_analysis", "data", basename(path)))
    hit <- alt[file.exists(alt)]
    if (length(hit) > 0) {
      path <- hit[1]
    } else {
      stop("DE table not found: ", path,
           "\n  Looked in:\n   - ", paste(alt, collapse = "\n   - "),
           "\n  Working directory: ", getwd(),
           "\n  -> put the three DE CSVs in downstream_analysis/data/ (run_all.R",
           " stages them automatically), or edit DE_PATHS in config_downstream.R",
           call. = FALSE)
    }
  }
  df <- utils::read.csv(path, row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)

  ## gene name: prefer explicit gene_name column, else row names
  gene <- if ("gene_name" %in% colnames(df)) as.character(df$gene_name) else rownames(df)

  ## tolerate either "stat"/"log2FoldChange" (DESeq2) or logFC/P.Value aliases
  get_col <- function(d, prefer) {
    for (nm in prefer) if (nm %in% colnames(d)) return(d[[nm]])
    return(rep(NA_real_, nrow(d)))
  }
  out <- data.frame(
    gene   = gene,
    stat   = get_col(df, c("stat")),
    log2FC = get_col(df, c("log2FoldChange", "logFC")),
    pvalue = get_col(df, c("pvalue", "P.Value")),
    padj   = get_col(df, c("padj", "adj.P.Val")),
    stringsAsFactors = FALSE
  )

  ## Drop make.unique suffixes (".1", ".2") only when the base symbol
  ## is not itself duplicated, so cross-dataset merging on symbol works.
  out$gene <- as.character(out$gene)

  ## collapse exact-duplicate gene rows by keeping the most significant
  out <- out[order(out$pvalue, na.last = TRUE), , drop = FALSE]
  out <- out[!duplicated(out$gene), , drop = FALSE]
  rownames(out) <- NULL
  out
}

## Load all configured DE tables into a named list of tibbles.
load_all_de <- function(paths = DE_PATHS) {
  lst <- lapply(paths, load_de_table)
  names(lst) <- names(paths)
  for (nm in names(lst)) {
    n_stat <- sum(!is.na(lst[[nm]]$stat))
    cat(sprintf("[INFO] %-6s: %d genes (%d with valid Wald stat)\n",
                nm, nrow(lst[[nm]]), n_stat))
  }
  lst
}

## ---------------------------------------------------------
## 8) SHARED LOOK + SANITY HELPERS
## ---------------------------------------------------------
## Kept here so every downstream script stays short and they all look and
## behave the same. Sourced by part4b / part5a / part5b / part5c.

FDR_CUT  <- 0.05
col_up   <- "#B2182B"   ## up in KAT8 KD/KO
col_down <- "#2166AC"   ## down in KAT8 KD/KO

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) %+replace%
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text    = ggplot2::element_text(color = "black"),
      axis.title   = ggplot2::element_text(color = "black", face = "bold"),
      plot.title   = ggplot2::element_text(face = "bold", hjust = 0),
      plot.caption = ggplot2::element_text(size = 7.5, color = "grey35", hjust = 0))
}

## Sanity ledger: each script calls init_sanity() once, log_sanity() as it goes,
## and write_sanity() at the end. Everything the script checked lands in one CSV
## so a run can be audited without re-reading the console.
.SANITY <- new.env(parent = emptyenv())
init_sanity <- function(script) {
  .SANITY$rows <- list(); .SANITY$script <- script
  cat("\n== ", script, " ==\n", sep = "")
}
log_sanity <- function(check, value, status = "INFO") {
  .SANITY$rows[[length(.SANITY$rows) + 1]] <- data.frame(
    script = .SANITY$script, check = check, value = as.character(value),
    status = status, stringsAsFactors = FALSE)
  cat(sprintf("  [%-4s] %-48s %s\n", status, check, value))
}
write_sanity <- function() {
  df <- do.call(rbind, .SANITY$rows)
  if (is.null(df)) return(invisible(NULL))
  f <- file.path(OUTDIR, paste0("SANITY_", .SANITY$script, "_", RUN_TAG, ".csv"))
  utils::write.csv(df, f, row.names = FALSE)
  nf <- sum(df$status == "FAIL"); nw <- sum(df$status == "WARN")
  cat("\n  ", nrow(df), " checks | ", nf, " FAIL | ", nw, " WARN  -> ", basename(f), "\n", sep = "")
  if (nf > 0) print(df[df$status == "FAIL", ], row.names = FALSE)
  invisible(df)
}

## Guard against the gene-name bug class (dplyr::filter dropping row names,
## which silently produced "1","2","3" identifiers and emptied enrichment).
## Every downstream script calls this on its inputs.
check_gene_names <- function(de_list) {
  for (nm in names(de_list)) {
    fr <- mean(grepl("^[0-9]+$", de_list[[nm]]$gene))
    log_sanity(paste0(nm, ": gene IDs that look numeric"),
               sprintf("%.1f%%", 100 * fr), ifelse(fr > 0.5, "FAIL", "OK"))
    if (fr > 0.5)
      stop("Gene names in '", nm, "' are numeric -- upstream lost gene symbols.",
           call. = FALSE)
  }
}

## GO gene sets from org.Mm.eg.db. Offline and independent of these data, so
## annotation is NOT circular (a seed list curated from this dataset could only
## ever re-find genes already annotated by eye). Cached per session.
.GOCACHE <- new.env(parent = emptyenv())
go_genes <- function(go_id) {
  if (!is.null(.GOCACHE[[go_id]])) return(.GOCACHE[[go_id]])
  out <- tryCatch({
    eg <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = go_id,
                                keytype = "GOALL", columns = "SYMBOL")
    unique(eg$SYMBOL[!is.na(eg$SYMBOL)])
  }, error = function(e) { message("  GO lookup failed for ", go_id, ": ",
                                   conditionMessage(e)); character(0) })
  .GOCACHE[[go_id]] <- out
  out
}
GO_IDS <- c(extracellular = "GO:0005576", ecm = "GO:0031012",
            cell_surface  = "GO:0009986", ligand = "GO:0048018",
            cytokine      = "GO:0005125", growth_factor = "GO:0008083",
            tf_activity   = "GO:0003700")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] Loaded config_downstream.R  | RUN_TAG =", RUN_TAG, "\n")
cat("[INFO] Output dir:", normalizePath(OUTDIR, mustWork = FALSE), "\n")
