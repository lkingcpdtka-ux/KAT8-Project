## =========================================================
## KAT8 downstream analysis - SHARED CONFIG
## =========================================================
## Sourced by:
##   part8_tf_activity.R           (TF activity inference)
##   part9_cell_tissue_concordance.R (direct/indirect, cell-autonomous)
##
## These two scripts consume the per-gene DESeq2 result tables that
## your existing pipeline already writes (they contain gene_name,
## log2FoldChange, stat, pvalue, padj). No raw counts required.
## =========================================================

## ---------------------------------------------------------
## 1) INPUT DE TABLES  (<<< EDIT THESE PATHS >>>)
## ---------------------------------------------------------
## Point each to the DESeq2 CSV your pipeline produced. Columns
## expected (any order): gene_name, log2FoldChange, stat, pvalue, padj.
## The first (unnamed) column of your CSVs is the gene name and is
## read as the row name automatically.
DE_PATHS <- list(
  cells = "data/DE_cells_KAT8KD_vs_CTL.csv",
  iWAT  = "data/DE_tissue_iWAT_KD_vs_CTL.csv",
  gWAT  = "data/DE_tissue_gWAT_KD_vs_CTL.csv"
)

## Which contrast is the clean, cell-autonomous reference?
CELL_KEY   <- "cells"
TISSUE_KEYS <- c("iWAT", "gWAT")

## ---------------------------------------------------------
## 2) OUTPUT
## ---------------------------------------------------------
OUTDIR <- file.path("downstream_analysis", "results")
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
## 5) TF ACTIVITY PARAMETERS (part8)
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

  ## Empirical null: permute gene labels to confirm the analytic
  ## p-values are calibrated (recommended at n = 4/group). 0 disables.
  n_permutations = 1000,

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
## 6) CONCORDANCE PARAMETERS (part9)
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
  ## Metabolic / adipocyte-identity = direct-candidate program.
  marker_metabolic = c("Ucp1", "Cidea", "Ppara", "Pparg", "Adipoq",
                       "Acox1", "Pck1", "Acacb", "Cd36", "Plin1"),
  ## Inflammatory / immune = candidate non-cell-autonomous program.
  marker_immune = c("Ccl2", "Ccl8", "Cxcl2", "Cxcl10", "Il1b",
                    "Il12b", "Itgax", "Ctss", "Cd44", "Gdf15"),

  ## Always keep Kat8 itself labelled.
  marker_anchor = c("Kat8"),

  plot_width  = 7.5,
  plot_height = 7.0
)

## ---------------------------------------------------------
## 7) SMALL SHARED HELPERS
## ---------------------------------------------------------

## Read one DESeq2 CSV into a standard tibble:
##   gene, stat, log2FC, pvalue, padj
load_de_table <- function(path) {
  if (!file.exists(path)) {
    stop("DE table not found: ", path,
         "\n  -> edit DE_PATHS in config_downstream.R", call. = FALSE)
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

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] Loaded config_downstream.R  | RUN_TAG =", RUN_TAG, "\n")
cat("[INFO] Output dir:", normalizePath(OUTDIR, mustWork = FALSE), "\n")
