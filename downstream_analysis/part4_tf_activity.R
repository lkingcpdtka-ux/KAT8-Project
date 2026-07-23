#!/usr/bin/env Rscript

## =========================================================
## KAT8 - PART 4: TRANSCRIPTION FACTOR ACTIVITY INFERENCE
##                (cells-focused; main analysis after DE/ORA/GSEA)
## =========================================================
## KAT8 is an H4K16ac co-activator, not a TF. When it is lost, WHICH
## transcription-factor programs shift? The clean read-out is the pure
## 3T3-L1 adipocytes (no immune infiltrate, no systemic context), so this
## script centres the CELLS and uses iWAT/gWAT tissue only as comparison.
##
## Method: decoupleR univariate linear model (ULM) over the signed CollecTRI
##   mouse regulon, scored on the GENOME-WIDE DESeq2 Wald statistic (every
##   tested gene, NOT the DEG list). ULM aggregates the coordinated shift of
##   each TF's whole target set, so a low cell DEG count does NOT limit it:
##   50 targets each drifting a little is a strong, testable signal even when
##   no single gene clears FDR.
##
## The concordance analysis (part4b) showed the adipocyte-identity/metabolic
## collapse seen in iWAT tissue is ABSENT in these cells. This script tests
## the regulator-level version of that: the adipogenic TFs (Pparg, Cebpa,
## Srebf1, ...) are expected to read strongly repressed in iWAT but FLAT in
## the cells -- i.e. that TF-program collapse is not cell-autonomous.
##
## Inputs : the DE CSVs configured in config_downstream.R (need gene_name + stat)
## Outputs: cells TF table + barplot + highlighted adipogenic/inflammatory
##   read-out; a TF x contrast activity heatmap; and a cells-vs-tissue TF
##   scatter (the cell-autonomy view).
##
## NOTE: fetching CollecTRI needs OmniPath reachable ONCE (then it is cached
##   to downstream_analysis/collectri_mouse.rds). Run this where you have
##   internet; offline machines can load a pre-saved cache.
## =========================================================

## 1) Packages ----------------------------------------------
required_cran <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "tibble", "pheatmap", "RColorBrewer")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("decoupleR", quietly = TRUE)) BiocManager::install("decoupleR", ask = FALSE, update = FALSE)
if (!requireNamespace("OmnipathR", quietly = TRUE)) BiocManager::install("OmnipathR", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(tibble); library(pheatmap); library(RColorBrewer); library(decoupleR)
})

## 2) Config -------------------------------------------------
source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
tp <- tf_params

## 3) CollecTRI network (cached) -----------------------------
## CollecTRI (Muller-Dott et al. 2023): the current best signed TF->target
## network. ~1000 mouse TFs; each edge signed +1 (activation) / -1 (repression),
## including complexes (NF-kB, AP-1, SREBF ...).
## Robust to the OmnipathR "static table" bug: try a local CSV, then
## get_collectri, then a cache-wipe retry, then a DoRothEA (A-C) fallback.
## The first success is cached to network_cache so later runs are offline.
valid_net <- function(x) {
  !is.null(x) && is.data.frame(x) && nrow(x) > 0 &&
    all(c("source", "target", "mor") %in% colnames(x))
}
get_network <- function() {
  if (file.exists(tp$network_cache)) {
    cat("[INFO] Loading cached network:", tp$network_cache, "\n")
    net <- readRDS(tp$network_cache); if (valid_net(net)) return(net)
    cat("[WARN] cached network looks invalid; refetching\n")
  }
  ## 1) local CSV escape hatch (set tf_params$network_csv)
  if (!is.null(tp$network_csv) && file.exists(tp$network_csv)) {
    cat("[INFO] Loading regulon from local CSV:", tp$network_csv, "\n")
    net <- utils::read.csv(tp$network_csv, stringsAsFactors = FALSE)
    if (!valid_net(net)) stop("network_csv must have columns: source, target, mor")
    saveRDS(net, tp$network_cache); return(net)
  }
  fetch_collectri <- function() decoupleR::get_collectri(organism = ORGANISM, split_complexes = FALSE)
  ## 2) CollecTRI via OmniPath
  cat("[INFO] Fetching CollecTRI (", ORGANISM, ") via OmniPath ...\n", sep = "")
  net <- tryCatch(fetch_collectri(), error = function(e) {
    cat("[WARN] get_collectri failed:", conditionMessage(e), "\n")
    ## 2b) OmnipathR often falls back to a buggy static table when its server
    ## is briefly down; wiping the cache and retrying usually hits the good path.
    cat("[INFO] wiping OmnipathR cache and retrying once...\n")
    try(OmnipathR::omnipath_cache_wipe(), silent = TRUE)
    tryCatch(fetch_collectri(), error = function(e2) {
      cat("[WARN] retry failed:", conditionMessage(e2), "\n"); NULL })
  })
  ## 3) DoRothEA (A-C) fallback (also OmniPath, different endpoint)
  if (!valid_net(net)) {
    cat("[INFO] falling back to DoRothEA (confidence A-C)...\n")
    net <- tryCatch(decoupleR::get_dorothea(organism = ORGANISM, levels = c("A","B","C")),
                    error = function(e) { cat("[WARN] DoRothEA failed:", conditionMessage(e), "\n"); NULL })
  }
  if (!valid_net(net)) {
    stop("Could not obtain a TF regulon.\n",
         "  Fixes: (1) BiocManager::install(c('OmnipathR','decoupleR')); ",
         "OmnipathR::omnipath_cache_wipe(); then retry.\n",
         "  Or (2) set tf_params$network_csv to a signed-regulon CSV ",
         "(columns: source, target, mor) in config_downstream.R.", call. = FALSE)
  }
  saveRDS(net, tp$network_cache)
  cat("[INFO] Cached network to:", tp$network_cache, "\n")
  net
}
net <- get_network()
cat("[INFO] Network:", length(unique(net$source)), "TFs,", nrow(net), "edges\n")

## 4) Load DE tables ----------------------------------------
de <- load_all_de()

## 5) Score TF activity per contrast ------------------------
## Run each contrast on its OWN expressed-gene universe (decoupleR intersects
## the network with the supplied genes and applies minsize automatically), so
## TFs are only scored where they actually have >= min_regulon_size targets
## measured in that data set.
score_contrast <- function(de_tbl, label) {
  v <- de_tbl %>%
    dplyr::filter(!is.na(.data[[RANK_STAT]]), !is.na(gene), gene != "") %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  mat <- matrix(v[[RANK_STAT]], ncol = 1, dimnames = list(v$gene, label))

  ulm <- decoupleR::run_ulm(mat, net, .source = "source", .target = "target",
                            .mor = "mor", minsize = tp$min_regulon_size) %>%
    dplyr::transmute(contrast = label, TF = source,
                     activity_ulm = score, p_ulm = p_value)

  out <- ulm
  if (isTRUE(tp$run_consensus)) {
    cons <- tryCatch(
      decoupleR::decouple(mat, net, .source = "source", .target = "target",
                          statistics = c("ulm", "mlm", "wsum"),
                          consensus_only = TRUE, minsize = tp$min_regulon_size) %>%
        dplyr::transmute(TF = source, activity_consensus = score),
      error = function(e) { cat("[WARN] consensus failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(cons)) out <- dplyr::left_join(out, cons, by = "TF")
  }
  out %>%
    dplyr::mutate(FDR = p.adjust(p_ulm, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(abs(activity_ulm)))
}

acts <- dplyr::bind_rows(lapply(names(de), function(nm) score_contrast(de[[nm]], nm)))

## 5b) Permutation calibration for the CELLS -----------------
## At n = 4/group, confirm the analytic p-values are not inflated: shuffle
## gene labels, re-score, and derive an empirical p from the null of max|score|.
if (tp$n_permutations > 0) {
  cat("[INFO] Permutation calibration on cells (", tp$n_permutations, " shuffles)...\n", sep = "")
  cv <- de[[CELL_KEY]] %>% dplyr::filter(!is.na(.data[[RANK_STAT]])) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  null_max <- numeric(tp$n_permutations)
  for (i in seq_len(tp$n_permutations)) {
    g <- sample(cv$gene)
    m <- matrix(cv[[RANK_STAT]], ncol = 1, dimnames = list(g, "perm"))
    a <- decoupleR::run_ulm(m, net, .source = "source", .target = "target",
                            .mor = "mor", minsize = tp$min_regulon_size)
    null_max[i] <- max(abs(a$score), na.rm = TRUE)
  }
  cell_idx <- acts$contrast == CELL_KEY
  acts$emp_p <- NA_real_
  acts$emp_p[cell_idx] <- vapply(abs(acts$activity_ulm[cell_idx]),
    function(x) (1 + sum(null_max >= x)) / (1 + tp$n_permutations), numeric(1))
}

## 6) Write tables ------------------------------------------
for (nm in unique(acts$contrast)) {
  tab <- acts %>% dplyr::filter(contrast == nm) %>% dplyr::arrange(FDR, dplyr::desc(abs(activity_ulm)))
  write.csv(tab, file.path(OUTDIR, paste0("TF_activity_", nm, "_", RUN_TAG, ".csv")), row.names = FALSE)
  cat(sprintf("[RESULT] %-6s: %d TFs scored, %d significant (FDR < %.2f)\n",
              nm, nrow(tab), sum(tab$FDR < tp$fdr_cutoff, na.rm = TRUE), tp$fdr_cutoff))
}
write.csv(acts, file.path(OUTDIR, paste0("TF_activity_ALL_", RUN_TAG, ".csv")), row.names = FALSE)

## 7) CELLS read-out: top activated / repressed TFs ---------
cells <- acts %>% dplyr::filter(contrast == CELL_KEY)
top <- cells %>% dplyr::arrange(dplyr::desc(abs(activity_ulm))) %>%
  dplyr::slice_head(n = tp$top_n_cells) %>%
  dplyr::mutate(TF = factor(TF, levels = TF[order(activity_ulm)]),
                sig = ifelse(FDR < tp$fdr_cutoff, "FDR < 0.05", "n.s."))
p_cells <- ggplot(top, aes(activity_ulm, TF, fill = activity_ulm)) +
  geom_col(aes(alpha = sig), color = "grey20", linewidth = 0.2) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.35), name = NULL) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B", midpoint = 0,
                       name = "TF activity") +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(title = "TF activity in KAT8-KD 3T3-L1 adipocytes (cell-autonomous)",
       subtitle = "decoupleR ULM x CollecTRI on the genome-wide DESeq2 Wald statistic",
       x = "Inferred TF activity (KAT8KD vs CTL)", y = NULL) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                   plot.title = element_text(face = "bold", size = 12))
ggsave(file.path(OUTDIR, paste0("TF_cells_barplot_", RUN_TAG, ".png")), p_cells,
       width = 7.5, height = 8, dpi = 300)

## 7b) The prediction: adipogenic vs inflammatory TFs -------
## Print where the key regulators land in CELLS. Prediction from part4b:
## adipogenic TFs ~ 0 (flat) in cells despite crashing in iWAT tissue.
hl <- function(set, name) {
  x <- cells %>% dplyr::filter(TF %in% set) %>% dplyr::arrange(dplyr::desc(activity_ulm)) %>%
    dplyr::select(TF, activity_ulm, FDR, dplyr::any_of("emp_p"))
  cat("\n--- ", name, " TFs in CELLS ---\n", sep = ""); if (nrow(x)) print(as.data.frame(x), row.names = FALSE)
  x
}
adip <- hl(tp$tf_adipogenic,  "Adipogenic/lipogenic")
infl <- hl(tp$tf_inflammatory, "Inflammatory/stress")

## 8) Cross-contrast heatmap (cells first) ------------------
wide <- acts %>% dplyr::select(TF, contrast, activity_ulm) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = activity_ulm)
strength <- acts %>% dplyr::group_by(TF) %>%
  dplyr::summarise(m = max(abs(activity_ulm), na.rm = TRUE), .groups = "drop") %>%
  dplyr::slice_max(m, n = tp$top_n_heatmap)
hm <- wide %>% dplyr::filter(TF %in% strength$TF) %>% tibble::column_to_rownames("TF") %>% as.matrix()
hm <- hm[, intersect(c(CELL_KEY, TISSUE_KEYS), colnames(hm)), drop = FALSE]
hm[is.na(hm)] <- 0
lim <- max(abs(hm), na.rm = TRUE)
png(file.path(OUTDIR, paste0("TF_heatmap_", RUN_TAG, ".png")), width = 1400, height = 2100, res = 220)
pheatmap(hm, color = colorRampPalette(c("#2166AC","#4393C3","grey97","#D6604D","#B2182B"))(100),
         breaks = seq(-lim, lim, length.out = 101), cluster_cols = FALSE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.1f", fontsize_number = 7, angle_col = 0,
         main = "TF activity across contrasts (cells = cell-autonomous)")
dev.off()

## 9) Cell-autonomy scatter: cells vs iWAT activity ---------
## TFs on the diagonal are cell-autonomous; TFs strong in iWAT but ~0 in
## cells (off the y=0 line, along the x axis) are non-cell-autonomous.
sc <- wide %>% dplyr::filter(!is.na(.data[[CELL_KEY]]), !is.na(.data[["iWAT"]])) %>%
  dplyr::mutate(class = dplyr::case_when(
    TF %in% tp$tf_adipogenic   ~ "adipogenic",
    TF %in% tp$tf_inflammatory ~ "inflammatory",
    TRUE ~ "other"))
p_sc <- ggplot(sc, aes(.data[["iWAT"]], .data[[CELL_KEY]])) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = .3) +
  geom_abline(slope = 1, linetype = 2, color = "grey60") +
  geom_point(aes(color = class), alpha = .6, size = 1.4) +
  ggrepel::geom_text_repel(data = dplyr::filter(sc, class != "other"),
                           aes(label = TF, color = class), size = 3, max.overlaps = 40) +
  scale_color_manual(values = c(adipogenic = "#2166AC", inflammatory = "#B2182B", other = "grey75")) +
  labs(title = "TF activity: cells vs iWAT (cell-autonomy view)",
       subtitle = "adipogenic TFs strong in iWAT but ~0 in cells => non-cell-autonomous",
       x = "iWAT TF activity", y = "cells TF activity") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                   plot.title = element_text(face = "bold"))
ggsave(file.path(OUTDIR, paste0("TF_cells_vs_iWAT_scatter_", RUN_TAG, ".png")), p_sc,
       width = 7.5, height = 6.5, dpi = 300)

cat("\n[DONE] Part 4 (TF activity) complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Read the cells table first; then check whether the adipogenic TFs above\n")
cat("are flat (predicted) vs their strong repression in the iWAT column of the heatmap.\n")
