#!/usr/bin/env Rscript

## =========================================================
## KAT8 - PART 8: TRANSCRIPTION FACTOR ACTIVITY INFERENCE
## =========================================================
## Question: KAT8 is an H4K16ac co-activator, not a TF. When it is
##   lost, WHICH transcription-factor programs shift? Doing this in
##   the 3T3-L1 cells (clean adipocyte system) gives the
##   cell-autonomous regulator response; doing it in iWAT/gWAT and
##   comparing tells us which TF shifts are cell-autonomous vs
##   likely driven by tissue composition (infiltration).
##
## Method: decoupleR univariate linear model (ULM) over the CollecTRI
##   mouse regulon network, using the genome-wide DESeq2 Wald
##   statistic as input. This aggregates the coordinated behaviour of
##   each TF's whole target set, so it does NOT require many DEGs.
##
## Inputs : the DE CSVs configured in config_downstream.R
## Outputs: per-contrast TF activity tables, a cells barplot, and a
##   TF x contrast activity heatmap (the cell-autonomy readout).
## =========================================================

## 1) Packages ----------------------------------------------
required_cran <- c("dplyr", "tidyr", "ggplot2", "tibble", "pheatmap", "RColorBrewer")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("decoupleR", quietly = TRUE)) BiocManager::install("decoupleR", ask = FALSE, update = FALSE)
## OmnipathR is used by decoupleR::get_collectri() to fetch the network
if (!requireNamespace("OmnipathR", quietly = TRUE)) BiocManager::install("OmnipathR", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(decoupleR)
})

## 2) Config -------------------------------------------------
source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)

## 3) CollecTRI network (cached) -----------------------------
## CollecTRI (Muller-Dott et al. 2023) is the current best signed
## TF -> target network: ~1000 TFs, each edge signed +1 (activation)
## or -1 (repression), including complexes (NF-kB, AP-1, SREBF...).
get_network <- function() {
  cache <- tf_params$network_cache
  if (file.exists(cache)) {
    cat("[INFO] Loading cached CollecTRI network:", cache, "\n")
    return(readRDS(cache))
  }
  cat("[INFO] Fetching CollecTRI (", ORGANISM, ") via OmniPath ...\n", sep = "")
  net <- decoupleR::get_collectri(organism = ORGANISM, split_complexes = FALSE)
  saveRDS(net, cache)
  cat("[INFO] Cached network to:", cache, "\n")
  net
}
net <- get_network()
cat("[INFO] Network: ", length(unique(net$source)), " TFs, ",
    nrow(net), " edges\n", sep = "")

## 4) Load DE tables ----------------------------------------
de <- load_all_de()

## 5) Score TF activity, one contrast at a time -------------
## Running per contrast (rather than one shared matrix) lets each
## contrast keep its full expressed-gene universe and avoids NA
## handling across data sets with different gene coverage.
score_contrast <- function(de_tbl, label) {
  v <- de_tbl %>%
    dplyr::filter(!is.na(.data[[RANK_STAT]]), !is.na(gene), gene != "") %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  mat <- matrix(v[[RANK_STAT]], ncol = 1,
                dimnames = list(v$gene, label))

  if (identical(tf_params$method, "consensus")) {
    acts <- decoupleR::decouple(
      mat, net,
      .source = "source", .target = "target",
      statistics = c("ulm", "mlm", "wsum"),
      consensus_only = TRUE, minsize = tf_params$min_regulon_size
    )
  } else {
    acts <- decoupleR::run_ulm(
      mat, net,
      .source = "source", .target = "target", .mor = "mor",
      minsize = tf_params$min_regulon_size
    )
  }

  acts %>%
    dplyr::transmute(
      contrast = label,
      TF       = source,
      activity = score,
      p_value  = p_value
    ) %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(abs(activity)))
}

acts_all <- dplyr::bind_rows(lapply(names(de), function(nm) score_contrast(de[[nm]], nm)))

## 5b) Optional permutation null (calibration check) --------
## Shuffle gene labels and re-score to confirm the analytic p-values
## are not inflated at n = 4/group. Reports an empirical FDR per TF
## for the cells contrast only (the one that matters most here).
if (tf_params$n_permutations > 0) {
  cat("[INFO] Permutation calibration on cells (",
      tf_params$n_permutations, " shuffles)...\n", sep = "")
  cell_tbl <- de[[CELL_KEY]] %>%
    dplyr::filter(!is.na(.data[[RANK_STAT]])) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
  obs <- acts_all %>% dplyr::filter(contrast == CELL_KEY)
  null_max <- numeric(tf_params$n_permutations)
  for (i in seq_len(tf_params$n_permutations)) {
    shuffled <- cell_tbl
    shuffled$gene <- sample(shuffled$gene)
    m <- matrix(shuffled[[RANK_STAT]], ncol = 1,
                dimnames = list(shuffled$gene, "perm"))
    a <- decoupleR::run_ulm(m, net, .source = "source", .target = "target",
                            .mor = "mor", minsize = tf_params$min_regulon_size)
    null_max[i] <- max(abs(a$score), na.rm = TRUE)
  }
  obs$emp_p <- vapply(abs(obs$activity),
                      function(x) (1 + sum(null_max >= x)) / (1 + tf_params$n_permutations),
                      numeric(1))
  write.csv(obs, file.path(OUTDIR, paste0("TF_cells_permutation_calibration_", RUN_TAG, ".csv")),
            row.names = FALSE)
  cat("[INFO] Wrote permutation calibration for cells.\n")
}

## 6) Write per-contrast tables -----------------------------
for (nm in unique(acts_all$contrast)) {
  tab <- acts_all %>% dplyr::filter(contrast == nm) %>% dplyr::arrange(FDR, dplyr::desc(abs(activity)))
  write.csv(tab, file.path(OUTDIR, paste0("TF_activity_", nm, "_", RUN_TAG, ".csv")),
            row.names = FALSE)
  n_sig <- sum(tab$FDR < tf_params$fdr_cutoff, na.rm = TRUE)
  cat(sprintf("[RESULT] %-6s: %d TFs scored, %d significant (FDR < %.2f)\n",
              nm, nrow(tab), n_sig, tf_params$fdr_cutoff))
}

## 7) Cells barplot: top activated / repressed TFs ----------
cells_top <- acts_all %>%
  dplyr::filter(contrast == CELL_KEY) %>%
  dplyr::arrange(dplyr::desc(abs(activity))) %>%
  dplyr::slice_head(n = tf_params$top_n_cells) %>%
  dplyr::mutate(
    TF = factor(TF, levels = TF[order(activity)]),
    sig = ifelse(FDR < tf_params$fdr_cutoff, "FDR < 0.05", "n.s.")
  )

p_cells <- ggplot(cells_top, aes(x = activity, y = TF, fill = activity)) +
  geom_col(aes(alpha = sig), color = "grey20", linewidth = 0.2) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.35), name = NULL) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 0, name = "TF activity") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  labs(title = "TF activity in KAT8-KD 3T3-L1 adipocytes (cell-autonomous)",
       subtitle = "decoupleR ULM x CollecTRI, ranked on DESeq2 Wald statistic",
       x = "Inferred TF activity (KAT8KD vs CTL)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(OUTDIR, paste0("TF_activity_cells_barplot_", RUN_TAG, ".png")),
       p_cells, width = 7.5, height = 8, dpi = 300)
cat("[INFO] Wrote cells TF barplot.\n")

## 8) Cross-contrast heatmap: the cell-autonomy readout -----
## A TF active in cells AND tissue = cell-autonomous regulator.
## A TF active only in tissue (e.g. inflammatory NF-kB/IRF/STAT)
## = candidate non-cell-autonomous (infiltration-driven).
wide <- acts_all %>%
  dplyr::select(TF, contrast, activity) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = activity)

## choose TFs that are strong in ANY contrast
strength <- acts_all %>%
  dplyr::group_by(TF) %>%
  dplyr::summarise(m = max(abs(activity), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(m)) %>%
  dplyr::slice_head(n = tf_params$top_n_heatmap)

hm_mat <- wide %>%
  dplyr::filter(TF %in% strength$TF) %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
## order columns cells first
col_order <- intersect(c(CELL_KEY, TISSUE_KEYS), colnames(hm_mat))
hm_mat <- hm_mat[, col_order, drop = FALSE]
## A TF may be unscored in a contrast (regulon < minsize there); treat as
## no inferred activity (0) so row clustering has no NA distances.
hm_mat[is.na(hm_mat)] <- 0

lim <- max(abs(hm_mat), na.rm = TRUE)
brk <- seq(-lim, lim, length.out = 101)
png(file.path(OUTDIR, paste0("TF_activity_heatmap_", RUN_TAG, ".png")),
    width = 1400, height = 2000, res = 220)
pheatmap(hm_mat,
         color = colorRampPalette(c("#2166AC", "#4393C3", "grey97", "#D6604D", "#B2182B"))(100),
         breaks = brk,
         cluster_cols = FALSE, cluster_rows = TRUE,
         display_numbers = TRUE, number_format = "%.1f", fontsize_number = 7,
         main = "TF activity across contrasts (cells = cell-autonomous)",
         angle_col = 0)
dev.off()
cat("[INFO] Wrote TF activity heatmap.\n")

## 9) Combined table ----------------------------------------
write.csv(acts_all %>% dplyr::arrange(TF, contrast),
          file.path(OUTDIR, paste0("TF_activity_ALL_contrasts_", RUN_TAG, ".csv")),
          row.names = FALSE)

cat("\n[DONE] Part 8 complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("       Interpretation guide:\n")
cat("       - Repressed adipogenic TFs in CELLS (Pparg/Cebpa/Srebf1/nuclear receptors)\n")
cat("         = cell-autonomous collapse of the adipocyte program (candidate DIRECT route).\n")
cat("       - Inflammatory TFs (Nfkb/Rel/Stat/Irf) active in TISSUE but not CELLS\n")
cat("         = candidate non-cell-autonomous (infiltration-driven, INDIRECT).\n")
