#!/usr/bin/env Rscript

## =========================================================
## KAT8 - PART 5: CELL-AUTONOMOUS MECHANISM (3T3-L1)
## =========================================================
## THE QUESTION
##   The AKO mouse has big phenotypes (lipoatrophy, inflammation, fibrosis),
##   but part4b showed the dominant tissue signature is NOT cell-autonomous:
##   in KAT8-KD 3T3-L1 the adipocyte-identity / lipid / thermogenic programs
##   are FLAT with tight CIs (Pparg 95% CI [-0.17,+0.05]) despite a stronger
##   Kat8 knockdown than the tissue has. The published result -- KAT8 is
##   required only when lost BEFORE differentiation -- explains why: these
##   cells are a MAINTENANCE experiment, and mature adipocyte identity does
##   not need ongoing KAT8.
##
##   So: what DOES the KAT8-null adipocyte do, and which TF drives it?
##
## THE MODULES
##   A. SECRETOME     - what the KAT8-null adipocyte broadcasts, and whether
##                      it replicates in vivo. This is the candidate
##                      NON-CELL-AUTONOMOUS mechanism for the tissue phenotype.
##   B. TF layer 1    - which TF genes themselves change   (zero assumptions)
##   C. TF layer 2    - decoupleR regulon activity          (curated regulon)
##   D. TF layer 3    - MSigDB C3:TFT target enrichment     (independent, ChIP-derived)
##   E. PROGENy       - upstream signalling pathway activity
##   F. KAT8 COMPLEX  - MSL vs NSL output: nuclear-encoded OXPHOS with
##                      mito-encoded + ribosomal genes as specificity controls
##   G. CONVERGENCE   - which TFs are supported by MULTIPLE independent layers
##
## WHY THREE TF LAYERS: regulon activity (C) is a re-weighted score on the same
## statistics, so on its own it partly restates the DE result and depends on a
## curated, largely non-adipose regulon. Layer B uses no regulon at all and
## layer D uses independent ChIP-derived target sets. Agreement across layers
## -- not any single layer -- is the evidence.
##
## INPUT : the DE CSVs configured in config_downstream.R (gene_name, stat,
##         log2FoldChange, padj). No counts, no new data.
## OUTPUT: downstream_analysis/results/  (tables + plots + a run log)
## =========================================================

## ---------------------------------------------------------
## 0) Packages
## ---------------------------------------------------------
required_cran <- c("dplyr","tidyr","ggplot2","ggrepel","tibble","scales","RColorBrewer")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (bp in c("org.Mm.eg.db","AnnotationDbi","GO.db","fgsea","decoupleR")) {
  if (!requireNamespace(bp, quietly = TRUE)) {
    tryCatch(BiocManager::install(bp, ask = FALSE, update = FALSE),
             error = function(e) message("Could not install ", bp, ": ", conditionMessage(e)))
  }
}
## msigdbr (layer D) and progeny (module E) are optional; modules skip cleanly.
if (!requireNamespace("msigdbr", quietly = TRUE))
  tryCatch(install.packages("msigdbr"), error = function(e) NULL)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(tibble); library(scales); library(RColorBrewer)
  library(org.Mm.eg.db); library(AnnotationDbi); library(GO.db)
})

## ---------------------------------------------------------
## 1) Config
## ---------------------------------------------------------
source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)

FDR_CUT <- 0.05
col_up   <- "#B2182B"   ## up in KAT8 KD
col_down <- "#2166AC"   ## down in KAT8 KD
col_ns   <- "grey75"

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(panel.grid.minor = element_blank(),
          axis.text  = element_text(color = "black"),
          axis.title = element_text(color = "black", face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          plot.caption = element_text(size = 7.5, color = "grey35", hjust = 0))
}

## Sanity-check ledger: every module appends, written out at the end so the
## whole run can be audited from one file.
SANITY <- list()
log_sanity <- function(module, check, value, status = "INFO") {
  SANITY[[length(SANITY) + 1]] <<- data.frame(
    module = module, check = check, value = as.character(value),
    status = status, stringsAsFactors = FALSE)
  cat(sprintf("  [%-4s] %-46s %s\n", status, check, value))
}

## ---------------------------------------------------------
## 2) Load data + input validation
## ---------------------------------------------------------
cat("\n===================== INPUT VALIDATION =====================\n")
de <- load_all_de()
stopifnot(CELL_KEY %in% names(de))

## Guard: the gene-name bug class (dplyr::filter dropping rownames) produced
## "1","2","3" identifiers. If that ever recurs upstream, fail loudly HERE
## rather than silently returning empty enrichment downstream.
for (nm in names(de)) {
  g <- de[[nm]]$gene
  frac_numeric <- mean(grepl("^[0-9]+$", g))
  log_sanity("input", paste0(nm, ": genes that look numeric"),
             sprintf("%.1f%%", 100 * frac_numeric),
             ifelse(frac_numeric > 0.5, "FAIL", "OK"))
  if (frac_numeric > 0.5)
    stop("Gene names in '", nm, "' are numeric -- upstream lost gene symbols. ",
         "Check the DE export before continuing.", call. = FALSE)
  log_sanity("input", paste0(nm, ": genes with finite Wald stat"),
             sum(is.finite(de[[nm]]$stat)), "OK")
}
cells <- de[[CELL_KEY]]
tissues <- intersect(TISSUE_KEYS, names(de))

## Knockdown sanity: Kat8 must be down in the cells.
if ("Kat8" %in% cells$gene) {
  k <- cells[cells$gene == "Kat8", ]
  log_sanity("input", "Kat8 log2FC in cells (expect < 0)",
             sprintf("%.2f (padj=%.1e)", k$log2FC[1], k$padj[1]),
             ifelse(isTRUE(k$log2FC[1] < 0), "OK", "FAIL"))
}

cell_sig <- cells %>% dplyr::filter(!is.na(padj), padj < FDR_CUT)
log_sanity("input", "cell DEGs at FDR<0.05 (no logFC gate)", nrow(cell_sig), "OK")
log_sanity("input", "median |log2FC| among cell DEGs",
           sprintf("%.2f", median(abs(cell_sig$log2FC), na.rm = TRUE)), "INFO")

## ---------------------------------------------------------
## 3) GO-based annotation helper (offline, NON-circular)
## ---------------------------------------------------------
## Deliberately NOT a hand-curated list built from this dataset: that would be
## circular (it could only re-find genes already annotated by eye). GO terms
## come from org.Mm.eg.db, are independent of these data, and need no network.
genes_in_go <- function(go_ids) {
  out <- tryCatch({
    eg <- AnnotationDbi::select(org.Mm.eg.db, keys = go_ids, keytype = "GOALL",
                                columns = c("SYMBOL"))
    unique(eg$SYMBOL[!is.na(eg$SYMBOL)])
  }, error = function(e) { message("  GO lookup failed: ", conditionMessage(e)); character(0) })
  out
}

GO_SETS <- list(
  extracellular = c("GO:0005576"),   ## extracellular region
  ecm           = c("GO:0031012"),   ## extracellular matrix
  cell_surface  = c("GO:0009986"),   ## cell surface
  ligand        = c("GO:0048018"),   ## receptor ligand activity
  cytokine      = c("GO:0005125"),   ## cytokine activity
  growth_factor = c("GO:0008083"),   ## growth factor activity
  tf_activity   = c("GO:0003700")    ## DNA-binding transcription factor activity
)
cat("\n-- building GO annotation sets --\n")
GO_GENES <- lapply(GO_SETS, genes_in_go)
for (nm in names(GO_GENES)) log_sanity("annotation", paste0("GO set '", nm, "' size"),
                                       length(GO_GENES[[nm]]),
                                       ifelse(length(GO_GENES[[nm]]) > 50, "OK", "WARN"))

secreted_genes <- unique(c(GO_GENES$extracellular, GO_GENES$ecm))
ligand_genes   <- unique(c(GO_GENES$ligand, GO_GENES$cytokine, GO_GENES$growth_factor))
surface_genes  <- GO_GENES$cell_surface
tf_genes       <- GO_GENES$tf_activity

## =========================================================
## MODULE A: SECRETOME  -- what does the KAT8-null adipocyte broadcast?
## =========================================================
cat("\n===================== A. SECRETOME =====================\n")

secretome <- cell_sig %>%
  dplyr::mutate(
    class = dplyr::case_when(
      gene %in% ligand_genes   ~ "Secreted ligand",
      gene %in% secreted_genes ~ "Secreted / ECM",
      gene %in% surface_genes  ~ "Cell surface",
      TRUE ~ "Other"),
    direction = ifelse(log2FC > 0, "Up in KAT8 KD", "Down in KAT8 KD")) %>%
  dplyr::filter(class != "Other") %>%
  dplyr::arrange(dplyr::desc(log2FC))

log_sanity("secretome", "communication-class DEGs (of all cell DEGs)",
           paste0(nrow(secretome), " / ", nrow(cell_sig)), "OK")

## Replication in vivo: does each secreted DEG move the SAME way in tissue?
## This is the cheapest and strongest validation available without new data.
for (t in tissues) {
  tt <- de[[t]] %>% dplyr::select(gene, tis_lfc = log2FC, tis_padj = padj)
  secretome <- secretome %>% dplyr::left_join(tt, by = "gene")
  names(secretome)[names(secretome) == "tis_lfc"]  <- paste0(t, "_log2FC")
  names(secretome)[names(secretome) == "tis_padj"] <- paste0(t, "_padj")
}
## Replication score: how many tissues agree in direction (and significantly)
dir_cols <- paste0(tissues, "_log2FC"); pad_cols <- paste0(tissues, "_padj")
if (length(tissues) == 0 || nrow(secretome) == 0) {
  secretome$n_tissue_same_dir <- integer(nrow(secretome))
  secretome$n_tissue_sig_same_dir <- integer(nrow(secretome))
  log_sanity("secretome", "tissue tables available for replication test",
             "0 - replication skipped", "WARN")
} else {
  secretome$n_tissue_same_dir <- rowSums(
    sign(as.matrix(secretome[, dir_cols, drop = FALSE])) == sign(secretome$log2FC), na.rm = TRUE)
  secretome$n_tissue_sig_same_dir <- rowSums(
    (sign(as.matrix(secretome[, dir_cols, drop = FALSE])) == sign(secretome$log2FC)) &
      (as.matrix(secretome[, pad_cols, drop = FALSE]) < FDR_CUT), na.rm = TRUE)
  log_sanity("secretome", "depots used for replication test",
             paste(tissues, collapse = ", "), "OK")
}

write.csv(secretome, file.path(OUTDIR, paste0("A_secretome_annotated_", RUN_TAG, ".csv")), row.names = FALSE)

## Formal test: are UP-secreted cell genes enriched for in-vivo replication
## more than DOWN-secreted ones? (Model 1 vs Model 2, decided by the data.)
up_s   <- secretome %>% dplyr::filter(log2FC > 0)
down_s <- secretome %>% dplyr::filter(log2FC < 0)
rep_rate <- function(d) if (nrow(d) == 0) NA_real_ else mean(d$n_tissue_sig_same_dir > 0)
log_sanity("secretome", "UP-secreted replicating in >=1 depot (sig)",
           sprintf("%d/%d (%.0f%%)", sum(up_s$n_tissue_sig_same_dir > 0), nrow(up_s), 100*rep_rate(up_s)), "INFO")
log_sanity("secretome", "DOWN-secreted replicating in >=1 depot (sig)",
           sprintf("%d/%d (%.0f%%)", sum(down_s$n_tissue_sig_same_dir > 0), nrow(down_s), 100*rep_rate(down_s)), "INFO")
if (nrow(up_s) > 0 && nrow(down_s) > 0) {
  ct <- matrix(c(sum(up_s$n_tissue_sig_same_dir > 0),  sum(up_s$n_tissue_sig_same_dir == 0),
                 sum(down_s$n_tissue_sig_same_dir > 0), sum(down_s$n_tissue_sig_same_dir == 0)),
               nrow = 2, byrow = TRUE)
  ft <- fisher.test(ct)
  log_sanity("secretome", "Fisher: UP vs DOWN replication asymmetry",
             sprintf("OR=%.2f p=%.3g", ft$estimate, ft$p.value), "INFO")
}

## PLOT: the replicating secreted core, cells vs each depot
core <- secretome %>% dplyr::filter(n_tissue_sig_same_dir > 0) %>%
  dplyr::arrange(dplyr::desc(log2FC)) %>% dplyr::slice_head(n = 40)
if (nrow(core) >= 3) {
  plot_df <- core %>%
    dplyr::select(gene, cells = log2FC, dplyr::all_of(dir_cols)) %>%
    tidyr::pivot_longer(-gene, names_to = "dataset", values_to = "log2FC") %>%
    dplyr::mutate(dataset = sub("_log2FC$", "", dataset),
                  gene = factor(gene, levels = rev(core$gene)))
  p <- ggplot(plot_df, aes(log2FC, gene, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75,
             color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
    labs(x = expression(log[2]~fold~change~(KAT8~KD/KO~vs~CTL)), y = NULL,
         title = "Secreted / surface DEGs that replicate in vivo",
         caption = paste0("Cell DEGs (FDR<", FDR_CUT, ") in GO extracellular/ECM/ligand/surface classes that move the SAME direction ",
                          "in >=1 depot at FDR<", FDR_CUT, ".\nCandidate non-cell-autonomous mechanism: signals the KAT8-null adipocyte ",
                          "sends to its neighbours. n=", nrow(core), " shown.")) +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("A_secretome_replicating_", RUN_TAG, ".png")),
         p, width = 9, height = max(5, nrow(core) * 0.26), dpi = 300, bg = "white")
}

## =========================================================
## MODULE B: TF LAYER 1 -- TF genes that themselves change
## =========================================================
cat("\n============ B. TF LAYER 1 (TF genes; no regulon) ============\n")
tf_deg <- cell_sig %>%
  dplyr::filter(gene %in% tf_genes) %>%
  dplyr::arrange(dplyr::desc(log2FC)) %>%
  dplyr::select(gene, log2FC, stat, padj)
log_sanity("TF_layer1", "TF genes changed in cells (FDR<0.05)", nrow(tf_deg), "OK")
write.csv(tf_deg, file.path(OUTDIR, paste0("B_TF_genes_changed_", RUN_TAG, ".csv")), row.names = FALSE)

if (nrow(tf_deg) >= 2) {
  d <- tf_deg %>% dplyr::mutate(gene = factor(gene, levels = gene[order(log2FC)]))
  p <- ggplot(d, aes(log2FC, gene, fill = log2FC > 0)) +
    geom_col(color = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                      labels = c(`TRUE` = "Up in KAT8 KD", `FALSE` = "Down in KAT8 KD"), name = NULL) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = expression(log[2]~fold~change), y = NULL,
         title = "Layer 1: transcription-factor genes changed in KAT8-KD cells",
         caption = paste0("TFs = GO:0003700 (DNA-binding transcription factor activity). ",
                          "FDR<", FDR_CUT, ". No regulon or inference involved.")) +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("B_TF_genes_changed_", RUN_TAG, ".png")),
         p, width = 7.5, height = max(3.5, nrow(tf_deg) * 0.28), dpi = 300, bg = "white")
}

## =========================================================
## MODULE C: TF LAYER 2 -- decoupleR regulon activity
## =========================================================
cat("\n============ C. TF LAYER 2 (regulon activity) ============\n")
tf_act <- NULL
run_layer2 <- function() {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    log_sanity("TF_layer2", "decoupleR available", "NO - module skipped", "WARN"); return(NULL)
  }
  ## Reuse the same regulon resolution as part4 (bundled CSV -> OmniPath).
  net <- NULL
  if (!is.null(tf_params$network_csv) && file.exists(tf_params$network_csv)) {
    net <- utils::read.csv(tf_params$network_csv, stringsAsFactors = FALSE)
    log_sanity("TF_layer2", "regulon source", basename(tf_params$network_csv), "OK")
  } else if (file.exists(tf_params$network_cache)) {
    net <- readRDS(tf_params$network_cache); log_sanity("TF_layer2", "regulon source", "cached RDS", "OK")
  } else {
    net <- tryCatch(decoupleR::get_collectri(organism = ORGANISM, split_complexes = FALSE),
                    error = function(e) NULL)
    log_sanity("TF_layer2", "regulon source", ifelse(is.null(net), "UNAVAILABLE", "OmniPath CollecTRI"),
               ifelse(is.null(net), "WARN", "OK"))
  }
  if (is.null(net) || !all(c("source","target","mor") %in% colnames(net))) return(NULL)
  log_sanity("TF_layer2", "regulon TFs / edges",
             paste0(length(unique(net$source)), " / ", nrow(net)), "OK")

  v <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE)
  mat <- matrix(v$stat, ncol = 1, dimnames = list(v$gene, "cells"))
  acts <- tryCatch(decoupleR::run_ulm(mat, net, .source = "source", .target = "target",
                                      .mor = "mor", minsize = tf_params$min_regulon_size),
                   error = function(e) { message("run_ulm failed: ", conditionMessage(e)); NULL })
  if (is.null(acts)) return(NULL)

  acts <- acts %>%
    dplyr::transmute(TF = source, activity = score, p_value = p_value) %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(dplyr::desc(abs(activity)))

  ## Cross-annotate with the TF's OWN mRNA -> direct vs indirect logic.
  own <- cells %>% dplyr::select(gene, own_log2FC = log2FC, own_padj = padj)
  acts <- acts %>% dplyr::left_join(own, by = c("TF" = "gene")) %>%
    dplyr::mutate(
      own_mrna = dplyr::case_when(
        is.na(own_padj)                      ~ "not tested",
        own_padj < FDR_CUT & own_log2FC < 0      ~ "TF mRNA down",
        own_padj < FDR_CUT & own_log2FC > 0      ~ "TF mRNA up",
        TRUE                                 ~ "TF mRNA unchanged"),
      interpretation = dplyr::case_when(
        activity < 0 & own_mrna == "TF mRNA down"      ~ "direct-target candidate (TF expression lost)",
        activity > 0 & own_mrna == "TF mRNA up"        ~ "TF induced and active",
        own_mrna == "TF mRNA unchanged"                ~ "post-translational / indirect",
        TRUE                                            ~ "see own_mrna"))
  write.csv(acts, file.path(OUTDIR, paste0("C_TF_activity_cells_", RUN_TAG, ".csv")), row.names = FALSE)
  log_sanity("TF_layer2", "TFs scored", nrow(acts), "OK")
  log_sanity("TF_layer2", "TFs significant (FDR<0.05)", sum(acts$FDR < FDR_CUT, na.rm = TRUE), "OK")

  ## Leading edge: which targets drove each top call.
  ## BUGFIX vs the earlier version of this analysis: build a NAMED stat vector.
  ## Previously `stat_vec <- de$stat` was unnamed, so
  ## `stat_vec[match(target, names(stat_vec))]` matched against NULL, every
  ## contribution became NA, and the leading-edge table was silently EMPTY.
  stat_named <- setNames(v$stat, v$gene)
  top_tfs <- acts %>% dplyr::slice_head(n = 20) %>% dplyr::pull(TF)
  le <- net %>%
    dplyr::filter(source %in% top_tfs) %>%
    dplyr::mutate(target_stat = unname(stat_named[target]),
                  contribution = mor * target_stat) %>%
    dplyr::filter(!is.na(contribution)) %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(dplyr::desc(abs(contribution)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 15) %>% dplyr::ungroup()
  write.csv(le, file.path(OUTDIR, paste0("C_TF_leading_edge_", RUN_TAG, ".csv")), row.names = FALSE)
  log_sanity("TF_layer2", "leading-edge rows (0 would indicate the old bug)",
             nrow(le), ifelse(nrow(le) > 0, "OK", "FAIL"))

  top <- acts %>% dplyr::slice_head(n = tf_params$top_n_cells) %>%
    dplyr::mutate(TF = factor(TF, levels = TF[order(activity)]),
                  sig = ifelse(FDR < FDR_CUT, "FDR < 0.05", "n.s."))
  p <- ggplot(top, aes(activity, TF, fill = activity)) +
    geom_col(aes(alpha = sig), color = "black", linewidth = 0.25) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.4), name = NULL) +
    scale_fill_gradient2(low = col_down, mid = "grey95", high = col_up, midpoint = 0,
                         name = "activity") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = "Inferred TF activity (decoupleR ULM)", y = NULL,
         title = "Layer 2: TF regulon activity in KAT8-KD cells",
         caption = paste0("Signed regulon x genome-wide DESeq2 Wald statistic. NOTE: this is a re-weighted score on the ",
                          "same statistics,\nso treat it as a summary -- the evidence is agreement with layers 1 and 3. ",
                          "See C_TF_leading_edge for per-target audit.")) +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("C_TF_activity_cells_", RUN_TAG, ".png")),
         p, width = 8, height = 7.5, dpi = 300, bg = "white")
  acts
}
tf_act <- tryCatch(run_layer2(), error = function(e) { message("Layer 2 failed: ", conditionMessage(e)); NULL })

## =========================================================
## MODULE D: TF LAYER 3 -- independent ChIP-derived target sets
## =========================================================
cat("\n====== D. TF LAYER 3 (MSigDB C3:TFT, independent) ======\n")
tft_res <- NULL
run_layer3 <- function() {
  if (!requireNamespace("msigdbr", quietly = TRUE) || !requireNamespace("fgsea", quietly = TRUE)) {
    log_sanity("TF_layer3", "msigdbr/fgsea available", "NO - module skipped", "WARN"); return(NULL)
  }
  ## C3:TFT = transcription-factor target sets (GTRD, ChIP-seq derived).
  ## Independent of the DoRothEA/CollecTRI regulon used in layer 2, which is
  ## the entire point: two different evidence bases, one conclusion.
  ## msigdbr >= 10 renamed category/subcategory -> collection/subcollection.
  ## Try both APIs, then fall back to the whole C3 collection, so a package
  ## upgrade cannot silently disable this layer.
  try_msig <- function(...) tryCatch(msigdbr::msigdbr(...), error = function(e) NULL)
  msig <- try_msig(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD")
  if (is.null(msig) || nrow(msig) == 0)
    msig <- try_msig(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD")
  if (is.null(msig) || nrow(msig) == 0)
    msig <- try_msig(species = "Mus musculus", collection = "C3")
  if (is.null(msig) || nrow(msig) == 0)
    msig <- try_msig(species = "Mus musculus", category = "C3")
  if (!is.null(msig) && nrow(msig) > 0 && "gs_name" %in% colnames(msig)) {
    ## keep only TFT sets if we had to fall back to all of C3
    if (any(grepl("^TFT|GTRD", msig$gs_name))) {
      keep <- grepl("^TFT|GTRD", msig$gs_name)
      if (sum(keep) > 0) msig <- msig[keep, , drop = FALSE]
    }
  }
  if (is.null(msig) || nrow(msig) == 0) {
    log_sanity("TF_layer3", "MSigDB C3 retrieved", "NO", "WARN"); return(NULL)
  }
  sym_col <- if ("gene_symbol" %in% colnames(msig)) "gene_symbol" else "human_gene_symbol"
  sets <- split(msig[[sym_col]], msig$gs_name)
  sets <- sets[vapply(sets, function(x) length(unique(x)) >= 10, logical(1))]
  log_sanity("TF_layer3", "TFT gene sets tested", length(sets), "OK")

  v <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(stat))
  ranks <- setNames(v$stat, v$gene)
  res <- tryCatch(fgsea::fgsea(sets, ranks, minSize = 10, maxSize = 2000, eps = 0),
                  error = function(e) { message("  fgsea failed: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  res <- res %>% as.data.frame() %>% dplyr::arrange(padj) %>%
    dplyr::mutate(leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ","), character(1)))
  write.csv(res, file.path(OUTDIR, paste0("D_TFT_enrichment_cells_", RUN_TAG, ".csv")), row.names = FALSE)
  log_sanity("TF_layer3", "TFT sets significant (FDR<0.05)", sum(res$padj < FDR_CUT, na.rm = TRUE), "OK")

  top <- res %>% dplyr::filter(padj < FDR_CUT) %>% dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = 25)
  if (nrow(top) >= 2) {
    top <- top %>% dplyr::mutate(pathway = factor(pathway, levels = pathway[order(NES)]))
    p <- ggplot(top, aes(NES, pathway, fill = NES > 0)) +
      geom_col(color = "black", linewidth = 0.25) +
      scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                        labels = c(`TRUE` = "targets up", `FALSE` = "targets down"), name = NULL) +
      geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
      labs(x = "NES", y = NULL,
           title = "Layer 3: TF target-set enrichment (MSigDB C3:TFT, ChIP-derived)",
           caption = paste0("Independent of the layer-2 regulon. FDR<", FDR_CUT, ". Agreement with layers 1-2 is the evidence.")) +
      theme_pub() + theme(axis.text.y = element_text(size = 6.5))
    ggsave(file.path(OUTDIR, paste0("D_TFT_enrichment_cells_", RUN_TAG, ".png")),
           p, width = 9.5, height = max(5, nrow(top) * 0.28), dpi = 300, bg = "white")
  }
  res
}
tft_res <- tryCatch(run_layer3(), error = function(e) { message("Layer 3 failed: ", conditionMessage(e)); NULL })

## =========================================================
## MODULE E: PROGENy -- upstream signalling pathway activity
## =========================================================
cat("\n============== E. PROGENy (pathway activity) ==============\n")
run_progeny <- function() {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    log_sanity("progeny", "decoupleR available", "NO - skipped", "WARN"); return(NULL) }
  net <- tryCatch(decoupleR::get_progeny(organism = ORGANISM, top = 500), error = function(e) NULL)
  if ((is.null(net) || nrow(net) < 50) && requireNamespace("progeny", quietly = TRUE)) {
    net <- tryCatch({
      m <- progeny::getModel(organism = "Mouse", top = 500)
      m$gene <- rownames(m)
      tidyr::pivot_longer(m, cols = -gene, names_to = "source", values_to = "weight") %>%
        dplyr::filter(weight != 0) %>% dplyr::transmute(source, target = gene, weight)
    }, error = function(e) NULL)
  }
  if (is.null(net) || nrow(net) < 50) {
    log_sanity("progeny", "PROGENy model retrieved", "NO - skipped (needs OmniPath or progeny pkg)", "WARN")
    return(NULL)
  }
  log_sanity("progeny", "PROGENy edges", nrow(net), "OK")
  v <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE)
  mat <- matrix(v$stat, ncol = 1, dimnames = list(v$gene, "cells"))
  acts <- tryCatch(decoupleR::run_mlm(mat, net, .source = "source", .target = "target",
                                      .mor = "weight", minsize = 5),
                   error = function(e) NULL)
  if (is.null(acts)) return(NULL)
  acts <- acts %>% dplyr::transmute(pathway = source, score, p_value) %>%
    dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>% dplyr::arrange(dplyr::desc(score))
  write.csv(acts, file.path(OUTDIR, paste0("E_PROGENy_cells_", RUN_TAG, ".csv")), row.names = FALSE)
  log_sanity("progeny", "pathways scored", nrow(acts), "OK")
  d <- acts %>% dplyr::mutate(pathway = factor(pathway, levels = pathway[order(score)]),
                              sig = ifelse(FDR < FDR_CUT, "FDR < 0.05", "n.s."))
  p <- ggplot(d, aes(score, pathway, fill = score > 0)) +
    geom_col(aes(alpha = sig), color = "black", linewidth = 0.25) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.4), name = NULL) +
    scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                      labels = c(`TRUE` = "active", `FALSE` = "suppressed"), name = NULL) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = "PROGENy pathway activity", y = NULL,
         title = "Upstream signalling activity in KAT8-KD cells",
         caption = "PROGENy uses perturbation-derived responsive genes -- a different gene set from the TF regulons.") +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("E_PROGENy_cells_", RUN_TAG, ".png")), p,
         width = 7.5, height = 5.5, dpi = 300, bg = "white")
  acts
}
progeny_res <- tryCatch(run_progeny(), error = function(e) { message("PROGENy failed: ", conditionMessage(e)); NULL })

## =========================================================
## MODULE F: KAT8 COMPLEX OUTPUT (MSL / NSL) + specificity controls
## =========================================================
cat("\n========== F. KAT8 COMPLEX (MSL/NSL) OUTPUT ==========\n")
## KAT8 is the catalytic H4K16ac subunit of two complexes:
##   MSL (Msl1/2/3)                    -> broad H4K16ac / transcriptional output
##   NSL/KANSL (Kansl1-3, Mcrs1, Phf20, Ogt, Wdr5, Hcfc1)
##                                     -> promoter-proximal; documented direct
##                                        regulator of nuclear-encoded OXPHOS
## siRNA removes the catalytic core but NOT the partners' mRNA, so the complex
## effect is invisible in membership and must be read from regulatory OUTPUT.
## Two specificity controls decide whether an OXPHOS shift is really NSL:
##   mito-encoded genes (mt-*) -- transcribed by mtDNA, OUTSIDE NSL control
##   ribosomal proteins        -- comparably abundant housekeeping genes
## If nuclear OXPHOS moves while BOTH controls stay flat, the shift is specific.
msl <- c("Msl1","Msl2","Msl3")
nsl <- c("Kansl1","Kansl2","Kansl3","Mcrs1","Phf20","Ogt","Wdr5","Hcfc1")
complex_members <- tibble::tibble(gene = c("Kat8", msl, nsl),
                                  complex = c("catalytic", rep("MSL", length(msl)), rep("NSL/KANSL", length(nsl))))
complex_expr <- complex_members %>%
  dplyr::left_join(cells %>% dplyr::select(gene, log2FC, padj, stat), by = "gene")
write.csv(complex_expr, file.path(OUTDIR, paste0("F_KAT8_complex_membership_", RUN_TAG, ".csv")), row.names = FALSE)
log_sanity("complex", "complex members found in cells",
           paste0(sum(!is.na(complex_expr$log2FC)), "/", nrow(complex_expr)), "OK")

nuclear_oxphos <- c("Ndufa1","Ndufa2","Ndufa3","Ndufa4","Ndufa5","Ndufa8","Ndufa9","Ndufa11","Ndufa13",
  "Ndufb3","Ndufb7","Ndufb8","Ndufb9","Ndufb10","Ndufb11","Ndufs2","Ndufs3","Ndufs6","Ndufs7","Ndufv1",
  "Sdhb","Sdhc","Sdhd","Uqcrb","Uqcrq","Uqcrh","Uqcr10","Uqcr11","Uqcrc1","Uqcrfs1","Cyc1","Cycs",
  "Cox4i1","Cox5b","Cox6a1","Cox6b1","Cox6c","Cox7a2","Cox7b","Cox7c","Cox10","Atp5f1a","Atp5f1b",
  "Atp5f1d","Atp5mc1","Atp5mc2","Atp5mc3","Atp5me","Atp5mf","Atp5pd")

complex_rows <- list()
for (nm in c(CELL_KEY, tissues)) {
  d <- de[[nm]]
  mito_g <- grep("^mt-", d$gene, value = TRUE)
  ribo_g <- grep("^Rp[sl][0-9]", d$gene, value = TRUE)
  prog <- function(gs, label) {
    p <- intersect(gs, d$gene); if (length(p) < 4) return(NULL)
    sub <- d[d$gene %in% p, ]; bg <- d[!d$gene %in% p, ]
    pv <- tryCatch(suppressWarnings(wilcox.test(sub$stat, bg$stat)$p.value), error = function(e) NA_real_)
    data.frame(dataset = nm, program = label, n = length(p),
               mean_log2FC = mean(sub$log2FC, na.rm = TRUE),
               mean_Wald = mean(sub$stat, na.rm = TRUE),
               p_vs_background = pv, stringsAsFactors = FALSE)
  }
  complex_rows <- c(complex_rows, list(
    prog(nuclear_oxphos, "Nuclear-encoded OXPHOS (NSL targets)"),
    prog(mito_g,         "Mito-encoded (control: outside NSL)"),
    prog(ribo_g,         "Ribosomal proteins (control: housekeeping)")))
}
complex_out <- dplyr::bind_rows(complex_rows)
write.csv(complex_out, file.path(OUTDIR, paste0("F_KAT8_complex_output_", RUN_TAG, ".csv")), row.names = FALSE)
cat("\n"); print(complex_out, row.names = FALSE)

## Data-driven verdict (never hardcode the conclusion in a caption)
cc <- complex_out %>% dplyr::filter(dataset == CELL_KEY)
nuc <- cc %>% dplyr::filter(grepl("^Nuclear", program))
ctl <- cc %>% dplyr::filter(grepl("control", program))
verdict <- "inconclusive"
if (nrow(nuc) == 1 && nrow(ctl) >= 1) {
  nuc_sig  <- isTRUE(nuc$p_vs_background < FDR_CUT)
  ctl_flat <- all(!is.finite(ctl$p_vs_background) | ctl$p_vs_background > FDR_CUT)
  verdict <- if (nuc_sig && ctl_flat)
    sprintf("SPECIFIC: nuclear OXPHOS shifts (mean Wald %+.2f, p=%.1e) while controls stay flat",
            nuc$mean_Wald, nuc$p_vs_background)
  else if (nuc_sig) "nuclear OXPHOS shifts BUT a control also shifts -> not NSL-specific"
  else "no nuclear-OXPHOS shift detected in cells"
}
log_sanity("complex", "NSL specificity verdict (cells)", verdict, "INFO")

plot_df <- complex_out %>% dplyr::filter(is.finite(mean_Wald))
p <- ggplot(plot_df, aes(x = mean_Wald, y = program, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
  labs(x = "mean DESeq2 Wald statistic vs background", y = NULL,
       title = "KAT8 complex regulatory output, with specificity controls",
       caption = paste0("NSL targets = nuclear-encoded OXPHOS. Controls: mito-encoded (mtDNA-transcribed, outside NSL) ",
                        "and ribosomal proteins.\nVerdict (cells): ", verdict)) +
  theme_pub()
ggsave(file.path(OUTDIR, paste0("F_KAT8_complex_output_", RUN_TAG, ".png")), p,
       width = 10, height = 4.6, dpi = 300, bg = "white")

## =========================================================
## MODULE G: CONVERGENCE -- which TFs have multi-layer support?
## =========================================================
cat("\n=============== G. CROSS-LAYER CONVERGENCE ===============\n")
conv <- data.frame(TF = character(), layer1_mRNA = character(),
                   layer2_activity = character(), layer3_TFT = character(),
                   n_layers = integer(), stringsAsFactors = FALSE)
cand <- unique(c(tf_deg$gene,
                 if (!is.null(tf_act)) tf_act$TF[tf_act$FDR < FDR_CUT] else character(0)))
if (length(cand) > 0) {
  conv <- data.frame(TF = cand, stringsAsFactors = FALSE) %>%
    dplyr::left_join(tf_deg %>% dplyr::transmute(TF = gene, l1 = sprintf("%+.2f", log2FC)), by = "TF")
  if (!is.null(tf_act)) {
    conv <- conv %>% dplyr::left_join(
      tf_act %>% dplyr::filter(FDR < FDR_CUT) %>% dplyr::transmute(TF, l2 = sprintf("%+.2f", activity)), by = "TF")
  } else conv$l2 <- NA_character_
  if (!is.null(tft_res)) {
    hits <- tft_res %>% dplyr::filter(padj < FDR_CUT)
    conv$l3 <- vapply(conv$TF, function(t) {
      m <- hits$pathway[grepl(paste0("(^|_)", toupper(t), "(_|$)"), toupper(hits$pathway))]
      if (length(m) == 0) NA_character_ else sprintf("%+.2f", hits$NES[hits$pathway == m[1]][1])
    }, character(1))
  } else conv$l3 <- NA_character_
  conv$n_layers <- rowSums(!is.na(conv[, c("l1","l2","l3")]))
  conv <- conv %>% dplyr::arrange(dplyr::desc(n_layers), TF) %>%
    dplyr::rename(layer1_mRNA_log2FC = l1, layer2_activity = l2, layer3_TFT_NES = l3)
  write.csv(conv, file.path(OUTDIR, paste0("G_TF_convergence_", RUN_TAG, ".csv")), row.names = FALSE)
  cat("\nTFs supported by >=2 independent layers:\n")
  multi <- conv %>% dplyr::filter(n_layers >= 2)
  if (nrow(multi) > 0) print(multi, row.names = FALSE) else cat("  (none)\n")
  log_sanity("convergence", "TFs with >=2 layers of support", nrow(multi), "INFO")
}

## ---------------------------------------------------------
## Write the sanity ledger + session info
## ---------------------------------------------------------
sanity_df <- dplyr::bind_rows(SANITY)
write.csv(sanity_df, file.path(OUTDIR, paste0("SANITY_part5_", RUN_TAG, ".csv")), row.names = FALSE)
n_fail <- sum(sanity_df$status == "FAIL"); n_warn <- sum(sanity_df$status == "WARN")
cat("\n==========================================================\n")
cat("PART 5 SANITY SUMMARY: ", nrow(sanity_df), " checks | ",
    n_fail, " FAIL | ", n_warn, " WARN\n", sep = "")
if (n_fail > 0) { cat("FAILED CHECKS:\n"); print(sanity_df[sanity_df$status == "FAIL", ], row.names = FALSE) }
cat("==========================================================\n")

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, paste0("sessionInfo_part5_", RUN_TAG, ".txt")))
cat("\n[DONE] Part 5 complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Read in this order: A (what the cell broadcasts) -> G (which TF drives it)\n")
