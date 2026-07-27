#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   TF analysis: 3 layers + convergence
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## PART 5b - TRANSCRIPTION FACTOR ANALYSIS
## Question: which TF drives the KAT8-null program?
## =========================================================
## KAT8 is an H4K16ac co-activator, not a TF. So the mechanistic question is
## which TF programs shift when it is lost -- especially in the pure 3T3-L1
## adipocytes, where there is no immune infiltrate and no systemic context.
##
## THREE LAYERS, DELIBERATELY INDEPENDENT
##   Layer 1  TF genes that themselves change      -- no regulon, no inference
##   Layer 2  regulon activity (decoupleR ULM)     -- curated TF->target network
##   Layer 3  MSigDB C3:TFT target enrichment      -- ChIP-derived, independent
##                                                    of the layer-2 regulon
##
## WHY NOT JUST LAYER 2: regulon activity is a re-weighted score computed from
## the SAME statistics as the DE result, so on its own it partly restates what
## you already know, and it depends on a curated, largely non-adipose regulon
## that is thin for adipocyte TFs. It is also blind to co-activators -- the very
## class KAT8 belongs to. Layer 1 assumes nothing; layer 3 uses a different
## evidence base entirely. AGREEMENT ACROSS LAYERS is the evidence, not any one
## layer. Section 6 tabulates that agreement.
##
## Also included: PROGENy (signalling one level above TFs) and a cross-contrast
## comparison so a TF active in cells AND tissue can be separated from one that
## is tissue-only (i.e. likely driven by infiltrating cells).
##
## Input : DE CSVs configured in config_downstream.R      Output: OUTDIR
## =========================================================

required_cran <- c("dplyr","tidyr","ggplot2","tibble","pheatmap","RColorBrewer")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (bp in c("org.Mm.eg.db","AnnotationDbi","GO.db","decoupleR","fgsea"))
  if (!requireNamespace(bp, quietly = TRUE))
    tryCatch(BiocManager::install(bp, ask = FALSE, update = FALSE), error = function(e) NULL)
if (!requireNamespace("msigdbr", quietly = TRUE))
  tryCatch(install.packages("msigdbr"), error = function(e) NULL)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(tibble)
  library(pheatmap); library(RColorBrewer)
  library(org.Mm.eg.db); library(AnnotationDbi); library(GO.db)
})

source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
init_sanity("part5b_tf")
`%||%` <- function(a, b) if (is.null(a)) b else a
tp <- tf_params

de <- load_all_de()
check_gene_names(de)
cells   <- de[[CELL_KEY]]
tissues <- intersect(TISSUE_KEYS, names(de))
cell_sig <- cells %>% dplyr::filter(!is.na(padj), padj < FDR_CUT)

## =========================================================
## LAYER 1 - TF genes that themselves change (zero assumptions)
## =========================================================
cat("\n---- Layer 1: TF genes ----\n")
tf_universe <- go_genes(GO_IDS[["tf_activity"]])   ## GO:0003700
log_sanity("GO TF universe size", length(tf_universe),
           ifelse(length(tf_universe) > 500, "OK", "WARN"))

tf_deg <- cell_sig %>% dplyr::filter(gene %in% tf_universe) %>%
  dplyr::arrange(dplyr::desc(log2FC)) %>%
  dplyr::select(gene, log2FC, stat, padj)
log_sanity("TF genes changed in cells (FDR<0.05)", nrow(tf_deg), "OK")
write.csv(tf_deg, file.path(OUTDIR, paste0("part5b_L1_TF_genes_", RUN_TAG, ".csv")), row.names = FALSE)

if (nrow(tf_deg) >= 2) {
  d <- tf_deg %>% dplyr::mutate(gene = factor(gene, levels = gene[order(log2FC)]))
  p <- ggplot(d, aes(log2FC, gene, fill = log2FC > 0)) +
    geom_col(color = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                      labels = c(`TRUE` = "Up in KAT8 KD", `FALSE` = "Down in KAT8 KD"), name = NULL) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = expression(log[2]~fold~change), y = NULL,
         title = "Layer 1: TF genes changed in KAT8-KD cells",
         caption = paste0("TFs = GO:0003700 (DNA-binding transcription factor activity), FDR<", FDR_CUT,
                          ". No regulon or inference involved -- these are measured directly.")) +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("part5b_L1_TF_genes_", RUN_TAG, ".png")),
         p, width = 7.5, height = max(3.5, nrow(tf_deg) * 0.28), dpi = 300, bg = "white")
}

## =========================================================
## LAYER 2 - regulon activity (decoupleR), all contrasts
## =========================================================
cat("\n---- Layer 2: regulon activity ----\n")
get_regulon <- function() {
  if (!is.null(tp$network_csv) && file.exists(tp$network_csv)) {
    log_sanity("regulon source", basename(tp$network_csv), "OK")
    return(utils::read.csv(tp$network_csv, stringsAsFactors = FALSE))
  }
  if (file.exists(tp$network_cache)) {
    log_sanity("regulon source", "cached RDS", "OK"); return(readRDS(tp$network_cache))
  }
  n <- tryCatch(decoupleR::get_collectri(organism = ORGANISM, split_complexes = FALSE),
                error = function(e) NULL)
  if (is.null(n)) n <- tryCatch(decoupleR::get_dorothea(organism = ORGANISM, levels = c("A","B","C")),
                                error = function(e) NULL)
  log_sanity("regulon source", ifelse(is.null(n), "UNAVAILABLE", "OmniPath"),
             ifelse(is.null(n), "WARN", "OK"))
  if (!is.null(n)) saveRDS(n, tp$network_cache)
  n
}
net <- get_regulon()
tf_act <- NULL

if (!is.null(net) && all(c("source","target","mor") %in% colnames(net))) {
  log_sanity("regulon TFs / edges", paste0(length(unique(net$source)), " / ", nrow(net)), "OK")

  score_one <- function(tbl, label) {
    v <- tbl %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE)
    m <- matrix(v$stat, ncol = 1, dimnames = list(v$gene, label))
    a <- tryCatch(decoupleR::run_ulm(m, net, .source = "source", .target = "target",
                                     .mor = "mor", minsize = tp$min_regulon_size),
                  error = function(e) NULL)
    if (is.null(a)) return(NULL)
    a %>% dplyr::transmute(contrast = label, TF = source, activity = score, p_value) %>%
      dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))
  }
  all_act <- dplyr::bind_rows(lapply(names(de), function(n) score_one(de[[n]], n)))
  write.csv(all_act, file.path(OUTDIR, paste0("part5b_L2_TF_activity_all_", RUN_TAG, ".csv")),
            row.names = FALSE)

  tf_act <- all_act %>% dplyr::filter(contrast == CELL_KEY) %>%
    dplyr::arrange(dplyr::desc(abs(activity)))
  log_sanity("TFs scored in cells", nrow(tf_act), "OK")
  log_sanity("TFs significant in cells (FDR<0.05)", sum(tf_act$FDR < FDR_CUT, na.rm = TRUE), "OK")

  ## Cross-annotate with the TF's OWN mRNA. For an activator's writer this is
  ## the direct/indirect discriminator: activity down AND the TF's own message
  ## down is consistent with KAT8 being needed to express that TF; activity
  ## changed while its mRNA is flat implies post-translational / indirect.
  own <- cells %>% dplyr::select(gene, own_log2FC = log2FC, own_padj = padj)
  tf_act <- tf_act %>% dplyr::left_join(own, by = c("TF" = "gene")) %>%
    dplyr::mutate(
      own_mrna = dplyr::case_when(
        is.na(own_padj)                  ~ "not tested",
        own_padj < FDR_CUT & own_log2FC < 0 ~ "TF mRNA down",
        own_padj < FDR_CUT & own_log2FC > 0 ~ "TF mRNA up",
        TRUE                             ~ "TF mRNA unchanged"),
      interpretation = dplyr::case_when(
        activity < 0 & own_mrna == "TF mRNA down" ~ "direct-target candidate (TF expression lost)",
        activity > 0 & own_mrna == "TF mRNA up"   ~ "TF induced and active",
        own_mrna == "TF mRNA unchanged"           ~ "post-translational / indirect",
        TRUE ~ "see own_mrna"))
  write.csv(tf_act, file.path(OUTDIR, paste0("part5b_L2_TF_activity_cells_", RUN_TAG, ".csv")),
            row.names = FALSE)

  ## Leading edge -- which targets drove each call, so every TF is auditable.
  ## NOTE: the stat vector MUST be named. An earlier version of this analysis
  ## used an unnamed vector, so the name-based lookup matched nothing and the
  ## leading-edge table came out silently EMPTY.
  v <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE)
  stat_named <- setNames(v$stat, v$gene)
  le <- net %>% dplyr::filter(source %in% head(tf_act$TF, 20)) %>%
    dplyr::mutate(target_stat = unname(stat_named[target]), contribution = mor * target_stat) %>%
    dplyr::filter(!is.na(contribution)) %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(dplyr::desc(abs(contribution)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 15) %>% dplyr::ungroup()
  write.csv(le, file.path(OUTDIR, paste0("part5b_L2_leading_edge_", RUN_TAG, ".csv")), row.names = FALSE)
  log_sanity("leading-edge rows (0 = the old empty-table bug)", nrow(le),
             ifelse(nrow(le) > 0, "OK", "FAIL"))

  ## Cells barplot
  top <- tf_act %>% dplyr::slice_head(n = tp$top_n_cells) %>%
    dplyr::mutate(TF = factor(TF, levels = TF[order(activity)]),
                  sig = ifelse(FDR < FDR_CUT, "FDR < 0.05", "n.s."))
  p <- ggplot(top, aes(activity, TF, fill = activity)) +
    geom_col(aes(alpha = sig), color = "black", linewidth = 0.25) +
    scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.4), name = NULL) +
    scale_fill_gradient2(low = col_down, mid = "grey95", high = col_up, midpoint = 0, name = "activity") +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(x = "Inferred TF activity (decoupleR ULM)", y = NULL,
         title = "Layer 2: TF regulon activity in KAT8-KD cells",
         caption = "A re-weighted score on the same statistics -- treat as a summary; the evidence is agreement with layers 1 and 3.") +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("part5b_L2_TF_activity_cells_", RUN_TAG, ".png")),
         p, width = 8, height = 7.5, dpi = 300, bg = "white")

  ## Cross-contrast heatmap: cell-autonomous vs tissue-only regulators.
  if (length(tissues) > 0) {
    wide <- all_act %>% dplyr::select(TF, contrast, activity) %>%
      tidyr::pivot_wider(names_from = contrast, values_from = activity)
    strong <- all_act %>% dplyr::group_by(TF) %>%
      dplyr::summarise(m = max(abs(activity), na.rm = TRUE), .groups = "drop") %>%
      dplyr::slice_max(m, n = tp$top_n_heatmap)
    hm <- wide %>% dplyr::filter(TF %in% strong$TF) %>%
      tibble::column_to_rownames("TF") %>% as.matrix()
    hm <- hm[, intersect(c(CELL_KEY, tissues), colnames(hm)), drop = FALSE]
    hm[is.na(hm)] <- 0
    lim <- max(abs(hm), na.rm = TRUE)
    png(file.path(OUTDIR, paste0("part5b_L2_TF_heatmap_", RUN_TAG, ".png")),
        width = 1400, height = 2100, res = 220)
    pheatmap(hm, color = colorRampPalette(c(col_down, "#4393C3", "grey97", "#D6604D", col_up))(100),
             breaks = seq(-lim, lim, length.out = 101), cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.1f", fontsize_number = 7, angle_col = "0",   ## must be a STRING: pheatmap match.arg()s it
             main = "TF activity by contrast (cells column = cell-autonomous)")
    dev.off()
    log_sanity("cross-contrast heatmap TFs", nrow(hm), "OK")
  }

  ## Optional permutation calibration at n = 4/group.
  ##
  ## PERFORMANCE: this used to call run_ulm once PER SHUFFLE (1000 separate
  ## calls), each re-intersecting the network with ~13k genes -- 10-30 minutes.
  ## decoupleR vectorises over the COLUMNS of the input matrix, so the whole
  ## null is computed in a handful of calls instead. Shuffling gene labels is
  ## equivalent to permuting the statistic vector with the gene names fixed,
  ## which is what each column below is.
  if (isTRUE(tp$n_permutations > 0)) {
    nperm <- tp$n_permutations
    chunk <- max(1, min(tp$perm_chunk %||% 100, nperm))
    cat("  permutation calibration: ", nperm, " shuffles in chunks of ", chunk, "\n", sep = "")
    t0 <- Sys.time(); nullmax <- numeric(0); done <- 0
    while (done < nperm) {
      k <- min(chunk, nperm - done)
      pm <- vapply(seq_len(k), function(i) sample(v$stat), numeric(nrow(v)))
      dimnames(pm) <- list(v$gene, paste0("perm", done + seq_len(k)))
      a <- decoupleR::run_ulm(pm, net, .source = "source", .target = "target",
                              .mor = "mor", minsize = tp$min_regulon_size)
      mx <- a %>% dplyr::group_by(condition) %>%
        dplyr::summarise(m = max(abs(score), na.rm = TRUE), .groups = "drop")
      nullmax <- c(nullmax, mx$m); done <- done + k
      cat("    ", done, "/", nperm, "  (", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s)\n", sep = "")
    }
    ## Empirical p against the null of the MAXIMUM |score| -> family-wise
    ## control across TFs, deliberately stricter than the per-TF BH FDR.
    tf_act$emp_p <- vapply(abs(tf_act$activity),
                           function(x) (1 + sum(nullmax >= x)) / (1 + length(nullmax)), numeric(1))
    write.csv(tf_act, file.path(OUTDIR, paste0("part5b_L2_TF_activity_cells_", RUN_TAG, ".csv")),
              row.names = FALSE)
    log_sanity("permutation time (s)",
               round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "INFO")
    log_sanity("TFs passing empirical p<0.05 (family-wise)",
               sum(tf_act$emp_p < 0.05, na.rm = TRUE), "INFO")
  }
} else {
  log_sanity("Layer 2", "no usable regulon - layer skipped", "WARN")
}

## =========================================================
## LAYER 3 - independent ChIP-derived TF target sets
## =========================================================
cat("\n---- Layer 3: MSigDB C3:TFT ----\n")
tft <- NULL
if (requireNamespace("msigdbr", quietly = TRUE) && requireNamespace("fgsea", quietly = TRUE)) {
  ## msigdbr >= 10 renamed category/subcategory -> collection/subcollection.
  try_msig <- function(...) tryCatch(msigdbr::msigdbr(...), error = function(e) NULL)
  msig <- try_msig(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD")
  if (is.null(msig) || nrow(msig) == 0) msig <- try_msig(species = "Mus musculus", category = "C3", subcategory = "TFT:GTRD")
  if (is.null(msig) || nrow(msig) == 0) msig <- try_msig(species = "Mus musculus", collection = "C3")
  if (is.null(msig) || nrow(msig) == 0) msig <- try_msig(species = "Mus musculus", category = "C3")
  if (!is.null(msig) && nrow(msig) > 0 && any(grepl("^TFT|GTRD", msig$gs_name)))
    msig <- msig[grepl("^TFT|GTRD", msig$gs_name), , drop = FALSE]

  if (!is.null(msig) && nrow(msig) > 0) {
    sc <- if ("gene_symbol" %in% colnames(msig)) "gene_symbol" else "human_gene_symbol"
    sets <- split(msig[[sc]], msig$gs_name)
    sets <- sets[vapply(sets, function(x) length(unique(x)) >= 10, logical(1))]
    log_sanity("TFT gene sets tested", length(sets), "OK")
    vv <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE) %>%
      dplyr::arrange(dplyr::desc(stat))
    tft <- tryCatch(fgsea::fgsea(sets, setNames(vv$stat, vv$gene), minSize = 10, maxSize = 2000, eps = 0),
                    error = function(e) NULL)
    if (!is.null(tft)) {
      tft <- as.data.frame(tft) %>% dplyr::arrange(padj) %>%
        dplyr::mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ","))
      write.csv(tft, file.path(OUTDIR, paste0("part5b_L3_TFT_enrichment_", RUN_TAG, ".csv")), row.names = FALSE)
      log_sanity("TFT sets significant (FDR<0.05)", sum(tft$padj < FDR_CUT, na.rm = TRUE), "OK")
      top <- tft %>% dplyr::filter(padj < FDR_CUT) %>%
        dplyr::arrange(dplyr::desc(abs(NES))) %>% dplyr::slice_head(n = 25)
      if (nrow(top) >= 2) {
        top <- top %>% dplyr::mutate(pathway = factor(pathway, levels = pathway[order(NES)]))
        p <- ggplot(top, aes(NES, pathway, fill = NES > 0)) +
          geom_col(color = "black", linewidth = 0.25) +
          scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                            labels = c(`TRUE` = "targets up", `FALSE` = "targets down"), name = NULL) +
          geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
          labs(x = "NES", y = NULL,
               title = "Layer 3: TF target-set enrichment (ChIP-derived)",
               caption = "MSigDB C3:TFT (GTRD). Independent of the layer-2 regulon -- convergence between them is the evidence.") +
          theme_pub() + theme(axis.text.y = element_text(size = 6.5))
        ggsave(file.path(OUTDIR, paste0("part5b_L3_TFT_enrichment_", RUN_TAG, ".png")),
               p, width = 9.5, height = max(5, nrow(top) * 0.28), dpi = 300, bg = "white")
      }
    }
  } else log_sanity("MSigDB C3 retrieved", "NO - layer skipped", "WARN")
} else log_sanity("msigdbr/fgsea available", "NO - layer skipped", "WARN")

## =========================================================
## PROGENy - signalling one level above the TFs
## =========================================================
cat("\n---- PROGENy ----\n")
if (requireNamespace("decoupleR", quietly = TRUE)) {
  pnet <- tryCatch(decoupleR::get_progeny(organism = ORGANISM, top = 500), error = function(e) NULL)
  if ((is.null(pnet) || nrow(pnet) < 50) && requireNamespace("progeny", quietly = TRUE)) {
    pnet <- tryCatch({
      m <- progeny::getModel(organism = "Mouse", top = 500); m$gene <- rownames(m)
      tidyr::pivot_longer(m, cols = -gene, names_to = "source", values_to = "weight") %>%
        dplyr::filter(weight != 0) %>% dplyr::transmute(source, target = gene, weight)
    }, error = function(e) NULL)
  }
  if (!is.null(pnet) && nrow(pnet) >= 50) {
    vv <- cells %>% dplyr::filter(is.finite(stat)) %>% dplyr::distinct(gene, .keep_all = TRUE)
    pa <- tryCatch(decoupleR::run_mlm(matrix(vv$stat, ncol = 1, dimnames = list(vv$gene, "cells")),
                                      pnet, .source = "source", .target = "target",
                                      .mor = "weight", minsize = 5), error = function(e) NULL)
    if (!is.null(pa)) {
      pa <- pa %>% dplyr::transmute(pathway = source, score, p_value) %>%
        dplyr::mutate(FDR = p.adjust(p_value, method = "BH")) %>% dplyr::arrange(dplyr::desc(score))
      write.csv(pa, file.path(OUTDIR, paste0("part5b_PROGENy_", RUN_TAG, ".csv")), row.names = FALSE)
      log_sanity("PROGENy pathways scored", nrow(pa), "OK")
      d <- pa %>% dplyr::mutate(pathway = factor(pathway, levels = pathway[order(score)]),
                                sig = ifelse(FDR < FDR_CUT, "FDR < 0.05", "n.s."))
      p <- ggplot(d, aes(score, pathway, fill = score > 0)) +
        geom_col(aes(alpha = sig), color = "black", linewidth = 0.25) +
        scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.4), name = NULL) +
        scale_fill_manual(values = c(`TRUE` = col_up, `FALSE` = col_down),
                          labels = c(`TRUE` = "active", `FALSE` = "suppressed"), name = NULL) +
        geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
        labs(x = "PROGENy pathway activity", y = NULL,
             title = "Upstream signalling activity in KAT8-KD cells",
             caption = "Perturbation-derived responsive genes -- a different gene set from the TF regulons.") +
        theme_pub()
      ggsave(file.path(OUTDIR, paste0("part5b_PROGENy_", RUN_TAG, ".png")), p,
             width = 7.5, height = 5.5, dpi = 300, bg = "white")
    }
  } else log_sanity("PROGENy model", "unavailable - skipped", "WARN")
}

## =========================================================
## 6) CONVERGENCE - which TFs have MULTI-LAYER support?
## =========================================================
cat("\n---- Convergence across layers ----\n")
cand <- unique(c(tf_deg$gene, if (!is.null(tf_act)) tf_act$TF[tf_act$FDR < FDR_CUT] else character(0)))
if (length(cand) > 0) {
  conv <- data.frame(TF = cand, stringsAsFactors = FALSE) %>%
    dplyr::left_join(tf_deg %>% dplyr::transmute(TF = gene, layer1_mRNA_log2FC = sprintf("%+.2f", log2FC)), by = "TF")
  ## Layer 2 counts only if the TF was SIGNIFICANT there.
  conv$layer2_activity <- NA_character_
  if (!is.null(tf_act)) {
    sig2 <- tf_act[!is.na(tf_act$FDR) & tf_act$FDR < FDR_CUT, c("TF", "activity")]
    i <- match(conv$TF, sig2$TF)
    conv$layer2_activity[!is.na(i)] <- sprintf("%+.2f", sig2$activity[i[!is.na(i)]])
  }
  conv$layer3_TFT_NES <- if (!is.null(tft)) {
    hits <- tft %>% dplyr::filter(padj < FDR_CUT)
    vapply(conv$TF, function(t) {
      m <- hits$pathway[grepl(paste0("(^|_)", toupper(t), "(_|$)"), toupper(hits$pathway))]
      if (length(m) == 0) NA_character_ else sprintf("%+.2f", hits$NES[hits$pathway == m[1]][1])
    }, character(1))
  } else NA_character_
  conv$n_layers <- rowSums(!is.na(conv[, c("layer1_mRNA_log2FC","layer2_activity","layer3_TFT_NES")]))
  conv <- conv %>% dplyr::arrange(dplyr::desc(n_layers), TF)
  write.csv(conv, file.path(OUTDIR, paste0("part5b_convergence_", RUN_TAG, ".csv")), row.names = FALSE)
  multi <- conv %>% dplyr::filter(n_layers >= 2)
  cat("\nTFs supported by >=2 independent layers:\n")
  if (nrow(multi) > 0) print(multi, row.names = FALSE) else cat("  (none)\n")
  log_sanity("TFs with >=2 layers of support", nrow(multi), "INFO")
}

write_sanity()
cat("\n[DONE] Part 5b complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Interpret the CONVERGENCE table first -- a TF supported by 2-3 layers is a\n")
cat("real candidate; a hit in layer 2 alone is a summary statistic, not evidence.\n")
