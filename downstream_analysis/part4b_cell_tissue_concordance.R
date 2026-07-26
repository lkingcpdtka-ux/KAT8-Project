#!/usr/bin/env Rscript

## =========================================================
## KAT8 - PART 4b: CELL <-> TISSUE CONCORDANCE
##                (direct/indirect & cell-autonomy)
## =========================================================
## The Adipoq-Cre KO tissue is a cell MIXTURE, so a tissue change can
## be adipocyte-intrinsic OR a shift in cell composition (infiltrating
## macrophages, fibrosis). The 3T3-L1 cells are a pure adipocyte system
## with no infiltrate, so overlaying the two separates:
##
##   concordant (cell & tissue, same sign)  -> CELL-AUTONOMOUS
##                                             (strongest direct-target candidates)
##   tissue-only, esp. inflammatory up      -> NON-CELL-AUTONOMOUS (indirect)
##   cell-only                              -> in-vitro-specific
##
## Because KAT8 is an ACTIVATOR, its direct consequence should be
## DOWN-regulation; so the concordant-down program is the prime
## direct/cell-autonomous candidate set.
##
## Everything here is THRESHOLD-FREE on the cell side (rank correlation,
## gene-set enrichment, directional tests that borrow power from the
## well-powered tissue DEGs), so a low cell DEG count does not limit it.
## =========================================================

## 1) Packages ----------------------------------------------
required_cran <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "scales")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea", ask = FALSE, update = FALSE)
have_qvalue <- requireNamespace("qvalue", quietly = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(scales)
  library(fgsea)
})

## 2) Config -------------------------------------------------
source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
cp <- concordance_params
de <- load_all_de()

summary_rows <- list()

## 3) Per tissue depot: compare against cells ---------------
analyse_pair <- function(cell_tbl, tis_tbl, tis_label) {
  cat("\n==============================================\n")
  cat("  CELLS  vs  ", tis_label, "\n", sep = "")
  cat("==============================================\n")

  ## ---- shared expressed-gene universe (merge on symbol) ----
  m <- dplyr::inner_join(
    cell_tbl %>% dplyr::select(gene, stat_cell = stat, lfc_cell = log2FC,
                               p_cell = pvalue, padj_cell = padj),
    tis_tbl  %>% dplyr::select(gene, stat_tis = stat, lfc_tis = log2FC,
                               p_tis = pvalue, padj_tis = padj),
    by = "gene"
  ) %>% dplyr::filter(!is.na(stat_cell), !is.na(stat_tis))
  cat("[INFO] shared genes:", nrow(m), "\n")

  ## ---- (a) genome-wide rank concordance of Wald stats ----
  rho_all <- suppressWarnings(cor(m$stat_cell, m$stat_tis, method = cp$cor_method))
  ## restricted to the well-powered tissue DEGs
  tis_sig <- m %>% dplyr::filter(!is.na(padj_tis), padj_tis < cp$tissue_fdr,
                                 abs(lfc_tis) >= cp$tissue_lfc)
  rho_sig <- if (nrow(tis_sig) > 5)
    suppressWarnings(cor(tis_sig$stat_cell, tis_sig$stat_tis, method = cp$cor_method)) else NA_real_
  cat(sprintf("[RESULT] Spearman(stat) all genes = %.3f | within tissue DEGs = %.3f (n=%d)\n",
              rho_all, rho_sig, nrow(tis_sig)))

  ## ---- (b) directional concordance, borrowing power from tissue ----
  ## Among tissue DEGs, how often does the cell effect share the sign?
  dir_test <- function(df, dirlabel) {
    if (nrow(df) < 3) return(NULL)
    same <- sum(sign(df$stat_cell) == sign(df$stat_tis))
    bt <- binom.test(same, nrow(df), p = 0.5, alternative = "greater")
    cat(sprintf("[RESULT] tissue %-12s DEGs: %d/%d same sign in cells (%.0f%%), binom p = %.2e\n",
                dirlabel, same, nrow(df), 100 * same / nrow(df), bt$p.value))
    data.frame(comparison = paste0("cells_vs_", tis_label), set = dirlabel,
               n = nrow(df), same_sign = same, frac = same / nrow(df), binom_p = bt$p.value)
  }
  r1 <- dir_test(tis_sig %>% dplyr::filter(lfc_tis < 0), "DOWN(direct?)")
  r2 <- dir_test(tis_sig %>% dplyr::filter(lfc_tis > 0), "UP(indirect?)")

  ## ---- (c) pi1: fraction of tissue DEGs truly non-null in cells ----
  pi1 <- NA_real_
  if (have_qvalue && nrow(tis_sig) > 20) {
    pv <- tis_sig$p_cell[!is.na(tis_sig$p_cell)]
    pi1 <- tryCatch(1 - qvalue::qvalue(pv)$pi0, error = function(e) NA_real_)
    cat(sprintf("[RESULT] pi1 (tissue DEGs non-null in cells) = %.2f\n", pi1))
  }

  ## ---- (d) signature GSEA: is the tissue PROGRAM shifted in cells? ----
  ## Rank ALL cell genes by Wald stat; test tissue up/down DEG sets.
  ranks <- m %>% dplyr::arrange(dplyr::desc(stat_cell))
  rank_vec <- setNames(ranks$stat_cell, ranks$gene)
  sets <- list(
    tissue_DOWN = tis_sig$gene[tis_sig$lfc_tis < 0],
    tissue_UP   = tis_sig$gene[tis_sig$lfc_tis > 0]
  )
  sets <- sets[vapply(sets, length, integer(1)) >= 5]
  gsea <- NULL
  if (length(sets) > 0) {
    gsea <- fgsea::fgsea(pathways = sets, stats = rank_vec,
                         minSize = 5, maxSize = 2000, eps = 0) %>%
      as.data.frame() %>%
      dplyr::mutate(comparison = paste0("cells_vs_", tis_label)) %>%
      dplyr::select(comparison, pathway, NES, pval, padj, size)
    print(gsea[, c("pathway", "NES", "pval", "padj", "size")])
    write.csv(gsea, file.path(OUTDIR, paste0("concordance_signatureGSEA_", tis_label, "_", RUN_TAG, ".csv")),
              row.names = FALSE)
  }

  ## ---- (e) labelled log2FC quadrant plot ----
  lab_genes <- c(cp$marker_metabolic, cp$marker_immune, cp$marker_anchor)
  m <- m %>% dplyr::mutate(
    quadrant = dplyr::case_when(
      lfc_cell < 0 & lfc_tis < 0 ~ "concordant down (cell-autonomous)",
      lfc_cell > 0 & lfc_tis > 0 ~ "concordant up",
      abs(lfc_cell) < 0.2 & abs(lfc_tis) > 0.5 ~ "tissue-only (non-cell-autonomous)",
      TRUE ~ "other / discordant"
    ),
    to_label = gene %in% lab_genes
  )
  qcol <- c("concordant down (cell-autonomous)" = "#2166AC",
            "concordant up" = "#B2182B",
            "tissue-only (non-cell-autonomous)" = "#F1A340",
            "other / discordant" = "grey75")
  rr <- if (is.finite(rho_all)) sprintf("%.2f", rho_all) else "NA"
  p <- ggplot(m, aes(lfc_tis, lfc_cell)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    geom_point(aes(color = quadrant), size = 1.1, alpha = 0.6) +
    geom_point(data = dplyr::filter(m, to_label), color = "black", size = 1.8) +
    ggrepel::geom_text_repel(data = dplyr::filter(m, to_label),
                             aes(label = gene), size = 3, max.overlaps = 30,
                             segment.color = "grey50") +
    scale_color_manual(values = qcol, name = NULL) +
    labs(title = paste0("KAT8-KD log2FC concordance: cells vs ", tis_label),
         subtitle = paste0("Spearman(Wald) = ", rr,
                           "  |  down-left quadrant = direct/cell-autonomous candidates"),
         x = paste0(tis_label, " log2FC (KO vs CTL)"),
         y = "3T3-L1 cells log2FC (KAT8KD vs CTL)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTDIR, paste0("concordance_scatter_", tis_label, "_", RUN_TAG, ".png")),
         p, width = cp$plot_width, height = cp$plot_height, dpi = 300)

  ## partition table for downstream use
  write.csv(m, file.path(OUTDIR, paste0("concordance_genes_", tis_label, "_", RUN_TAG, ".csv")),
            row.names = FALSE)

  summary_rows[[tis_label]] <<- data.frame(
    comparison = paste0("cells_vs_", tis_label),
    n_shared = nrow(m), rho_all = rho_all, rho_tissueDEGs = rho_sig,
    n_tissueDEGs = nrow(tis_sig), pi1 = pi1,
    NES_tissueDOWN = if (!is.null(gsea) && "tissue_DOWN" %in% gsea$pathway) gsea$NES[gsea$pathway == "tissue_DOWN"] else NA,
    padj_tissueDOWN = if (!is.null(gsea) && "tissue_DOWN" %in% gsea$pathway) gsea$padj[gsea$pathway == "tissue_DOWN"] else NA,
    NES_tissueUP = if (!is.null(gsea) && "tissue_UP" %in% gsea$pathway) gsea$NES[gsea$pathway == "tissue_UP"] else NA,
    padj_tissueUP = if (!is.null(gsea) && "tissue_UP" %in% gsea$pathway) gsea$padj[gsea$pathway == "tissue_UP"] else NA
  )
  dplyr::bind_rows(r1, r2)
}

dir_tests <- dplyr::bind_rows(lapply(TISSUE_KEYS, function(k)
  analyse_pair(de[[CELL_KEY]], de[[k]], k)))

## 4) Write summaries ---------------------------------------
if (length(dir_tests) > 0)
  write.csv(dir_tests, file.path(OUTDIR, paste0("concordance_directional_tests_", RUN_TAG, ".csv")),
            row.names = FALSE)
summ <- dplyr::bind_rows(summary_rows)
write.csv(summ, file.path(OUTDIR, paste0("concordance_SUMMARY_", RUN_TAG, ".csv")), row.names = FALSE)

cat("\n=== CONCORDANCE SUMMARY ===\n"); print(summ)

## =========================================================
## 5) PROGRAM-LEVEL VIEW  (cells vs iWAT vs gWAT)
## =========================================================
## The gene-by-gene concordance above can hide the headline. This section asks
## the same question one level up: for each curated biological PROGRAM, how big
## is the KAT8 effect in the pure adipocytes versus each depot?
##
## Each program gets a COMPETITIVE test (Wilcoxon of its genes' Wald statistics
## against all other genes in that dataset), so "the program moved" is a
## statistical claim, not just a mean. This is what shows -- in one figure --
## that the adipocyte-identity / thermogenic / lipid programs collapse in iWAT
## but are FLAT in the cells, while ECM is the program that reproduces.
PROGRAMS <- list(
  "Adipocyte identity" = c("Pparg","Adipoq","Plin1","Cd36","Fabp4","Lep","Cfd","Lipe",
                           "Pnpla2","Cav1","Retn"),
  "Thermogenic"        = c("Ucp1","Cidea","Ppargc1a","Cox7a1","Elovl3","Dio2","Cox8b","Prdm16"),
  "Lipid metabolism"   = c("Ppara","Acox1","Acacb","Pck1","Scd1","Fasn","Acaca","Cpt1b","Lpl"),
  "Inflammatory"       = c("Ccl2","Ccl8","Cxcl2","Cxcl10","Il1b","Il12b","Itgax","Ctss",
                           "Adgre1","Tnf","Cd68","Lyz2"),
  "ECM / fibrosis"     = c("Col1a1","Col3a1","Col6a1","Fn1","Timp1","Sparc","Col1a2","Lox")
)

prog_rows <- list()
for (nm in names(de)) {
  d <- de[[nm]]
  for (pname in names(PROGRAMS)) {
    present <- intersect(PROGRAMS[[pname]], d$gene)
    if (length(present) < 3) next
    sub <- d[d$gene %in% present, ]
    bg  <- d[!d$gene %in% present, ]
    pv <- tryCatch(suppressWarnings(wilcox.test(sub$stat, bg$stat)$p.value),
                   error = function(e) NA_real_)
    prog_rows[[length(prog_rows) + 1]] <- data.frame(
      dataset      = nm,
      program      = pname,
      n_genes      = length(present),
      mean_log2FC  = mean(sub$log2FC, na.rm = TRUE),
      mean_Wald    = mean(sub$stat, na.rm = TRUE),
      n_sig        = sum(sub$padj < cp$tissue_fdr, na.rm = TRUE),
      p_vs_background = pv,
      stringsAsFactors = FALSE)
  }
}
prog <- dplyr::bind_rows(prog_rows)
if (nrow(prog) > 0) {
  write.csv(prog, file.path(OUTDIR, paste0("program_level_concordance_", RUN_TAG, ".csv")),
            row.names = FALSE)
  cat("\n=== PROGRAM-LEVEL SUMMARY (mean log2FC) ===\n")
  print(tidyr::pivot_wider(prog[, c("program","dataset","mean_log2FC")],
                           names_from = dataset, values_from = mean_log2FC),
        n = 100)

  ord <- c(CELL_KEY, TISSUE_KEYS)
  prog$dataset <- factor(prog$dataset, levels = ord[ord %in% unique(prog$dataset)])
  prog$program <- factor(prog$program, levels = rev(names(PROGRAMS)))
  p_prog <- ggplot(prog, aes(x = mean_log2FC, y = program, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             color = "black", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
    labs(x = expression(mean~log[2]~fold~change~(KD/KO~vs~CTL)), y = NULL,
         title = "Program-level KAT8 response: pure adipocytes vs tissue",
         subtitle = "flat in cells + strong in tissue = NOT cell-autonomous",
         caption = paste0("Curated programs; each tested competitively against all other genes ",
                          "in that dataset (Wilcoxon on Wald stats).\nSee program_level_concordance_",
                          RUN_TAG, ".csv for n, p-values and significant-gene counts.")) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(size = 7.5, colour = "grey35", hjust = 0))
  ggsave(file.path(OUTDIR, paste0("program_level_concordance_", RUN_TAG, ".png")),
         p_prog, width = 9, height = 4.8, dpi = 300, bg = "white")
  cat("[OK] Wrote program-level table + figure.\n")
}

cat("\n[DONE] Part 4b (concordance) complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Read-out:\n")
cat(" - High Spearman + positive NES(tissue_DOWN) in cells = the metabolic/adipocyte\n")
cat("   down-program is CELL-AUTONOMOUS (direct-candidate).\n")
cat(" - tissue_UP inflammatory set NOT enriched in cells = NON-cell-autonomous (indirect,\n")
cat("   consistent with the macrophage infiltration seen histologically).\n")
cat(" - Program-level figure: any program strong in tissue but FLAT in cells is\n")
cat("   non-cell-autonomous; a program that reproduces in cells is the direct candidate.\n")
