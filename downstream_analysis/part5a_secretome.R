#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   secretome: what the KAT8-null adipocyte broadcasts
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## PART 5a - SECRETOME
## Question: what does the KAT8-null adipocyte BROADCAST?
## =========================================================
## Why this matters
##   Part 4b showed the dominant tissue signature is NOT cell-autonomous: in
##   KAT8-KD 3T3-L1 the adipocyte-marker / lipid / thermogenic programs are
##   flat with tight confidence intervals, even though the knockdown is
##   stronger than in the tissue. The published result -- KAT8 is required only
##   when lost BEFORE differentiation -- explains it: these cells are a
##   MAINTENANCE experiment, and the mature adipocyte programme does not need
##   ongoing KAT8.
##
##   So if the KAT8-null adipocyte still looks and behaves like an adipocyte,
##   how does the tissue end
##   up inflamed and fibrotic? The obvious candidate is what the cell SENDS to
##   its neighbours. This script finds the secreted/surface genes that change
##   in the cells and tests whether each one reproduces in vivo.
##
## Read-out
##   Genes changed in cells AND moving the same direction in a depot are the
##   candidate non-cell-autonomous signals -- the mechanistic bridge between an
##   apparently healthy KAT8-null adipocyte and a diseased tissue.
##
## Input : DE CSVs configured in config_downstream.R      Output: OUTDIR
## =========================================================

required_cran <- c("dplyr","tidyr","ggplot2")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (bp in c("org.Mm.eg.db","AnnotationDbi","GO.db"))
  if (!requireNamespace(bp, quietly = TRUE))
    tryCatch(BiocManager::install(bp, ask = FALSE, update = FALSE), error = function(e) NULL)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(org.Mm.eg.db); library(AnnotationDbi); library(GO.db)
})

source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
init_sanity("part5a_secretome")

## ---- 1) Data ----
de <- load_all_de()
check_gene_names(de)
cells   <- de[[CELL_KEY]]
tissues <- intersect(TISSUE_KEYS, names(de))
cell_sig <- cells %>% dplyr::filter(!is.na(padj), padj < FDR_CUT)
log_sanity("cell DEGs at FDR<0.05 (no logFC gate)", nrow(cell_sig), "OK")
log_sanity("depots available for replication test",
           ifelse(length(tissues) == 0, "NONE", paste(tissues, collapse = ", ")),
           ifelse(length(tissues) == 0, "WARN", "OK"))

## ---- 2) Annotate by GO class (offline, non-circular) ----
secreted <- unique(c(go_genes(GO_IDS[["extracellular"]]), go_genes(GO_IDS[["ecm"]])))
ligands  <- unique(c(go_genes(GO_IDS[["ligand"]]), go_genes(GO_IDS[["cytokine"]]),
                     go_genes(GO_IDS[["growth_factor"]])))
surface  <- go_genes(GO_IDS[["cell_surface"]])
log_sanity("GO secreted/ECM genes", length(secreted), ifelse(length(secreted) > 500, "OK", "WARN"))
log_sanity("GO ligand genes",       length(ligands),  ifelse(length(ligands)  > 100, "OK", "WARN"))

sec <- cell_sig %>%
  dplyr::mutate(class = dplyr::case_when(
    gene %in% ligands  ~ "Secreted ligand",
    gene %in% secreted ~ "Secreted / ECM",
    gene %in% surface  ~ "Cell surface",
    TRUE ~ "Other")) %>%
  dplyr::filter(class != "Other") %>%
  dplyr::arrange(dplyr::desc(log2FC))
log_sanity("communication-class DEGs",
           paste0(nrow(sec), " / ", nrow(cell_sig)), "OK")

## ---- 3) Does each signal replicate in vivo? ----
for (t in tissues) {
  tt <- de[[t]] %>% dplyr::select(gene, l = log2FC, p = padj)
  names(tt) <- c("gene", paste0(t, "_log2FC"), paste0(t, "_padj"))
  sec <- dplyr::left_join(sec, tt, by = "gene")
}
dir_cols <- paste0(tissues, "_log2FC"); pad_cols <- paste0(tissues, "_padj")
if (length(tissues) > 0 && nrow(sec) > 0) {
  same <- sign(as.matrix(sec[, dir_cols, drop = FALSE])) == sign(sec$log2FC)
  sig  <- as.matrix(sec[, pad_cols, drop = FALSE]) < FDR_CUT
  sec$n_depot_same_dir     <- rowSums(same, na.rm = TRUE)
  sec$n_depot_sig_same_dir <- rowSums(same & sig, na.rm = TRUE)
} else {
  sec$n_depot_same_dir <- 0L; sec$n_depot_sig_same_dir <- 0L
}
write.csv(sec, file.path(OUTDIR, paste0("part5a_secretome_annotated_", RUN_TAG, ".csv")),
          row.names = FALSE)

## Which direction replicates better? Let the data decide, rather than assuming
## the KAT8-null cell "gains" a signal or "loses" one.
up_s <- sec %>% dplyr::filter(log2FC > 0); dn_s <- sec %>% dplyr::filter(log2FC < 0)
log_sanity("UP-secreted replicating in >=1 depot",
           sprintf("%d/%d", sum(up_s$n_depot_sig_same_dir > 0), nrow(up_s)), "INFO")
log_sanity("DOWN-secreted replicating in >=1 depot",
           sprintf("%d/%d", sum(dn_s$n_depot_sig_same_dir > 0), nrow(dn_s)), "INFO")
if (nrow(up_s) > 2 && nrow(dn_s) > 2) {
  ct <- matrix(c(sum(up_s$n_depot_sig_same_dir > 0), sum(up_s$n_depot_sig_same_dir == 0),
                 sum(dn_s$n_depot_sig_same_dir > 0), sum(dn_s$n_depot_sig_same_dir == 0)),
               nrow = 2, byrow = TRUE)
  ft <- fisher.test(ct)
  log_sanity("Fisher: UP vs DOWN replication asymmetry",
             sprintf("OR=%.2f p=%.3g", ft$estimate, ft$p.value), "INFO")
}

## ---- 4) Figure: the replicating secreted core ----
core <- sec %>% dplyr::filter(n_depot_sig_same_dir > 0) %>%
  dplyr::arrange(dplyr::desc(log2FC)) %>% dplyr::slice_head(n = 40)
log_sanity("secreted genes replicating in >=1 depot", nrow(core),
           ifelse(nrow(core) > 0, "OK", "WARN"))
if (nrow(core) >= 3) {
  pd <- core %>%
    dplyr::select(gene, cells = log2FC, dplyr::all_of(dir_cols)) %>%
    tidyr::pivot_longer(-gene, names_to = "dataset", values_to = "log2FC") %>%
    dplyr::mutate(dataset = sub("_log2FC$", "", dataset),
                  gene = factor(gene, levels = rev(core$gene)))
  p <- ggplot(pd, aes(log2FC, gene, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75,
             color = "black", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
    labs(x = expression(log[2]~fold~change~(KD/KO~vs~CTL)), y = NULL,
         title = "Secreted / surface DEGs that replicate in vivo",
         caption = paste0("Cell DEGs (FDR<", FDR_CUT, ") in GO extracellular / ECM / ligand / surface classes ",
                          "that move the SAME direction in >=1 depot at FDR<", FDR_CUT, ".\n",
                          "Candidate non-cell-autonomous mechanism: what the KAT8-null adipocyte sends to its ",
                          "neighbours. n=", nrow(core), " shown.")) +
    theme_pub()
  ggsave(file.path(OUTDIR, paste0("part5a_secretome_replicating_", RUN_TAG, ".png")),
         p, width = 9, height = max(5, nrow(core) * 0.26), dpi = 300, bg = "white")
  cat("[OK] wrote replicating-secretome figure\n")
}

write_sanity()
cat("\n[DONE] Part 5a complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Next: part5b_tf_analysis.R (which TF drives this program?)\n")
