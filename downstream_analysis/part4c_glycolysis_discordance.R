#!/usr/bin/env Rscript

## -------------------------------------------------------------------------
## SCRIPT VERSION: 2026-07-27
##   glycolysis discordance: cells vs iWAT vs gWAT
##   All pipeline scripts should carry the SAME version date. run_all.R prints
##   them at pre-flight -- a date that differs from the rest means that file is
##   a stale copy and should be re-downloaded before you trust its output.
## -------------------------------------------------------------------------
## =========================================================
## KAT8 - PART 4c: GLYCOLYSIS DISCORDANCE
##                (pure adipocytes vs iWAT vs gWAT)
## =========================================================
## Part 4b asks the cell-autonomy question for whole programs. This script
## asks it for ONE pathway, at gene resolution, because glycolysis is the
## place the three datasets disagree most sharply -- and the disagreement is
## the result, not noise to be averaged away.
##
## Three claims are tested, each with a statistic rather than an eyeball:
##
##   (1) CELLS vs iWAT.   In the pure 3T3-L1 adipocytes the glycolytic
##       cassette is flat-to-UP under KAT8-KD. In iWAT the SAME genes fall
##       coherently. Since 3T3-L1 has no immune infiltrate and no stromal
##       compartment, a tissue-only collapse is not a direct KAT8
##       transcriptional effect -- it is secondary (glucose handling /
##       insulin sensitivity, or a shift in cell composition).
##
##   (2) iWAT vs gWAT.    The two depots also disagree WITH EACH OTHER.
##       Pgk1 is the cleanest case: significantly DOWN in iWAT and
##       significantly UP in gWAT. So there is no single "adipose tissue"
##       glycolytic response to compare the cells against.
##
##   (3) Is Hk2/Pfkl special, or is it the whole cascade? The competitive
##       Wilcoxon (below) tests the pathway as a unit against the rest of
##       the transcriptome in each dataset, so the answer does not rest on
##       whichever gene was looked at first.
##
## NOTE ON DIRECTION: KAT8 is an activator, so its DIRECT consequence should
## be down-regulation. Glycolysis going UP in the cells is therefore the
## opposite of a direct-target signature, which is the second independent
## reason to read the iWAT collapse as indirect.
##
## Everything here is threshold-free on the cell side (Wald statistics, no
## fold-change cut), so the low cell DEG count does not limit it.
## =========================================================

## 1) Packages ----------------------------------------------
required_cran <- c("dplyr", "tidyr", "ggplot2", "scales")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

## 2) Config -------------------------------------------------
source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
init_sanity("part4c_glycolysis")

de <- load_all_de()
check_gene_names(de)

DATASETS <- c(CELL_KEY, TISSUE_KEYS)

## 3) Gene sets ----------------------------------------------
## Ordered as the PATHWAY, not alphabetically: glucose in at the top,
## lactate out at the bottom. Reading the figure top-to-bottom is reading
## the cascade, so a break in the chain is visible as a break.
##
## GLYCOLYSIS_CORE is the ten-step backbone and is what the competitive
## test uses -- it is a textbook set, fixed in advance, with no gene chosen
## because of how it behaves in these data. The regulators/transporters are
## shown in the figure but tested separately, because PFKFB2/3/4 and the
## GLUTs are hormonally controlled and would otherwise smuggle an insulin-
## signalling effect into a claim about glycolytic enzymes.
## Aldoa is DELIBERATELY kept in this list even though it is absent from all
## three DE tables (Aldob and Aldoc are present in the tissue data; Aldoa is
## not, and the cells table contains no aldolase at all). Dropping it would
## hide that the canonical glycolytic aldolase is unmeasured here; leaving it
## in draws an explicit "n.d." row and logs it as a coverage gap. It
## contributes to no statistic either way.
GLYCOLYSIS_CORE <- c(
  "Hk1", "Hk2", "Hk3", "Gpi1",           ## glucose -> F6P
  "Pfkl", "Pfkm", "Pfkp",                ## committed step
  "Aldoa", "Tpi1",                       ## cleavage / isomerisation
  "Gapdh", "Pgk1", "Pgam1", "Eno1",      ## payoff phase
  "Pkm"                                  ## final step
)

GLYCOLYSIS_REGULATORY <- c(
  "Slc2a1", "Slc2a4",                    ## glucose uptake (GLUT1 / GLUT4)
  "Pfkfb1", "Pfkfb2", "Pfkfb3", "Pfkfb4",## F2,6BP -- allosteric PFK control
  "Ldha", "Ldhb"                         ## pyruvate <-> lactate
)

GLYCOLYSIS_ALL <- c(GLYCOLYSIS_CORE, GLYCOLYSIS_REGULATORY)

## Report coverage rather than silently dropping genes.
for (nm in DATASETS) {
  miss <- setdiff(GLYCOLYSIS_ALL, de[[nm]]$gene)
  log_sanity(paste0(nm, ": glycolytic genes measured"),
             sprintf("%d/%d", length(GLYCOLYSIS_ALL) - length(miss), length(GLYCOLYSIS_ALL)),
             ifelse(length(miss) > length(GLYCOLYSIS_ALL) / 3, "WARN", "OK"),
             if (length(miss)) paste0("not measured: ", paste(miss, collapse = ", ")) else "")
}

## 4) Competitive pathway test per dataset -------------------
## Same design as the program-level test in part4b: Wilcoxon of the set's
## Wald statistics against ALL OTHER genes in that dataset. Competitive (not
## self-contained), so "the pathway moved" means it moved relative to the
## transcriptome, which is the claim worth making.
competitive_test <- function(d, set_genes, set_label, dataset_label) {
  present <- intersect(set_genes, d$gene)
  sub <- d[d$gene %in% present & !is.na(d$stat), ]
  bg  <- d[!d$gene %in% present & !is.na(d$stat), ]
  if (nrow(sub) < 3) return(NULL)
  pv <- tryCatch(suppressWarnings(wilcox.test(sub$stat, bg$stat)$p.value),
                 error = function(e) NA_real_)
  data.frame(
    dataset      = dataset_label,
    gene_set     = set_label,
    n_genes      = nrow(sub),
    mean_log2FC  = mean(sub$log2FC, na.rm = TRUE),
    mean_Wald    = mean(sub$stat, na.rm = TRUE),
    n_up_sig     = sum(sub$padj < FDR_CUT & sub$log2FC > 0, na.rm = TRUE),
    n_down_sig   = sum(sub$padj < FDR_CUT & sub$log2FC < 0, na.rm = TRUE),
    p_vs_background = pv,
    stringsAsFactors = FALSE)
}

path_rows <- list()
for (nm in DATASETS) {
  path_rows[[length(path_rows) + 1]] <-
    competitive_test(de[[nm]], GLYCOLYSIS_CORE, "Glycolysis (core cascade)", nm)
  path_rows[[length(path_rows) + 1]] <-
    competitive_test(de[[nm]], GLYCOLYSIS_REGULATORY, "Glycolysis (regulators/transporters)", nm)
}
pathway <- dplyr::bind_rows(path_rows)

cat("\n=== COMPETITIVE PATHWAY TEST (Wilcoxon vs all other genes) ===\n")
print(pathway, row.names = FALSE, digits = 3)

for (i in seq_len(nrow(pathway))) {
  r <- pathway[i, ]
  log_sanity(sprintf("%s / %s: direction", r$dataset, r$gene_set),
             sprintf("mean Wald %+.2f (%s), p = %.2e",
                     r$mean_Wald, ifelse(r$mean_Wald < 0, "DOWN", "UP"), r$p_vs_background),
             ifelse(is.na(r$p_vs_background) || r$p_vs_background > 0.05, "INFO", "OK"))
}

write.csv(pathway, file.path(OUTDIR, paste0("glycolysis_pathway_tests_", RUN_TAG, ".csv")),
          row.names = FALSE)

## 5) Gene-level matrix across the three datasets ------------
## Long table first (one row per gene per dataset), because both the
## discordance tests and the figure want it in that shape.
gl <- dplyr::bind_rows(lapply(DATASETS, function(nm) {
  de[[nm]] %>%
    dplyr::filter(gene %in% GLYCOLYSIS_ALL) %>%
    dplyr::transmute(gene, dataset = nm, log2FC, stat, padj)
}))

gl <- gl %>% dplyr::mutate(
  block = ifelse(gene %in% GLYCOLYSIS_CORE, "Core cascade", "Regulators / transporters"),
  sig   = dplyr::case_when(
    is.na(padj)   ~ "",
    padj < 0.001  ~ "***",
    padj < 0.01   ~ "**",
    padj < 0.05   ~ "*",
    TRUE          ~ ""))

write.csv(gl %>% dplyr::arrange(block, gene, dataset),
          file.path(OUTDIR, paste0("glycolysis_gene_level_", RUN_TAG, ".csv")),
          row.names = FALSE)

## 6) Discordance tests --------------------------------------
## (a) SIGN DISCORDANCE between a pair of datasets, counted only over genes
##     that are significant in AT LEAST ONE of the two. Restricting to genes
##     with real signal somewhere is what stops the count being dominated by
##     genes whose sign is a coin flip because nothing happened to them.
## (b) PAIRED Wilcoxon on the Wald statistics of the same genes in the two
##     datasets. This is the direct test of "the pathway responds
##     differently here than there", and unlike (a) it uses effect size.
discordance <- function(a, b) {
  m <- dplyr::inner_join(
    de[[a]] %>% dplyr::filter(gene %in% GLYCOLYSIS_ALL) %>%
      dplyr::select(gene, stat_a = stat, lfc_a = log2FC, padj_a = padj),
    de[[b]] %>% dplyr::filter(gene %in% GLYCOLYSIS_ALL) %>%
      dplyr::select(gene, stat_b = stat, lfc_b = log2FC, padj_b = padj),
    by = "gene") %>%
    dplyr::filter(!is.na(stat_a), !is.na(stat_b))
  if (nrow(m) < 5) return(NULL)

  either_sig <- m %>% dplyr::filter((!is.na(padj_a) & padj_a < FDR_CUT) |
                                    (!is.na(padj_b) & padj_b < FDR_CUT))
  n_disc <- sum(sign(either_sig$lfc_a) != sign(either_sig$lfc_b))
  ## Both-significant AND opposite sign: the strongest, least arguable form.
  both_opposite <- m %>%
    dplyr::filter(!is.na(padj_a), !is.na(padj_b), padj_a < FDR_CUT, padj_b < FDR_CUT,
                  sign(lfc_a) != sign(lfc_b))

  p_paired <- tryCatch(
    suppressWarnings(wilcox.test(m$stat_a, m$stat_b, paired = TRUE)$p.value),
    error = function(e) NA_real_)
  rho <- suppressWarnings(cor(m$stat_a, m$stat_b, method = "spearman"))

  cat(sprintf("\n-- %s vs %s --\n", a, b))
  cat(sprintf("   shared glycolytic genes      : %d\n", nrow(m)))
  cat(sprintf("   sig in either, opposite sign : %d/%d (%.0f%%)\n",
              n_disc, nrow(either_sig),
              100 * n_disc / max(nrow(either_sig), 1)))
  cat(sprintf("   sig in BOTH, opposite sign   : %d%s\n", nrow(both_opposite),
              if (nrow(both_opposite)) paste0("  [", paste(both_opposite$gene, collapse = ", "), "]") else ""))
  cat(sprintf("   Spearman(Wald)               : %+.3f\n", rho))
  cat(sprintf("   paired Wilcoxon (Wald a vs b): p = %.2e\n", p_paired))

  data.frame(
    comparison        = paste0(a, "_vs_", b),
    n_shared          = nrow(m),
    n_sig_either      = nrow(either_sig),
    n_opposite_sign   = n_disc,
    frac_opposite     = n_disc / max(nrow(either_sig), 1),
    genes_sig_both_opposite = paste(both_opposite$gene, collapse = ";"),
    spearman_Wald     = rho,
    p_paired_wilcoxon = p_paired,
    stringsAsFactors  = FALSE)
}

cat("\n=== PAIRWISE DISCORDANCE (glycolytic genes only) ===\n")
pairs <- list(c(CELL_KEY, TISSUE_KEYS[1]),
              c(CELL_KEY, TISSUE_KEYS[2]),
              c(TISSUE_KEYS[1], TISSUE_KEYS[2]))
disc <- dplyr::bind_rows(lapply(pairs, function(p) discordance(p[1], p[2])))
write.csv(disc, file.path(OUTDIR, paste0("glycolysis_discordance_tests_", RUN_TAG, ".csv")),
          row.names = FALSE)

for (i in seq_len(nrow(disc))) {
  r <- disc[i, ]
  log_sanity(paste0(r$comparison, ": Spearman(Wald) on glycolysis"),
             sprintf("%+.3f | paired p = %.1e | %d/%d opposite sign",
                     r$spearman_Wald, r$p_paired_wilcoxon,
                     r$n_opposite_sign, r$n_sig_either),
             "INFO")
}

## 7) Figure --------------------------------------------------
## A dot-heatmap, not a bar chart: the claim is about the SIGN PATTERN across
## three columns, and colour reads that faster than length does. Genes stay in
## cascade order (top = glucose entry) and significance is printed IN the tile,
## so nothing about the figure requires the caption to interpret.
## A gene that was NOT MEASURED in a dataset must not render as an empty gap:
## an absent tile is white, and white on a diverging scale reads as "log2FC =
## 0", i.e. the exact opposite of "no data". Hk3 is missing from the cells and
## Aldoa from all three, so this is not hypothetical. The grid is completed
## explicitly and the blanks are drawn in grey and labelled "n.d.".
plot_df <- expand.grid(gene = GLYCOLYSIS_ALL, dataset = DATASETS,
                       stringsAsFactors = FALSE) %>%
  dplyr::left_join(gl, by = c("gene", "dataset")) %>%
  dplyr::mutate(
    measured = !is.na(log2FC),
    sig      = ifelse(is.na(sig), "", sig),
    label    = ifelse(measured, sig, "n.d."),
    block    = ifelse(gene %in% GLYCOLYSIS_CORE, "Core cascade", "Regulators / transporters"),
    gene     = factor(gene, levels = rev(GLYCOLYSIS_ALL)),
    dataset  = factor(dataset, levels = DATASETS),
    block    = factor(block, levels = c("Core cascade", "Regulators / transporters")))

n_nd <- sum(!plot_df$measured)
log_sanity("figure: gene x dataset cells with no measurement",
           sprintf("%d/%d drawn as n.d.", n_nd, nrow(plot_df)), "INFO",
           "grey tiles = not measured; NOT the same as log2FC = 0")

## Symmetric colour limits so "same distance from white" means "same effect
## size" in both directions. Clipped at the 99th percentile of |log2FC| so one
## extreme gene (Slc2a4) does not flatten every other tile to near-white.
lim <- as.numeric(quantile(abs(plot_df$log2FC), 0.99, na.rm = TRUE))
lim <- max(lim, 1)

p_heat <- ggplot(plot_df, aes(x = dataset, y = gene)) +
  ## grey underlay first, so a not-measured cell is visibly grey, not blank
  geom_tile(fill = "grey88", color = "white", linewidth = 0.6) +
  geom_tile(data = dplyr::filter(plot_df, measured),
            aes(fill = log2FC), color = "white", linewidth = 0.6) +
  geom_text(aes(label = label,
                fontface = ifelse(measured, "plain", "italic"),
                size = ifelse(measured, 3.6, 2.6)),
            vjust = 0.78, color = "grey15", show.legend = FALSE) +
  scale_size_identity() +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = col_down, mid = "white", high = col_up,
                       midpoint = 0, limits = c(-lim, lim),
                       oob = scales::squish,
                       name = expression(log[2]~FC)) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Glycolysis is discordant between pure adipocytes and adipose tissue",
    subtitle = paste0("KAT8 KD (3T3-L1) vs KAT8 KO (iWAT, gWAT)\n",
                      "core cascade UP in cells, DOWN in iWAT, UP in gWAT"),
    x = NULL, y = NULL,
    caption = paste0(
      "Tile = log2 fold change (KD/KO vs CTL); * padj<0.05, ** <0.01, *** <0.001 (DESeq2 Benjamini-Hochberg).\n",
      "Grey + 'n.d.' = not measured in that dataset, which is NOT the same as no change.\n",
      "Genes ordered as the pathway, glucose entry at top. Colour clipped at +/-", sprintf("%.1f", lim),
      " so one extreme gene does not flatten the rest.\n",
      "Statistics in glycolysis_pathway_tests_", RUN_TAG, ".csv and glycolysis_discordance_tests_", RUN_TAG, ".csv")) +
  theme_pub(11) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "italic"),
        legend.position = "right")

ggsave(file.path(OUTDIR, paste0("glycolysis_discordance_heatmap_", RUN_TAG, ".png")),
       p_heat, width = 7.2, height = 7.4, dpi = 300, bg = "white")

## Companion panel: the pathway-level summary the heatmap implies, with the
## competitive p-value printed on each bar so the figure carries its own test.
bar_df <- pathway %>%
  dplyr::mutate(dataset = factor(dataset, levels = DATASETS),
                lab = ifelse(is.na(p_vs_background), "",
                             ifelse(p_vs_background < 0.001, "p<0.001",
                                    sprintf("p=%.3f", p_vs_background))))

p_bar <- ggplot(bar_df, aes(x = dataset, y = mean_Wald, fill = mean_Wald > 0)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.68) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = lab,
                vjust = ifelse(mean_Wald > 0, -0.6, 1.5)), size = 3.1) +
  facet_wrap(~ gene_set, nrow = 1) +
  scale_fill_manual(values = c("TRUE" = col_up, "FALSE" = col_down), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.18)) +
  labs(title = "Glycolytic pathway response, tested as a unit",
       subtitle = "Mean DESeq2 Wald statistic; competitive Wilcoxon vs all other genes in that dataset",
       x = NULL, y = "mean Wald statistic (KD/KO vs CTL)",
       caption = "Opposite signs between cells and iWAT = the tissue effect is not reproduced in pure adipocytes.") +
  theme_pub(11) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(face = "bold"))

ggsave(file.path(OUTDIR, paste0("glycolysis_pathway_bars_", RUN_TAG, ".png")),
       p_bar, width = 7.6, height = 4.0, dpi = 300, bg = "white")

write_sanity()

cat("\n[DONE] Part 4c (glycolysis discordance) complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Read-out:\n")
cat(" - Cells flat/UP while iWAT is DOWN  = the glycolytic collapse is NOT cell-autonomous.\n")
cat(" - iWAT vs gWAT opposite signs       = there is no single adipose glycolytic response;\n")
cat("                                       depot must be stated in any claim about it.\n")
cat(" - KAT8 is an ACTIVATOR, so a direct target set should go DOWN in the cells.\n")
cat("   Glycolysis going UP there is the opposite of a direct-target signature.\n")
