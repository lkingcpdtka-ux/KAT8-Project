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
## Part 4b asks the cell-autonomy question for whole programs. This asks it
## for ONE pathway at gene resolution, because glycolysis is where the three
## datasets disagree most sharply -- and the disagreement is the result.
##
## WHAT IS ACTUALLY CLAIMED, in order of how well it survives scrutiny:
##
##   1. DISCORDANCE (the robust result). The glycolytic response in the pure
##      adipocytes differs from iWAT, and iWAT differs from gWAT, by more than
##      expression-matched random gene sets do (both p = 2e-04). Four genes are
##      significant in BOTH cells and iWAT with OPPOSITE signs, and five in
##      both iWAT and gWAT. Those per-gene calls carry their own FDR and depend
##      on no gene set at all, so they are the claim to lead with.
##
##   2. DIRECTION IN CELLS / gWAT (holds). The core cascade is UP in the cells
##      (p = 2e-04) and UP in gWAT (p = 0.009) against an expression-matched
##      null. KAT8 is an ACTIVATOR, so its direct consequence should be DOWN-
##      regulation; glycolysis going UP in the pure adipocytes is the opposite
##      of a direct-target signature.
##
##   3. DIRECTION IN iWAT (WEAKER THAN IT LOOKS -- read this before presenting).
##      Naively the cascade is clearly down in iWAT (Wilcoxon vs all other
##      genes, p = 0.011). But glycolytic genes sit around the 73rd-90th
##      expression percentile, and in iWAT highly expressed genes trend DOWN
##      generally (expression-matched null mean Wald = -0.65, not 0). Against
##      that matched null the cascade gives p = 0.054 -- borderline, not the
##      clean collapse the naive test implies. The REGULATORS/transporters
##      (GLUT4, PFKFB1/3) do survive matching (p = 0.006), so what is solid in
##      iWAT is the glucose-uptake/allosteric-control arm, not the enzymes.
##
## Hence both tests are computed and BOTH are reported. The naive Wilcoxon is
## kept only so the gap between the two is visible rather than hidden.
##
## GENE SET PROVENANCE: the primary set is GO:0006096 (glycolytic process)
## read offline from org.Mm.eg.db -- curated externally, with no knowledge of
## these data, so it cannot have been chosen to fit the result. The hand-
## ordered cascade below is used for FIGURE ROW ORDER and as a fallback if
## org.Mm.eg.db is unavailable; it is also tested separately so the two can be
## compared. If the conclusion held only for the curated list it would be a
## gene-picking artefact, and that is exactly what the comparison exposes.
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
N_PERM   <- 10000   ## permutations for the expression-matched null
N_BINS   <- 20      ## expression strata used for matching

## 3) Gene sets ----------------------------------------------
## Cascade order: glucose in at the top, pyruvate out at the bottom, so
## reading the figure top-to-bottom is reading the pathway.
##
## Aldoa is kept in this list even though it is absent from ALL THREE DE
## tables. That absence is real, not a naming variant: Aldob and Aldoc are
## present in the tissue data, no Ensembl-ID row corresponds to Aldoa, and the
## cells table contains no aldolase at all. Keeping it draws an explicit
## "n.d." row instead of silently omitting the canonical glycolytic aldolase.
## It contributes to no statistic either way.
CASCADE_ORDER <- c(
  "Hk1", "Hk2", "Hk3", "Gpi1",           ## glucose -> F6P
  "Pfkl", "Pfkm", "Pfkp",                ## committed step
  "Aldoa", "Tpi1",                       ## cleavage / isomerisation
  "Gapdh", "Pgk1", "Pgam1", "Eno1",      ## payoff phase
  "Pkm"                                  ## final step
)

## Kept SEPARATE from the enzymes on purpose. GLUT1/GLUT4 and PFKFB2/3/4 are
## hormonally controlled, and in iWAT this arm is dominated by Slc2a4
## (log2FC ~ -2.1). Pooling them would let an insulin-signalling effect be
## reported as a claim about glycolytic enzymes.
REGULATORY_ORDER <- c(
  "Slc2a1", "Slc2a4",                    ## glucose uptake (GLUT1 / GLUT4)
  "Pfkfb1", "Pfkfb2", "Pfkfb3", "Pfkfb4",## F2,6BP -- allosteric PFK control
  "Ldha", "Ldhb"                         ## pyruvate <-> lactate
)

FIGURE_GENES <- c(CASCADE_ORDER, REGULATORY_ORDER)

## Externally curated primary set (offline, via the go_genes() helper that
## config_downstream.R already provides for part5a).
go_glyc <- tryCatch(go_genes("GO:0006096"), error = function(e) character(0))
have_go <- length(go_glyc) >= 10
log_sanity("primary gene set: GO:0006096 (glycolytic process)",
           if (have_go) paste0(length(go_glyc), " genes from org.Mm.eg.db")
           else "UNAVAILABLE - falling back to curated cascade",
           if (have_go) "OK" else "WARN",
           if (have_go) "" else
             "install org.Mm.eg.db to remove reliance on the hand-curated list")

## The sets that get tested. GO first when present, so the headline number is
## never the hand-picked one.
TEST_SETS <- list()
if (have_go) TEST_SETS[["GO:0006096 glycolytic process"]] <- go_glyc
TEST_SETS[["Curated core cascade"]]   <- CASCADE_ORDER
TEST_SETS[["Curated regulators/transporters"]] <- REGULATORY_ORDER

for (nm in DATASETS) {
  miss <- setdiff(FIGURE_GENES, de[[nm]]$gene)
  log_sanity(paste0(nm, ": figure genes measured"),
             sprintf("%d/%d", length(FIGURE_GENES) - length(miss), length(FIGURE_GENES)),
             ifelse(length(miss) > length(FIGURE_GENES) / 3, "WARN", "OK"),
             if (length(miss)) paste0("not measured: ", paste(miss, collapse = ", ")) else "")
}

## 4) baseMean, for expression matching ----------------------
## config's load_de_table() does not carry baseMean, and the matched null
## needs it. Read it straight from the same CSVs, mirroring config's
## working-directory fallback so this behaves identically when the script is
## sourced from the project root or from inside downstream_analysis/.
resolve_de_path <- function(path) {
  if (file.exists(path)) return(path)
  alt <- c(file.path("downstream_analysis", "data", basename(path)),
           file.path("data", basename(path)),
           file.path("..", "downstream_analysis", "data", basename(path)))
  hit <- alt[file.exists(alt)]
  if (length(hit) > 0) hit[1] else NA_character_
}

load_basemean <- function(path) {
  p <- resolve_de_path(path)
  if (is.na(p)) return(NULL)
  df <- utils::read.csv(p, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"baseMean" %in% colnames(df)) return(NULL)
  gene <- if ("gene_name" %in% colnames(df)) as.character(df$gene_name) else rownames(df)
  out <- data.frame(gene = gene, baseMean = as.numeric(df$baseMean),
                    stringsAsFactors = FALSE)
  out[!duplicated(out$gene) & !is.na(out$baseMean), , drop = FALSE]
}

for (nm in DATASETS) {
  bmt <- load_basemean(DE_PATHS[[nm]])
  de[[nm]] <- if (is.null(bmt)) {
    de[[nm]]$baseMean <- NA_real_; de[[nm]]
  } else dplyr::left_join(de[[nm]], bmt, by = "gene")
}
have_bm <- vapply(DATASETS, function(nm) any(!is.na(de[[nm]]$baseMean)), logical(1))
log_sanity("baseMean available for expression matching",
           paste(DATASETS[have_bm], collapse = ", "),
           ifelse(all(have_bm), "OK", "WARN"),
           if (all(have_bm)) "" else "matched test skipped where absent; naive Wilcoxon still reported")

## 5) Two competitive tests ----------------------------------
## NAIVE: Wilcoxon of the set's Wald statistics against every other gene.
## This is what part4b's program test does, and it is CONFOUNDED here:
## glycolytic genes are highly expressed, highly expressed genes are better
## powered, and in iWAT they trend down as a class. So it is reported, but it
## is not the headline.
naive_wilcox <- function(d, genes) {
  present <- intersect(genes, d$gene)
  sub <- d[d$gene %in% present & !is.na(d$stat), ]
  bg  <- d[!d$gene %in% present & !is.na(d$stat), ]
  if (nrow(sub) < 3) return(NA_real_)
  tryCatch(suppressWarnings(wilcox.test(sub$stat, bg$stat)$p.value),
           error = function(e) NA_real_)
}

## MATCHED: compare the set's mean Wald statistic against random gene sets
## drawn to have the SAME expression profile (same baseMean stratum counts).
## The null is then "genes as well measured as these ones", which is the
## comparison the naive test should have been making all along.
matched_perm <- function(d, genes, n_perm = N_PERM, n_bins = N_BINS) {
  d <- d[!is.na(d$stat) & !is.na(d$baseMean), ]
  present <- intersect(genes, d$gene)
  if (length(present) < 3 || nrow(d) < 200)
    return(list(obs = NA_real_, null_mean = NA_real_, p = NA_real_, n = length(present)))
  d <- d[order(d$baseMean), ]
  bin <- pmin(seq_len(nrow(d)) %/% max(1, floor(nrow(d) / n_bins)) + 1, n_bins)
  names(bin) <- d$gene
  stat <- setNames(d$stat, d$gene)
  obs <- mean(stat[present])
  want <- table(bin[present])
  pool <- split(d$gene, bin)
  null <- replicate(n_perm, {
    s <- unlist(lapply(names(want), function(b) {
      k <- want[[b]]; cand <- pool[[b]]
      if (length(cand) >= k) sample(cand, k) else sample(cand, k, replace = TRUE)
    }), use.names = FALSE)
    mean(stat[s])
  })
  p_one <- (sum(if (obs >= 0) null >= obs else null <= obs) + 1) / (n_perm + 1)
  list(obs = obs, null_mean = mean(null), p = min(1, 2 * p_one), n = length(present))
}

rows <- list()
for (sname in names(TEST_SETS)) {
  for (nm in DATASETS) {
    d <- de[[nm]]
    present <- intersect(TEST_SETS[[sname]], d$gene)
    if (length(present) < 3) next
    sub <- d[d$gene %in% present, ]
    mp  <- matched_perm(d, TEST_SETS[[sname]])
    rows[[length(rows) + 1]] <- data.frame(
      gene_set        = sname,
      dataset         = nm,
      n_genes         = length(present),
      mean_log2FC     = mean(sub$log2FC, na.rm = TRUE),
      mean_Wald       = mean(sub$stat, na.rm = TRUE),
      n_up_sig        = sum(sub$padj < FDR_CUT & sub$log2FC > 0, na.rm = TRUE),
      n_down_sig      = sum(sub$padj < FDR_CUT & sub$log2FC < 0, na.rm = TRUE),
      p_naive_wilcox  = naive_wilcox(d, TEST_SETS[[sname]]),
      matched_null_mean = mp$null_mean,
      p_matched       = mp$p,
      stringsAsFactors = FALSE)
  }
}
pathway <- dplyr::bind_rows(rows)

cat("\n=== PATHWAY TESTS: naive vs expression-matched ===\n")
cat("(p_naive compares to ALL other genes; p_matched compares to random sets\n")
cat(" with the same expression profile -- the honest one)\n\n")
print(pathway, row.names = FALSE, digits = 3)

## Flag any set where matching changes the verdict at 0.05. That disagreement
## is a finding about the METHOD and belongs in the sanity ledger, not buried.
flip <- pathway %>% dplyr::filter(!is.na(p_naive_wilcox), !is.na(p_matched),
                                  p_naive_wilcox < 0.05, p_matched >= 0.05)
if (nrow(flip) > 0) {
  for (i in seq_len(nrow(flip))) {
    r <- flip[i, ]
    log_sanity(paste0(r$dataset, " / ", r$gene_set, ": significance lost on matching"),
               sprintf("naive p=%.3f -> matched p=%.3f", r$p_naive_wilcox, r$p_matched),
               "WARN",
               "driven by these genes being highly expressed, not by the pathway; do not present the naive p")
  }
} else {
  log_sanity("sets losing significance once expression-matched", "none", "OK")
}

write.csv(pathway, file.path(OUTDIR, paste0("glycolysis_pathway_tests_", RUN_TAG, ".csv")),
          row.names = FALSE)

## 6) Leave-one-out robustness -------------------------------
## "Did one gene carry the whole set?" Recompute the set mean dropping each
## gene in turn; if no single drop flips the sign, no single gene is the result.
loo_rows <- list()
for (sname in names(TEST_SETS)) {
  for (nm in DATASETS) {
    d <- de[[nm]]; present <- intersect(TEST_SETS[[sname]], d$gene)
    st <- d$stat[match(present, d$gene)]; st <- st[!is.na(st)]
    if (length(st) < 4) next
    full <- mean(st)
    loo  <- vapply(seq_along(st), function(i) mean(st[-i]), numeric(1))
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      gene_set = sname, dataset = nm, mean_Wald = full,
      loo_min = min(loo), loo_max = max(loo),
      any_sign_flip = any((loo >= 0) != (full >= 0)),
      stringsAsFactors = FALSE)
  }
}
loo <- dplyr::bind_rows(loo_rows)
cat("\n=== LEAVE-ONE-OUT (does a single gene carry the set?) ===\n")
print(loo, row.names = FALSE, digits = 3)
log_sanity("leave-one-out sign flips across all sets/datasets",
           sum(loo$any_sign_flip), ifelse(any(loo$any_sign_flip), "WARN", "OK"),
           ifelse(any(loo$any_sign_flip),
                  "a single gene flips a set's direction -- that set is not a program-level result", ""))
write.csv(loo, file.path(OUTDIR, paste0("glycolysis_leave_one_out_", RUN_TAG, ".csv")),
          row.names = FALSE)

## 7) Gene-level table ---------------------------------------
gl <- dplyr::bind_rows(lapply(DATASETS, function(nm) {
  de[[nm]] %>%
    dplyr::filter(gene %in% FIGURE_GENES) %>%
    dplyr::transmute(gene, dataset = nm, log2FC, stat, padj)
}))
gl <- gl %>% dplyr::mutate(
  block = ifelse(gene %in% CASCADE_ORDER, "Core cascade", "Regulators / transporters"),
  sig   = dplyr::case_when(is.na(padj) ~ "", padj < 0.001 ~ "***",
                           padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))
write.csv(gl %>% dplyr::arrange(block, gene, dataset),
          file.path(OUTDIR, paste0("glycolysis_gene_level_", RUN_TAG, ".csv")),
          row.names = FALSE)

## 8) Discordance between datasets ---------------------------
## The per-gene "significant in BOTH, opposite sign" count is the part of this
## script that needs no gene set and no set-level null: each gene carries its
## own FDR in each dataset. It is the claim to lead with.
##
## The set-level version uses the SAME expression-matched permutation as
## above, applied to the per-gene DIFFERENCE in Wald statistic between the two
## datasets, so "these genes disagree more than comparable genes do" is tested
## rather than asserted.
disc_matched <- function(a, b, genes, n_perm = N_PERM, n_bins = N_BINS) {
  m <- dplyr::inner_join(
    de[[a]] %>% dplyr::select(gene, sa = stat, bm = baseMean),
    de[[b]] %>% dplyr::select(gene, sb = stat),
    by = "gene") %>% dplyr::filter(!is.na(sa), !is.na(sb), !is.na(bm))
  present <- intersect(genes, m$gene)
  if (length(present) < 3 || nrow(m) < 200) return(list(obs = NA_real_, null_mean = NA_real_, p = NA_real_))
  m <- m[order(m$bm), ]
  bin <- pmin(seq_len(nrow(m)) %/% max(1, floor(nrow(m) / n_bins)) + 1, n_bins)
  names(bin) <- m$gene
  dif <- setNames(m$sa - m$sb, m$gene)
  obs <- mean(dif[present])
  want <- table(bin[present]); pool <- split(m$gene, bin)
  null <- replicate(n_perm, {
    s <- unlist(lapply(names(want), function(bb) {
      k <- want[[bb]]; cand <- pool[[bb]]
      if (length(cand) >= k) sample(cand, k) else sample(cand, k, replace = TRUE)
    }), use.names = FALSE)
    mean(dif[s])
  })
  p_one <- (sum(if (obs >= 0) null >= obs else null <= obs) + 1) / (n_perm + 1)
  list(obs = obs, null_mean = mean(null), p = min(1, 2 * p_one))
}

discordance <- function(a, b) {
  m <- dplyr::inner_join(
    de[[a]] %>% dplyr::filter(gene %in% FIGURE_GENES) %>%
      dplyr::select(gene, sa = stat, la = log2FC, pa = padj),
    de[[b]] %>% dplyr::filter(gene %in% FIGURE_GENES) %>%
      dplyr::select(gene, sb = stat, lb = log2FC, pb = padj),
    by = "gene") %>% dplyr::filter(!is.na(sa), !is.na(sb))
  if (nrow(m) < 5) return(NULL)

  either <- m %>% dplyr::filter((!is.na(pa) & pa < FDR_CUT) | (!is.na(pb) & pb < FDR_CUT))
  n_disc <- sum(sign(either$la) != sign(either$lb))
  both_opp <- m %>% dplyr::filter(!is.na(pa), !is.na(pb), pa < FDR_CUT, pb < FDR_CUT,
                                  sign(la) != sign(lb))
  rho <- suppressWarnings(cor(m$sa, m$sb, method = "spearman"))
  mp  <- disc_matched(a, b, FIGURE_GENES)

  cat(sprintf("\n-- %s vs %s --\n", a, b))
  cat(sprintf("   sig in BOTH, opposite sign  : %d%s\n", nrow(both_opp),
              if (nrow(both_opp)) paste0("  [", paste(both_opp$gene, collapse = ", "), "]") else ""))
  cat(sprintf("   sig in either, opposite sign: %d/%d\n", n_disc, nrow(either)))
  cat(sprintf("   Spearman(Wald)              : %+.3f\n", rho))
  cat(sprintf("   mean Wald difference        : %+.2f (matched null %+.2f), p = %.2e\n",
              mp$obs, mp$null_mean, mp$p))

  data.frame(comparison = paste0(a, "_vs_", b), n_shared = nrow(m),
             n_sig_either = nrow(either), n_opposite_sign = n_disc,
             genes_sig_both_opposite = paste(both_opp$gene, collapse = ";"),
             spearman_Wald = rho, mean_Wald_difference = mp$obs,
             matched_null_mean = mp$null_mean, p_matched = mp$p,
             stringsAsFactors = FALSE)
}

cat("\n=== DISCORDANCE BETWEEN DATASETS ===\n")
disc <- dplyr::bind_rows(lapply(
  list(c(CELL_KEY, TISSUE_KEYS[1]), c(CELL_KEY, TISSUE_KEYS[2]), c(TISSUE_KEYS[1], TISSUE_KEYS[2])),
  function(p) discordance(p[1], p[2])))
write.csv(disc, file.path(OUTDIR, paste0("glycolysis_discordance_tests_", RUN_TAG, ".csv")),
          row.names = FALSE)
for (i in seq_len(nrow(disc))) {
  r <- disc[i, ]
  log_sanity(paste0(r$comparison, ": glycolytic discordance"),
             sprintf("%d genes sig-both-opposite | matched p = %.1e",
                     lengths(regmatches(r$genes_sig_both_opposite,
                                        gregexpr(";", r$genes_sig_both_opposite))) +
                       ifelse(nzchar(r$genes_sig_both_opposite), 1, 0),
                     r$p_matched), "INFO")
}

## 9) Figure --------------------------------------------------
## A gene that was NOT MEASURED must not render as an empty gap: an absent
## tile is white, and white on a diverging scale reads as log2FC = 0, the
## opposite of "no data". Grid is completed explicitly and blanks drawn grey.
plot_df <- expand.grid(gene = FIGURE_GENES, dataset = DATASETS, stringsAsFactors = FALSE) %>%
  dplyr::left_join(gl, by = c("gene", "dataset")) %>%
  dplyr::mutate(
    measured = !is.na(log2FC),
    label    = ifelse(measured, ifelse(is.na(sig), "", sig), "n.d."),
    block    = ifelse(gene %in% CASCADE_ORDER, "Core cascade", "Regulators / transporters"),
    gene     = factor(gene, levels = rev(FIGURE_GENES)),
    dataset  = factor(dataset, levels = DATASETS),
    block    = factor(block, levels = c("Core cascade", "Regulators / transporters")))

## Symmetric limits so equal distance from white means equal effect size,
## clipped at the 99th percentile so Slc2a4 does not flatten every other tile.
lim <- max(as.numeric(quantile(abs(plot_df$log2FC), 0.99, na.rm = TRUE)), 1)

p_heat <- ggplot(plot_df, aes(x = dataset, y = gene)) +
  geom_tile(fill = "grey88", color = "white", linewidth = 0.6) +
  geom_tile(data = dplyr::filter(plot_df, measured),
            aes(fill = log2FC), color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, size = ifelse(measured, 3.5, 2.5)),
            vjust = 0.78, color = "grey15", show.legend = FALSE) +
  scale_size_identity() +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = col_down, mid = "white", high = col_up, midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish,
                       name = expression(log[2]~FC)) +
  scale_x_discrete(position = "top") +
  labs(title = "Glycolytic gene expression after KAT8 loss",
       subtitle = "3T3-L1 adipocytes (KD) and two adipose depots (KO), vs their controls",
       x = NULL, y = NULL,
       caption = paste0("* padj<0.05  ** <0.01  *** <0.001.   Grey = not measured.   ",
                        "Rows in pathway order.   Colour clipped at +/-", sprintf("%.1f", lim), ".")) +
  theme_pub(11) +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey92", color = NA),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "italic"),
        legend.position = "right")

ggsave(file.path(OUTDIR, paste0("glycolysis_discordance_heatmap_", RUN_TAG, ".png")),
       p_heat, width = 7.0, height = 7.2, dpi = 300, bg = "white")

## Companion: pathway-level effect with BOTH p-values on each bar, so the
## place where matching changes the answer is visible in the figure itself.
bar_df <- pathway %>%
  dplyr::mutate(dataset = factor(dataset, levels = DATASETS),
                gene_set = factor(gene_set, levels = names(TEST_SETS)),
                lab = ifelse(is.na(p_matched), "matched p = NA",
                             ifelse(p_matched < 0.001, "matched p < 0.001",
                                    sprintf("matched p = %.3f", p_matched))))

p_bar <- ggplot(bar_df, aes(x = dataset, y = mean_Wald, fill = mean_Wald > 0)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.68) +
  geom_point(aes(y = matched_null_mean), shape = 18, size = 2.6, colour = "grey25") +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = lab, vjust = ifelse(mean_Wald > 0, -0.7, 1.6)), size = 2.7) +
  facet_wrap(~ gene_set, nrow = 1, labeller = label_wrap_gen(width = 26)) +
  scale_fill_manual(values = c("TRUE" = col_up, "FALSE" = col_down), guide = "none") +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  labs(title = "Glycolytic response tested as a set",
       subtitle = "Bar = mean Wald statistic; diamond = mean of expression-matched random sets",
       x = NULL, y = "mean Wald statistic",
       caption = paste0("p from 10,000 permutations of gene sets matched on expression level. ",
                        "The diamond is the null the bar must beat --\nwhere it sits away from zero, ",
                        "a test against all genes would have overstated the effect.")) +
  theme_pub(11) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text = element_text(face = "bold", size = 8),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        axis.text.x = element_text(face = "bold"))

ggsave(file.path(OUTDIR, paste0("glycolysis_pathway_bars_", RUN_TAG, ".png")),
       p_bar, width = 8.4, height = 4.2, dpi = 300, bg = "white")

write_sanity()

cat("\n[DONE] Part 4c complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Read-out:\n")
cat(" - Lead with the per-gene 'significant in BOTH, opposite sign' genes: no gene set,\n")
cat("   no set-level null, each gene carries its own FDR.\n")
cat(" - Cells UP while iWAT is DOWN = the glycolytic response is NOT cell-autonomous.\n")
cat(" - iWAT vs gWAT also disagree, so 'adipose tissue' is not one response; name the depot.\n")
cat(" - Check p_matched, never p_naive_wilcox, before claiming a set moved.\n")
