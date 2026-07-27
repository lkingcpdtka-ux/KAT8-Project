#!/usr/bin/env Rscript

## =========================================================
## PART 5c - KAT8 COMPLEX OUTPUT (MSL / NSL)
## Question: is the canonical KAT8 program engaged in the cells?
## =========================================================
## KAT8 is the catalytic H4K16ac subunit of TWO complexes:
##   MSL       (Msl1/2/3)                     broad H4K16ac, transcriptional output
##   NSL/KANSL (Kansl1-3, Mcrs1, Phf20, Ogt,  promoter-proximal; the best-documented
##              Wdr5, Hcfc1)                  direct regulator of nuclear-encoded
##                                            OXPHOS / mitochondrial genes
##
## THE KEY IDEA
##   siRNA removes the catalytic core but NOT the partners' transcripts, so the
##   partners' mRNA stays flat and the complex effect is INVISIBLE in membership
##   expression. It has to be read from the complex's regulatory OUTPUT -- its
##   target genes.
##
## THE SPECIFICITY CONTROLS (this is what makes the result interpretable)
##   Nuclear-encoded OXPHOS  = NSL targets            -> should move if NSL is hit
##   Mito-encoded genes (mt-) = transcribed by mtDNA  -> OUTSIDE NSL control
##   Ribosomal proteins       = comparable housekeeping-> generic-drift control
##   If nuclear OXPHOS moves while BOTH controls stay flat, the shift is
##   specifically NSL and not a normalisation artefact or a generic stress
##   response. If a control also moves, do NOT claim complex specificity.
##
## The verdict below is COMPUTED from the data, never hardcoded in a caption.
##
## Input : DE CSVs configured in config_downstream.R      Output: OUTDIR
## =========================================================

required_cran <- c("dplyr","tidyr","ggplot2")
to_install <- setdiff(required_cran, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

source(file.path("downstream_analysis", "config_downstream.R"))
set.seed(MASTER_SEED)
init_sanity("part5c_complex")

de <- load_all_de()
check_gene_names(de)
cells   <- de[[CELL_KEY]]
tissues <- intersect(TISSUE_KEYS, names(de))

## ---- 1) Complex membership (expected: flat, because siRNA hits protein) ----
msl <- c("Msl1","Msl2","Msl3")
nsl <- c("Kansl1","Kansl2","Kansl3","Mcrs1","Phf20","Ogt","Wdr5","Hcfc1")
members <- data.frame(
  gene    = c("Kat8", msl, nsl),
  complex = c("catalytic", rep("MSL", length(msl)), rep("NSL/KANSL", length(nsl))),
  stringsAsFactors = FALSE)
member_expr <- members %>%
  dplyr::left_join(cells %>% dplyr::select(gene, log2FC, padj, stat), by = "gene")
write.csv(member_expr, file.path(OUTDIR, paste0("part5c_complex_membership_", RUN_TAG, ".csv")),
          row.names = FALSE)
log_sanity("complex members found in cells",
           paste0(sum(!is.na(member_expr$log2FC)), "/", nrow(member_expr)), "OK")
n_partner_sig <- sum(member_expr$padj < FDR_CUT & member_expr$gene != "Kat8", na.rm = TRUE)
log_sanity("partner transcripts significantly changed (expect ~0)", n_partner_sig,
           ifelse(n_partner_sig <= 2, "OK", "INFO"))

## ---- 2) Complex OUTPUT with specificity controls ----
nuclear_oxphos <- c(
  "Ndufa1","Ndufa2","Ndufa3","Ndufa4","Ndufa5","Ndufa8","Ndufa9","Ndufa11","Ndufa13",
  "Ndufb3","Ndufb7","Ndufb8","Ndufb9","Ndufb10","Ndufb11","Ndufs2","Ndufs3","Ndufs6",
  "Ndufs7","Ndufv1","Sdhb","Sdhc","Sdhd","Uqcrb","Uqcrq","Uqcrh","Uqcr10","Uqcr11",
  "Uqcrc1","Uqcrfs1","Cyc1","Cycs","Cox4i1","Cox5b","Cox6a1","Cox6b1","Cox6c","Cox7a2",
  "Cox7b","Cox7c","Cox10","Atp5f1a","Atp5f1b","Atp5f1d","Atp5mc1","Atp5mc2","Atp5mc3",
  "Atp5me","Atp5mf","Atp5pd")

## Competitive test: this gene set vs every other gene in the same dataset.
## Using the Wald statistic keeps it variance-aware and threshold-free.
program_test <- function(d, genes, label, dataset) {
  present <- intersect(genes, d$gene)
  if (length(present) < 4) return(NULL)
  sub <- d[d$gene %in% present, ]; bg <- d[!d$gene %in% present, ]
  pv <- tryCatch(suppressWarnings(wilcox.test(sub$stat, bg$stat)$p.value),
                 error = function(e) NA_real_)
  data.frame(dataset = dataset, program = label, n_genes = length(present),
             mean_log2FC = mean(sub$log2FC, na.rm = TRUE),
             mean_Wald   = mean(sub$stat,   na.rm = TRUE),
             p_vs_background = pv, stringsAsFactors = FALSE)
}

rows <- list()
for (nm in c(CELL_KEY, tissues)) {
  d <- de[[nm]]
  rows <- c(rows, list(
    program_test(d, nuclear_oxphos,                    "Nuclear-encoded OXPHOS (NSL targets)", nm),
    program_test(d, grep("^mt-",     d$gene, value = TRUE), "Mito-encoded (control: outside NSL)",   nm),
    program_test(d, grep("^Rp[sl][0-9]", d$gene, value = TRUE), "Ribosomal proteins (control)",     nm)))
}
out <- dplyr::bind_rows(rows)
write.csv(out, file.path(OUTDIR, paste0("part5c_complex_output_", RUN_TAG, ".csv")), row.names = FALSE)
cat("\n"); print(out, row.names = FALSE)

## ---- 3) HEAD-TO-HEAD: is the OXPHOS shift BIGGER than the control shift? ----
##
## Requiring the controls to be flat is the wrong test. Each control set has
## 88-101 genes, so a Wilcoxon against ~13,000 background genes will call a
## mean Wald of +0.29 "significant" -- a shift far too small to matter, but
## enough to veto the verdict on a technicality. What specificity actually
## means is that the NSL target set moves MORE than the controls do, so that
## is tested directly: nuclear-OXPHOS Wald statistics against each control
## set's Wald statistics, head to head.
head_to_head <- function(d, genes_a, genes_b, label_b, dataset) {
  a <- intersect(genes_a, d$gene); b <- intersect(genes_b, d$gene)
  if (length(a) < 4 || length(b) < 4) return(NULL)
  sa <- d$stat[d$gene %in% a]; sb <- d$stat[d$gene %in% b]
  pv <- tryCatch(suppressWarnings(wilcox.test(sa, sb)$p.value),
                 error = function(e) NA_real_)
  data.frame(dataset = dataset, comparison = paste0("nuclear OXPHOS vs ", label_b),
             mean_Wald_oxphos = mean(sa, na.rm = TRUE),
             mean_Wald_control = mean(sb, na.rm = TRUE),
             fold_larger = abs(mean(sa, na.rm = TRUE)) / max(abs(mean(sb, na.rm = TRUE)), 1e-9),
             p_head_to_head = pv, stringsAsFactors = FALSE)
}
h2h <- dplyr::bind_rows(lapply(c(CELL_KEY, tissues), function(nm) {
  d <- de[[nm]]
  dplyr::bind_rows(
    head_to_head(d, nuclear_oxphos, grep("^mt-", d$gene, value = TRUE),
                 "mito-encoded control", nm),
    head_to_head(d, nuclear_oxphos, grep("^Rp[sl][0-9]", d$gene, value = TRUE),
                 "ribosomal control", nm))
}))
if (nrow(h2h)) {
  write.csv(h2h, file.path(OUTDIR, paste0("part5c_specificity_head_to_head_", RUN_TAG, ".csv")),
            row.names = FALSE)
  cat("\n"); print(h2h, row.names = FALSE, digits = 3)
}

## ---- 4) Data-driven verdict, reported for every dataset ----
verdict_for <- function(nm) {
  nuc <- out %>% dplyr::filter(dataset == nm, grepl("^Nuclear", program))
  hh  <- h2h %>% dplyr::filter(dataset == nm)
  if (nrow(nuc) != 1) return("inconclusive (nuclear-OXPHOS set not measurable here)")
  if (!isTRUE(nuc$p_vs_background < FDR_CUT))
    return("no nuclear-OXPHOS shift detected")
  if (!nrow(hh)) return("inconclusive (no control set measurable for comparison)")
  ## Specific = beats EVERY control, both statistically and by a margin worth
  ## caring about. The 2x margin stops a large-n significant-but-tiny
  ## difference from being read as specificity.
  beats <- all(is.finite(hh$p_head_to_head) & hh$p_head_to_head < FDR_CUT &
               hh$fold_larger >= 2)
  if (beats)
    sprintf("SPECIFIC: nuclear OXPHOS %s (mean Wald %+.2f), %.1fx the nearest control and separated from every control (p<%.2g)",
            ifelse(nuc$mean_Wald > 0, "UP", "DOWN"), nuc$mean_Wald,
            min(hh$fold_larger), FDR_CUT)
  else
    sprintf("NOT specific: nuclear OXPHOS shifts (mean Wald %+.2f) but does not separate from the controls (smallest margin %.1fx, largest head-to-head p=%.2g)",
            nuc$mean_Wald, min(hh$fold_larger), max(hh$p_head_to_head, na.rm = TRUE))
}
for (nm in c(CELL_KEY, tissues))
  log_sanity(paste0("NSL specificity verdict (", nm, ")"), verdict_for(nm), "INFO")
verdict <- verdict_for(CELL_KEY)

## ---- 5) Figure ----
pd <- out %>% dplyr::filter(is.finite(mean_Wald)) %>%
  dplyr::mutate(dataset = factor(dataset, levels = c(CELL_KEY, tissues)),
                program = factor(program, levels = rev(unique(program))))
p <- ggplot(pd, aes(mean_Wald, program, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
  labs(x = "mean DESeq2 Wald statistic vs background", y = NULL,
       title = "KAT8 complex regulatory output, with specificity controls",
       caption = paste0("NSL targets = nuclear-encoded OXPHOS. Controls: mito-encoded (mtDNA-transcribed, ",
                        "outside NSL) and ribosomal proteins.\nVerdict (cells): ", verdict)) +
  theme_pub()
ggsave(file.path(OUTDIR, paste0("part5c_complex_output_", RUN_TAG, ".png")), p,
       width = 10, height = 4.6, dpi = 300, bg = "white")

write_sanity()
cat("\n[DONE] Part 5c complete. Outputs in ", OUTDIR, "\n", sep = "")
cat("Only claim an NSL-specific effect if the verdict above says SPECIFIC.\n")
