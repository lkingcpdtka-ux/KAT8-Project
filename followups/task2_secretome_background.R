#!/usr/bin/env Rscript
## =========================================================
## TASK 2 - Secretome reproduction vs background rate
## =========================================================
## Tests whether secreted/surface genes reproduce in vivo more often than
## other cell-significant genes do. Writes NEW timestamped files; touches
## nothing that already exists.
##
## Reproduction rule (identical to part5a): a gene reproduces in a depot if
##   tissue padj < 0.05  AND  sign(tissue log2FC) == sign(cell log2FC)
##
## Background: genes significant in cells (padj < 0.05) that are NOT in the
## secretome set. The comparison is only meaningful against this, because
## "45 of 78 reproduce" means nothing until you know what an arbitrary
## cell-significant gene does.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")
##   source("followups/task2_secretome_background.R")
## =========================================================

RUN      <- "savepoints/RUN_20260727_172121"
SECRETO  <- "handoff/part5a_secretome_annotated_20260728_000412.csv"
OUTDIR   <- "followups/output"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP    <- format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(12345)

find1 <- function(pat) {
  f <- list.files(RUN, pattern = pat, recursive = TRUE, full.names = TRUE)
  f <- f[grepl("_[0-9]{8}_[0-9]{6}\\.csv$", basename(f))]
  if (!length(f)) stop("STOP: no file matching ", pat, " under ", RUN, call. = FALSE)
  f[order(basename(f))][length(f)]
}
if (!file.exists(SECRETO)) stop("STOP: secretome file not found: ", SECRETO, call. = FALSE)

cells <- read.csv(find1("^DE_cells_"), stringsAsFactors = FALSE)
iwat  <- read.csv(find1("^DE_tissue_iWAT_"), stringsAsFactors = FALSE)
gwat  <- read.csv(find1("^DE_tissue_gWAT_"), stringsAsFactors = FALSE)
sec   <- read.csv(SECRETO, stringsAsFactors = FALSE)

## column names differ between cells (adj.P.Val) and tissue (padj)
padj_of <- function(d) if ("padj" %in% names(d)) d$padj else d$adj.P.Val
gene_of <- function(d) if ("gene_name" %in% names(d)) d$gene_name else d[[1]]

## secretome gene symbols: take whichever column overlaps the cell genes most
cand <- names(sec)[vapply(sec, is.character, logical(1))]
if (!length(cand)) stop("STOP: no character column in the secretome file", call. = FALSE)
ov   <- vapply(cand, function(c) length(intersect(sec[[c]], gene_of(cells))), integer(1))
sec_genes <- unique(sec[[cand[which.max(ov)]]])
cat("[INFO] secretome gene column: '", cand[which.max(ov)], "' -> ",
    length(sec_genes), " unique symbols\n", sep = "")

cg  <- gene_of(cells); cp <- padj_of(cells); cl <- cells$log2FoldChange
sig <- !is.na(cp) & cp < 0.05
cat("[INFO] genes significant in cells (padj < 0.05): ", sum(sig), "\n", sep = "")

## reproduction flags, evaluated on the cell-significant set
lookup <- function(tis, genes) {
  t <- tis[!duplicated(gene_of(tis)), ]
  i <- match(genes, gene_of(t))
  data.frame(padj = padj_of(t)[i], lfc = t$log2FoldChange[i])
}
gs <- cg[sig]; ls <- cl[sig]
ti <- lookup(iwat, gs); tg <- lookup(gwat, gs)
rep_i <- !is.na(ti$padj) & ti$padj < 0.05 & sign(ti$lfc) == sign(ls)
rep_g <- !is.na(tg$padj) & tg$padj < 0.05 & sign(tg$lfc) == sign(ls)
in_sec <- gs %in% sec_genes

tab <- data.frame(gene = gs, cell_log2FC = ls, in_secretome = in_sec,
                  repro_iWAT = rep_i, repro_gWAT = rep_g,
                  repro_any = rep_i | rep_g, repro_both = rep_i & rep_g,
                  stringsAsFactors = FALSE)
write.csv(tab, file.path(OUTDIR, paste0("task2_reproduction_per_gene_", STAMP, ".csv")),
          row.names = FALSE)

ft <- function(hit_s, n_s, hit_b, n_b) {
  m <- matrix(c(hit_s, n_s - hit_s, hit_b, n_b - hit_b), nrow = 2, byrow = TRUE)
  f <- fisher.test(m)
  list(or = unname(f$estimate), p = f$p.value,
       ps = hit_s / n_s, pb = hit_b / n_b, m = m)
}
n_s <- sum(in_sec); n_b <- sum(!in_sec)
res <- list(
  any  = ft(sum(tab$repro_any  &  in_sec), n_s, sum(tab$repro_any  & !in_sec), n_b),
  both = ft(sum(tab$repro_both &  in_sec), n_s, sum(tab$repro_both & !in_sec), n_b))

## ---- optional refinement: background matched on baseMean and |cell log2FC| ----
## Confirms the both-depot result is not simply an expression-level effect:
## highly expressed genes are easier to call significant in tissue too.
bm <- cells$baseMean[sig]
strata <- function(x, k = 5) cut(rank(x, ties.method = "first"), k, labels = FALSE)
key <- paste(strata(bm), strata(abs(ls)))
matched <- logical(length(gs))
for (k in unique(key[in_sec])) {
  need <- sum(key == k & in_sec)
  pool <- which(key == k & !in_sec)
  if (length(pool)) matched[sample(pool, min(need, length(pool)))] <- TRUE
}
n_m <- sum(matched)
res$both_matched <- ft(sum(tab$repro_both & in_sec), n_s, sum(tab$repro_both & matched), n_m)
res$any_matched  <- ft(sum(tab$repro_any  & in_sec), n_s, sum(tab$repro_any  & matched), n_m)

fmt <- function(lbl, r, n_s, n_b) sprintf(
  "| %s | %d/%d (%.0f%%) | %d/%d (%.0f%%) | %.2f | %.4g | %s |",
  lbl, round(r$ps * n_s), n_s, r$ps * 100, round(r$pb * n_b), n_b, r$pb * 100,
  r$or, r$p, ifelse(r$p < 0.05, "**significant**", "not significant"))

md <- c(
  "# Task 2 - Secretome reproduction vs background", "",
  paste0("Generated ", Sys.time(), " | seed 12345"), "",
  "Reproduction: tissue padj < 0.05 AND sign(tissue log2FC) == sign(cell log2FC).",
  paste0("Cell-significant genes: **", sum(sig), "**; secretome **", n_s,
         "**; background **", n_b, "**."), "",
  "| Comparison | Secretome | Background | OR | Fisher p | Verdict |",
  "|---|---|---|---|---|---|",
  fmt(">=1 depot", res$any, n_s, n_b),
  fmt("BOTH depots", res$both, n_s, n_b),
  fmt(">=1 depot (matched bg)", res$any_matched, n_s, n_m),
  fmt("BOTH depots (matched bg)", res$both_matched, n_s, n_m), "",
  "## Verdict", "",
  paste0("- `>=1 depot`: ", ifelse(res$any$p < 0.05, "SURVIVES", "**DOES NOT SURVIVE**"),
         " testing against background. ",
         ifelse(res$any$p < 0.05, "",
                "The manuscript's '45 of 78 reproduce' is at background rate and must not be quoted as evidence.")),
  paste0("- `BOTH depots`: ", ifelse(res$both$p < 0.05, "SURVIVES", "does not survive"),
         ". This is the defensible statistic."),
  paste0("- Expression-matched background: both-depot p = ", signif(res$both_matched$p, 3),
         " -> conclusion ", ifelse((res$both_matched$p < 0.05) == (res$both$p < 0.05),
                                   "UNCHANGED", "**CHANGED - investigate**"), "."), "",
  "## sessionInfo", "```", capture.output(sessionInfo()), "```")
writeLines(md, file.path(OUTDIR, paste0("task2_report_", STAMP, ".md")))
cat(paste(md[1:16], collapse = "\n"), "\n")
cat("\n[DONE] -> ", OUTDIR, "/task2_report_", STAMP, ".md\n", sep = "")
