#!/usr/bin/env Rscript
## =========================================================
## PANEL: plasma membrane / ECM  -- same genes, this pipeline
## =========================================================
## Takes the nine-gene panel from the standalone 8-sample script and asks it
## against THIS study: both tissue depots and the 3T3-L1 cells, using the
## pipeline's own VST matrix, colour scheme and z-score convention.
##
## The point is not to redraw someone else's figure. It is to answer whether
## the panel reproduces here -- so the DE table across all three contrasts is
## the real output, and the heatmap is the illustration.
##
## Differences from the standalone script, all deliberate:
##   - VST (DESeq2 variance-stabilised), not log2(CPM+1). The pipeline uses VST
##     everywhere for visualisation; mixing the two makes figures incomparable.
##   - Column names resolved for BOTH DESeq2 (padj) and limma (adj.P.Val)
##     spellings, so it survives the cells-table rename either way.
##   - Gene symbols resolved through an alias table. Ccn2/Ctgf is the live
##     example: a renamed gene silently drops out of a hand-curated panel.
##   - A missing panel gene is FATAL, not a warning. A heatmap quietly missing
##     a row is worse than one that refuses to draw.
##   - Effect size is printed next to every row, because a per-gene z-score
##     stretches a 5% change and a 4-fold change to the same colour.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")
##   source("followups/panel_plasma_membrane_ecm.R")
## =========================================================

RUN    <- "savepoints/RUN_20260727_172121"
OUTDIR <- "followups/output"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

suppressPackageStartupMessages({
  library(dplyr); library(ComplexHeatmap); library(circlize); library(grid)
})
source("parameters.R")

PANEL <- c("Car12","C1qtnf6","Tnfaip6","Ccn2","Plaur","Gpr137b","Pcolce2","Fxyd5","Itgb3")

## Renamed symbols. Ccn2 was Ctgf; annotations built before the CCN
## renaming still carry the old name, and a panel gene that silently fails to
## match is indistinguishable from one that is genuinely absent.
ALIASES <- list(Ccn2 = c("Ccn2","Ctgf"), C1qtnf6 = c("C1qtnf6","Ctrp6"),
                Plaur = c("Plaur","Upar"), Fxyd5 = c("Fxyd5","Dysadherin"))
resolve <- function(g, pool) {
  for (nm in c(ALIASES[[g]], g)) if (!is.null(nm) && nm %in% pool) return(nm)
  NA_character_
}

find1 <- function(pat, dir = RUN) {
  f <- list.files(dir, pattern = pat, recursive = TRUE, full.names = TRUE)
  if (!length(f)) return(NA_character_)
  f[order(basename(f))][length(f)]
}
padj_of <- function(d) if ("padj" %in% names(d)) d$padj else d$adj.P.Val
gene_of <- function(d) if ("gene_name" %in% names(d)) d$gene_name else d[[1]]

## ---- 1) DE across all three contrasts -------------------------------------
contrasts <- c(iWAT  = "^DE_tissue_iWAT_.*_[0-9]{8}_[0-9]{6}\\.csv$",
               gWAT  = "^DE_tissue_gWAT_.*_[0-9]{8}_[0-9]{6}\\.csv$",
               cells = "^DE_cells_.*_[0-9]{8}_[0-9]{6}\\.csv$")
de <- lapply(contrasts, function(p) { f <- find1(p); if (is.na(f)) NULL else
  read.csv(f, stringsAsFactors = FALSE) })
if (all(vapply(de, is.null, logical(1)))) stop("STOP: no DE tables found under ", RUN, call. = FALSE)

rows <- list()
for (cn in names(de)) {
  d <- de[[cn]]; if (is.null(d)) next
  sym <- vapply(PANEL, resolve, character(1), pool = gene_of(d))
  i   <- match(sym, gene_of(d))
  rows[[cn]] <- data.frame(
    contrast = cn, panel_gene = PANEL, matched_as = sym,
    log2FC = d$log2FoldChange[i], lfcSE = d$lfcSE[i], padj = padj_of(d)[i],
    stringsAsFactors = FALSE)
}
de_tab <- dplyr::bind_rows(rows) %>%
  dplyr::mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "*", ""))
write.csv(de_tab, file.path(OUTDIR, paste0("panel_pm_ecm_DE_", STAMP, ".csv")), row.names = FALSE)

cat("\n=== panel across all three contrasts ===\n")
print(de_tab, row.names = FALSE, digits = 3)

lost <- unique(de_tab$panel_gene[is.na(de_tab$matched_as)])
if (length(lost))
  cat("\n[NOTE] not measured in at least one contrast: ", paste(lost, collapse = ", "),
      "\n       (absent from the filtered matrix = below the expression floor,\n",
      "        which is itself a result -- not a failure to look)\n", sep = "")

## ---- 2) Heatmap on the pipeline's own VST ---------------------------------
qf <- find1("^QC_data_.*\\.rds$")
if (is.na(qf)) {
  cat("\n[SKIP] no QC_data RDS found; DE table written, heatmap skipped\n")
} else {
  qc  <- readRDS(qf)
  vst <- qc$vst_mat_qc; si <- qc$sample_info
  si  <- si[match(colnames(vst), si$Sample), , drop = FALSE]

  sym <- vapply(PANEL, resolve, character(1), pool = rownames(vst))
  if (anyNA(sym))
    stop("STOP: panel genes absent from the VST matrix: ",
         paste(PANEL[is.na(sym)], collapse = ", "),
         "\n      Checked aliases too. Fix the panel rather than dropping rows.",
         call. = FALSE)

  m <- vst[sym, , drop = FALSE]
  rownames(m) <- ifelse(sym == PANEL, PANEL, paste0(PANEL, " (", sym, ")"))

  ## z-score per gene, clipped -- the pipeline's convention, from parameters.R
  clip <- heatmap_params$lfc_clip
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  z <- pmax(pmin(z, clip), -clip)

  col_fun <- circlize::colorRamp2(
    seq(-clip, clip, length.out = length(heatmap_params$custom_gradient)),
    heatmap_params$custom_gradient)

  ## Effect size shown beside the map, because a per-gene z-score cannot
  ## distinguish a 5% change from a 4-fold one -- every row is stretched to
  ## the same range by construction.
  lab <- function(cn) {
    x <- de_tab[de_tab$contrast == cn, ]
    sprintf("%+.2f%s", x$log2FC[match(PANEL, x$panel_gene)],
            x$sig[match(PANEL, x$panel_gene)])
  }
  ra <- rowAnnotation(
    `iWAT log2FC` = anno_text(lab("iWAT"), gp = gpar(fontsize = 9)),
    `gWAT log2FC` = anno_text(lab("gWAT"), gp = gpar(fontsize = 9)),
    `cells log2FC` = anno_text(lab("cells"), gp = gpar(fontsize = 9)))

  ha <- HeatmapAnnotation(
    Genotype = si$Genotype, Sex = si$Sex,
    col = list(Genotype = c(CTL = "#0072B2", KAT8KD = "#D55E00"),
               Sex      = c(F = "#E7298A", M = "#7570B3")))

  ht <- Heatmap(
    z, name = "z-score", col = col_fun,
    top_annotation = ha, right_annotation = ra,
    column_split = si$Depot, cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(title = "z-score (VST)"))

  png(file.path(OUTDIR, paste0("panel_pm_ecm_heatmap_", STAMP, ".png")),
      width = 11, height = 4.6, units = "in", res = 320)
  draw(ht, column_title = "Plasma membrane / ECM panel  -- VST z-score, this study",
       column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()

  pdf(file.path(OUTDIR, paste0("panel_pm_ecm_heatmap_", STAMP, ".pdf")),
      width = 11, height = 4.6)
  draw(ht, column_title = "Plasma membrane / ECM panel  -- VST z-score, this study")
  dev.off()
  cat("\n[OK] heatmap -> ", OUTDIR, "/panel_pm_ecm_heatmap_", STAMP, ".png\n", sep = "")
}

cat("[DONE] DE table -> ", OUTDIR, "/panel_pm_ecm_DE_", STAMP, ".csv\n", sep = "")
