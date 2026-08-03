#!/usr/bin/env Rscript
## =========================================================
## Plasma membrane / ECM panel -- cells only, pipeline-style heatmap
## =========================================================
## Same nine genes and same eight samples (JS_41-48) as the standalone
## 8-sample script, but drawn with this pipeline's own DESeq2 VST matrix and
## heatmap conventions, so it can be compared directly against the other
## version of the figure.
##
## The DE table (already verified: all 9 genes significant, all up,
## FDR 6.4e-9 to 0.031) is not recomputed here -- it is read from
## DE_cells_KAT8KD_vs_CTL, which came from the SAME DESeq2 fit used elsewhere
## in this pipeline. This script only rebuilds VST for the heatmap, because
## the cells step does not currently save its VST matrix to disk.
##
## Differences from the standalone script, all deliberate:
##   - VST, not log2(CPM+1) -- matches every other heatmap in this pipeline.
##   - log2FC and FDR from the pipeline's own DE table are printed beside each
##     row. A per-gene z-score stretches every row to the same amplitude, so
##     a 1.7-fold and a 3.5-fold change would otherwise look identical.
##   - A panel gene absent from the count matrix is FATAL, not a warning.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")
##   source("followups/panel_pm_ecm_cells_heatmap.R")
## =========================================================

OUTDIR <- "followups/output"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(ComplexHeatmap); library(circlize); library(grid)
})
source("parameters.R")

PANEL <- c("Car12","C1qtnf6","Tnfaip6","Ccn2","Plaur","Gpr137b","Pcolce2","Fxyd5","Itgb3")
CELL_SAMPLES <- paste0("JS_", sprintf("%02d", 41:48))

## ---- 1) DE numbers, from the pipeline's own cells DE table (not recomputed) ----
de_file <- list.files("savepoints/RUN_20260727_172121", pattern = "^DE_cells_.*\\.csv$",
                      recursive = TRUE, full.names = TRUE)
if (!length(de_file))
  de_file <- "downstream_analysis/data/DE_cells_KAT8KD_vs_CTL.csv"
de_file <- de_file[order(basename(de_file))][length(de_file)]
de <- read.csv(de_file, stringsAsFactors = FALSE)
padj_col <- if ("padj" %in% names(de)) "padj" else "adj.P.Val"

de_panel <- de %>%
  dplyr::filter(gene_name %in% PANEL) %>%
  dplyr::transmute(gene_name, log2FC = log2FoldChange, padj = .data[[padj_col]]) %>%
  { .[match(PANEL, .$gene_name), ] }
if (anyNA(de_panel$gene_name))
  stop("STOP: panel genes missing from DE table: ",
       paste(PANEL[is.na(de_panel$gene_name)], collapse = ", "), call. = FALSE)

cat("\n=== panel DE (from pipeline's cells DESeq2 fit, ", basename(de_file), ") ===\n", sep = "")
print(de_panel, row.names = FALSE, digits = 3)

## ---- 2) VST for the heatmap only. The cells DE step does not save a VST ----
## matrix to disk, so it is rebuilt here from the same counts.txt and the
## same 8 samples -- NOT re-run to get new statistics, only to get a
## per-sample matrix to colour.
counts_raw <- read.delim("counts.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
missing_s <- setdiff(CELL_SAMPLES, colnames(counts_raw))
if (length(missing_s)) stop("STOP: cell sample columns missing from counts.txt: ",
                            paste(missing_s, collapse = ", "), call. = FALSE)

cc <- counts_raw[, c("gene_id", "gene", CELL_SAMPLES)]
gn <- cc$gene; gn[is.na(gn) | gn == ""] <- cc$gene_id[is.na(gn) | gn == ""]
gn <- make.unique(gn)
rownames(cc) <- gn
mat <- as.matrix(cc[, CELL_SAMPLES])
storage.mode(mat) <- "integer"

missing_panel <- setdiff(PANEL, rownames(mat))
if (length(missing_panel))
  stop("STOP: panel genes absent from counts.txt: ", paste(missing_panel, collapse = ", "),
       call. = FALSE)

sample_annot <- data.frame(Sample = CELL_SAMPLES,
                           Genotype = factor(c(rep("CTL", 4), rep("KAT8KD", 4)),
                                             levels = c("CTL", "KAT8KD")),
                           row.names = CELL_SAMPLES)
keep <- rowSums(mat >= 10) >= 4          ## same prefilter rule as the rest of the pipeline
dds  <- DESeqDataSetFromMatrix(mat[keep, , drop = FALSE], colData = sample_annot, design = ~ Genotype)
vst_mat <- assay(vst(dds, blind = FALSE))

m <- vst_mat[PANEL, , drop = FALSE]

clip <- heatmap_params$lfc_clip
z <- t(scale(t(m)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, clip), -clip)

col_fun <- circlize::colorRamp2(
  seq(-clip, clip, length.out = length(heatmap_params$custom_gradient)),
  heatmap_params$custom_gradient)

lab <- sprintf("%+.2f%s", de_panel$log2FC,
              ifelse(!is.na(de_panel$padj) & de_panel$padj < 0.05, "*", ""))
ra <- rowAnnotation(`log2FC (padj)` = anno_text(lab, gp = gpar(fontsize = 9)))

ha <- HeatmapAnnotation(Genotype = sample_annot$Genotype,
                        col = list(Genotype = c(CTL = "#0072B2", KAT8KD = "#D55E00")))

ht <- Heatmap(z, name = "z-score", col = col_fun,
             top_annotation = ha, right_annotation = ra,
             column_split = sample_annot$Genotype, cluster_columns = FALSE, cluster_rows = FALSE,
             show_column_names = TRUE, column_names_gp = gpar(fontsize = 9),
             row_names_gp = gpar(fontsize = 10, fontface = "bold"),
             heatmap_legend_param = list(title = "z-score (VST)"))

png(file.path(OUTDIR, paste0("panel_pm_ecm_cells_heatmap_", STAMP, ".png")),
    width = 7.5, height = 4.2, units = "in", res = 320)
draw(ht, column_title = "Plasma membrane / ECM panel -- cells (JS_41-48), this pipeline's VST",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

pdf(file.path(OUTDIR, paste0("panel_pm_ecm_cells_heatmap_", STAMP, ".pdf")), width = 7.5, height = 4.2)
draw(ht, column_title = "Plasma membrane / ECM panel -- cells (JS_41-48), this pipeline's VST")
dev.off()

write.csv(de_panel, file.path(OUTDIR, paste0("panel_pm_ecm_cells_DE_", STAMP, ".csv")), row.names = FALSE)
cat("\n[OK] heatmap -> ", OUTDIR, "/panel_pm_ecm_cells_heatmap_", STAMP, ".png\n", sep = "")
cat("[OK] DE table -> ", OUTDIR, "/panel_pm_ecm_cells_DE_", STAMP, ".csv\n", sep = "")
