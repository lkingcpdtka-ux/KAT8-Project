#!/usr/bin/env Rscript
## =========================================================
## TASK 3 - Sensitivity check: gWAT down-regulated ORA
## =========================================================
## The gWAT DOWN direction produced no significant ORA enrichment at the
## primary threshold, and that absence is being reported as a finding. The
## primary threshold discards 81% of significant gWAT genes (6,195 reach
## padj < 0.05; 1,159 survive |log2FC| > 1), so the absence has to be shown
## not to be an artefact of the gate.
##
## Relaxed DOWN set: padj < 0.05 AND log2FoldChange < -0.5   (n = 1,399)
##
## PARAMETERS ARE INHERITED, NOT RETYPED. This script sources parameters.R and
## uses ora_params, so the cutoffs, set-size bounds and simplify threshold are
## by construction identical to the primary run. The universe is built the same
## way part2_ora.R builds it: every gene tested by DESeq2, mapped SYMBOL ->
## ENTREZID through org.Mm.eg.db. Using all annotated genes instead would
## inflate every term, which is the mistake this pipeline already avoids.
##
## SENSITIVITY CHECK ONLY. Nothing here substitutes into a primary analysis,
## figure or DEG count.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")
##   source("followups/task3_relaxed_gwat_ora.R")
## =========================================================

RUN    <- "savepoints/RUN_20260727_172121"
OUTDIR <- "followups/output"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

suppressPackageStartupMessages({
  library(dplyr); library(clusterProfiler); library(org.Mm.eg.db); library(AnnotationDbi)
})
source("parameters.R")
set.seed(if (exists("MASTER_SEED")) MASTER_SEED else 12345)

RELAXED_LFC <- -0.5
FDR_REPORT  <- 0.05

f <- list.files(RUN, pattern = "^DE_tissue_gWAT_.*_[0-9]{8}_[0-9]{6}\\.csv$",
                recursive = TRUE, full.names = TRUE)
if (!length(f)) stop("STOP: gWAT DE table not found under ", RUN, call. = FALSE)
de <- read.csv(f[order(basename(f))][length(f)], row.names = 1, stringsAsFactors = FALSE)
gn <- if ("gene_name" %in% names(de)) de$gene_name else rownames(de)

cat("[INFO] DE table: ", basename(f[length(f)]), "  (", nrow(de), " genes tested)\n", sep = "")
cat("[INFO] padj<0.05 any direction : ", sum(de$padj < 0.05, na.rm = TRUE), "\n", sep = "")
cat("[INFO] primary DOWN (lfc < -1) : ", sum(de$padj < 0.05 & de$log2FoldChange < -1, na.rm = TRUE), "\n", sep = "")

sel <- !is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange < RELAXED_LFC
cat("[INFO] RELAXED DOWN (lfc < ", RELAXED_LFC, ") : ", sum(sel), "\n", sep = "")

## ---- universe: every gene tested, exactly as part2_ora.R does it ----
uni <- suppressWarnings(bitr(gn, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
universe_entrez <- unique(uni$ENTREZID[!is.na(uni$ENTREZID)])
gsel <- suppressWarnings(bitr(gn[sel], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db))
entrez <- unique(gsel$ENTREZID[!is.na(gsel$ENTREZID)])
cat("[INFO] universe Entrez: ", length(universe_entrez),
    "  | relaxed-set Entrez: ", length(entrez), "\n", sep = "")
cat("[INFO] universe = all genes tested by DESeq2 (matches the primary run)\n")

ego <- enrichGO(gene = entrez, OrgDb = org.Mm.eg.db, ont = "BP",
                pvalueCutoff = ora_params$pvalue_cutoff,
                qvalueCutoff = ora_params$qvalue_cutoff,
                readable = TRUE,
                minGSSize = ora_params$min_gs_size,
                maxGSSize = ora_params$max_gs_size,
                universe = universe_entrez)
ego_s <- if (!is.null(ego) && nrow(ego@result) > 0)
  simplify(ego, cutoff = ora_params$simplify_cutoff, by = "p.adjust", select_fun = min) else NULL

ekg <- enrichKEGG(gene = entrez, organism = "mmu",
                  pvalueCutoff = ora_params$pvalue_cutoff,
                  qvalueCutoff = ora_params$qvalue_cutoff,
                  minGSSize = ora_params$min_gs_size,
                  maxGSSize = ora_params$max_gs_size,
                  universe = universe_entrez)

grab <- function(x) if (is.null(x) || !nrow(as.data.frame(x))) data.frame() else as.data.frame(x)
gores <- grab(ego_s); kgres <- grab(ekg)
nsig <- function(d) if (!nrow(d)) 0L else sum(d$p.adjust < FDR_REPORT, na.rm = TRUE)

if (nrow(gores)) write.csv(gores, file.path(OUTDIR, paste0("task3_RELAXED_gobp_gWAT_Down_", STAMP, ".csv")), row.names = FALSE)
if (nrow(kgres)) write.csv(kgres, file.path(OUTDIR, paste0("task3_RELAXED_kegg_gWAT_Down_", STAMP, ".csv")), row.names = FALSE)

top15 <- function(d) {
  if (!nrow(d)) return("_no terms returned_")
  d <- d[order(d$p.adjust), ][seq_len(min(15, nrow(d))), ]
  c("| Term | GeneRatio | p.adjust |", "|---|---|---|",
    sprintf("| %s | %s | %.3g |", d$Description, d$GeneRatio, d$p.adjust))
}

md <- c(
  "# Task 3 - Sensitivity check: gWAT down-regulated ORA", "",
  paste0("Generated ", Sys.time()), "",
  "## Setup", "",
  paste0("- Relaxed set: padj < 0.05 AND log2FC < ", RELAXED_LFC, " -> **", sum(sel), " genes**"),
  paste0("- Primary set for comparison: padj < 0.05 AND log2FC < -1 -> ",
         sum(de$padj < 0.05 & de$log2FoldChange < -1, na.rm = TRUE), " genes"),
  paste0("- Universe: **all ", nrow(de), " genes tested by DESeq2**, ",
         length(universe_entrez), " mapped to Entrez - identical construction to the primary run"),
  paste0("- Parameters inherited from `ora_params` in parameters.R: p ", ora_params$pvalue_cutoff,
         ", q ", ora_params$qvalue_cutoff, ", set size ", ora_params$min_gs_size, "-",
         ora_params$max_gs_size, ", simplify ", ora_params$simplify_cutoff),
  paste0("- org.Mm.eg.db ", as.character(packageVersion("org.Mm.eg.db")),
         " | clusterProfiler ", as.character(packageVersion("clusterProfiler"))), "",
  "## Result", "",
  paste0("- GO:BP terms returned: ", nrow(gores), "; significant at FDR < ", FDR_REPORT, ": **", nsig(gores), "**"),
  paste0("- KEGG terms returned: ", nrow(kgres), "; significant at FDR < ", FDR_REPORT, ": **", nsig(kgres), "**"), "",
  "### Top GO:BP", top15(gores), "", "### Top KEGG", top15(kgres), "",
  "## Verdict", "",
  if (nsig(gores) == 0 && nsig(kgres) == 0)
    paste0("Relaxing the gate from |log2FC| > 1 to < ", RELAXED_LFC, " (", sum(sel),
           " genes, up from ", sum(de$padj < 0.05 & de$log2FoldChange < -1, na.rm = TRUE),
           ") still returns no significant enrichment. **The absence is not a threshold artefact** ",
           "and can be reported as a real finding.")
  else
    paste0("The relaxed set DOES return significant enrichment (", nsig(gores), " GO:BP, ",
           nsig(kgres), " KEGG). **The reported absence is at least partly a threshold artefact** ",
           "and the claim must be softened: state that no enrichment survives the primary ",
           "|log2FC| > 1 gate, not that gWAT down-regulation is unstructured. ",
           "Do NOT substitute these terms into the primary results."), "",
  "_Sensitivity check only. No primary analysis, figure or DEG count is altered._", "",
  "## sessionInfo", "```", capture.output(sessionInfo()), "```")
writeLines(md, file.path(OUTDIR, paste0("task3_report_", STAMP, ".md")))
cat("\n", paste(md[c(6:10, 14:17)], collapse = "\n"), "\n", sep = "")
cat("\n[DONE] -> ", OUTDIR, "/task3_report_", STAMP, ".md\n", sep = "")
