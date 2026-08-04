#!/usr/bin/env Rscript
## =========================================================
## DIAGNOSE MISSING / SUPPRESSED GENES IN counts.txt
## =========================================================
## Part 4c found that Aldoa -- the canonical glycolytic aldolase -- is absent
## from all three DE tables, while Aldob and Aldoc are present in the tissue
## data and the cells table contains no aldolase at all. That is not a symbol
## mismatch (Hprt appears as "Hprt1", so the pipeline does carry variant
## symbols) and it is not an Ensembl-ID fallback (no ENSMUS* row corresponds
## to it). So the gene is lost UPSTREAM of the DE tables.
##
## There are exactly two possibilities and they have different fixes:
##
##   A. Aldoa is IN counts.txt but was removed by the prefilter
##      (rowSums(counts >= 10) >= 4). For a gene of Aldoa's expected
##      abundance that would mean its counts are near zero, which points at
##      QUANTIFICATION -- most likely multimapping reads being discarded,
##      since Aldoa has a large processed-pseudogene family in mouse.
##
##   B. Aldoa is NOT in counts.txt at all -> the ANNOTATION (GTF) used for
##      counting either lacks it or names it something else.
##
## This script decides which, and then checks whether the same thing has
## happened to other genes, because the answer changes how much of the
## pipeline is affected. Two observations from the DE tables motivate that
## wider check:
##   - Actg1 has baseMean 11 in iWAT while its paralogue Actb has 4069.
##   - Pgk1 has baseMean 27 in iWAT, which is very low for a glycolytic
##     housekeeping gene, and Pgk1 is one of the genes carrying the
##     iWAT-vs-gWAT discordance result. If its counts are unreliable, that
##     specific claim needs re-checking.
## Both are multimapping-prone. Neither is proof on its own -- hence this.
##
## READ-ONLY. Writes a report to followups/output/; changes nothing.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")     # the folder holding counts.txt
##   source("followups/diagnose_missing_genes.R")
## =========================================================

OUTDIR <- file.path("followups", "output")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

## Must match part1_main_analysis.R exactly, or "it was filtered" is not a
## statement about this pipeline.
PREFILTER_MIN_COUNT   <- 10
PREFILTER_MIN_SAMPLES <- 4

## 1) Load counts -------------------------------------------
if (!file.exists("counts.txt"))
  stop("counts.txt not found in ", getwd(),
       "\n  Run this from the project root, the same directory part1 runs from.",
       call. = FALSE)

counts_raw <- read.delim("counts.txt", header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
cat("counts.txt: ", nrow(counts_raw), " rows x ", ncol(counts_raw), " cols\n", sep = "")
cat("columns: ", paste(head(colnames(counts_raw), 8), collapse = ", "),
    if (ncol(counts_raw) > 8) " ..." else "", "\n", sep = "")

id_col   <- if ("gene_id" %in% colnames(counts_raw)) "gene_id" else colnames(counts_raw)[1]
name_col <- if ("gene" %in% colnames(counts_raw)) "gene" else
            if ("gene_name" %in% colnames(counts_raw)) "gene_name" else NA_character_
cat("using id column: ", id_col, " | symbol column: ",
    ifelse(is.na(name_col), "(none)", name_col), "\n", sep = "")

sample_cols <- setdiff(colnames(counts_raw), c(id_col, name_col))
sample_cols <- sample_cols[vapply(counts_raw[sample_cols], is.numeric, logical(1))]
cat("numeric sample columns: ", length(sample_cols), "\n", sep = "")
if (length(sample_cols) == 0) stop("No numeric sample columns found in counts.txt", call. = FALSE)

cm  <- as.matrix(counts_raw[, sample_cols, drop = FALSE])
sym <- if (is.na(name_col)) counts_raw[[id_col]] else as.character(counts_raw[[name_col]])
ids <- as.character(counts_raw[[id_col]])

## The prefilter part1 applies, recomputed here per gene.
n_pass  <- rowSums(cm >= PREFILTER_MIN_COUNT, na.rm = TRUE)
kept    <- n_pass >= PREFILTER_MIN_SAMPLES
tot     <- rowSums(cm, na.rm = TRUE)
cat("\nprefilter (counts >= ", PREFILTER_MIN_COUNT, " in >= ", PREFILTER_MIN_SAMPLES,
    " samples): ", sum(kept), " / ", length(kept), " genes kept\n", sep = "")

## 2) The Aldoa question ------------------------------------
report <- list()
say <- function(...) { line <- paste0(...); cat(line, "\n"); report[[length(report) + 1]] <<- line; invisible(NULL) }

say("")
say("=========================================================")
say(" ALDOA")
say("=========================================================")

## Search symbol AND id, and allow variant spellings, so "absent" is a
## conclusion rather than a failure to look properly.
pat <- "^Aldo|^ENSMUSG00000030695"
hits <- which(grepl(pat, sym, ignore.case = TRUE) | grepl(pat, ids, ignore.case = TRUE))

if (length(hits) == 0) {
  say("VERDICT: no Aldo* feature of any kind in counts.txt.")
  say("  -> Case B: the annotation used for counting does not contain it.")
  say("  -> Check the GTF/annotation the quantifier was given.")
} else {
  say(sprintf("%d Aldo*/Aldoa-ID feature(s) in counts.txt:", length(hits)))
  ald <- data.frame(id = ids[hits], symbol = sym[hits],
                    total_count = tot[hits], n_samples_ge10 = n_pass[hits],
                    passed_prefilter = kept[hits], stringsAsFactors = FALSE)
  ald <- ald[order(-ald$total_count), ]
  print(ald, row.names = FALSE)
  for (i in seq_len(nrow(ald)))
    report[[length(report) + 1]] <- paste(capture.output(print(ald[i, ], row.names = FALSE)), collapse = " ")

  exact <- which(tolower(sym) == "aldoa" | grepl("^ENSMUSG00000030695", ids))
  if (length(exact) == 0) {
    say("")
    say("VERDICT: Aldo* features exist, but NONE is Aldoa itself.")
    say("  -> Case B: annotation gap specific to Aldoa.")
  } else {
    say("")
    say("Aldoa IS present in counts.txt. Per-sample counts:")
    pr <- cm[exact, , drop = FALSE]
    rownames(pr) <- sym[exact]
    print(pr)
    report[[length(report) + 1]] <- paste(capture.output(print(pr)), collapse = "\n")
    if (any(kept[exact])) {
      say("")
      say("VERDICT: Aldoa PASSES the prefilter, so it should be in the DE tables.")
      say("  -> The loss is between prefiltering and the DE CSVs: check the")
      say("     symbol-mapping / deduplication step in part1, not the counts.")
    } else {
      say("")
      say(sprintf("VERDICT: Aldoa is present but FAILS the prefilter (max %d samples >= %d; total %d reads).",
                  max(n_pass[exact]), PREFILTER_MIN_COUNT, max(tot[exact])))
      say("  -> Case A: near-zero counts for a gene that should be abundant.")
      say("     Consistent with multimapping reads being discarded (Aldoa has a")
      say("     large processed-pseudogene family). Check the quantifier's")
      say("     multimapping setting; a method that assigns multireads by EM")
      say("     (Salmon/RSEM/STAR --quantMode with EM) would recover it.")
    }
  }
}

## 3) Is Aldoa alone? ---------------------------------------
## A single missing gene is a curiosity. A class of missing genes is a
## pipeline problem, and the two need different responses.
say("")
say("=========================================================")
say(" IS IT JUST ALDOA? canonical abundant genes")
say("=========================================================")

## Genes that should be unambiguously well expressed in adipocytes/adipose.
## Split by whether they have large paralogue/processed-pseudogene families,
## because that is the mechanism under suspicion.
CANON_MULTIMAP_PRONE <- c("Aldoa", "Gapdh", "Pgk1", "Actg1", "Eef1a1", "Eef2",
                           "Hspa8", "Ptma", "Npm1", "Hmgb1", "Rpl13a", "Ldha",
                           "Pkm", "Eno1", "Tpi1")
CANON_CONTROL        <- c("Actb", "Vim", "Fn1", "Sdha", "Calr", "Canx", "Ywhaz",
                           "Hnrnpk", "Rpl4", "Atp5f1b", "Ppia", "Tbp")

audit <- function(genes, label) {
  rows <- lapply(genes, function(g) {
    i <- which(tolower(sym) == tolower(g))
    if (length(i) == 0)
      return(data.frame(group = label, gene = g, in_counts = FALSE, total_count = NA_real_,
                        n_samples_ge10 = NA_integer_, passed_prefilter = NA,
                        stringsAsFactors = FALSE))
    i <- i[which.max(tot[i])]
    data.frame(group = label, gene = g, in_counts = TRUE, total_count = tot[i],
               n_samples_ge10 = n_pass[i], passed_prefilter = kept[i],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
aud <- rbind(audit(CANON_MULTIMAP_PRONE, "multimap-prone"),
             audit(CANON_CONTROL, "control"))
print(aud, row.names = FALSE)
report[[length(report) + 1]] <- paste(capture.output(print(aud, row.names = FALSE)), collapse = "\n")

lost <- aud[!is.na(aud$in_counts) & (!aud$in_counts | !aud$passed_prefilter), ]
say("")
if (nrow(lost) == 0) {
  say("All canonical genes checked are present and pass the prefilter.")
} else {
  say(sprintf("%d canonical gene(s) absent or filtered out: %s",
              nrow(lost), paste(lost$gene, collapse = ", ")))
  by_grp <- table(lost$group)
  say(paste0("  by group: ", paste(names(by_grp), by_grp, sep = "=", collapse = ", ")))
  say("  If these are concentrated in 'multimap-prone', the quantification step")
  say("  is discarding multireads and the effect is systematic, not gene-specific.")
}

## 4) Genes that look suppressed rather than absent ---------
## The damaging case is not a gene that vanishes -- that is visible. It is a
## gene that survives with a fraction of its true counts, because it then
## produces a real-looking but under-powered DE result. Pgk1 is the reason
## this matters here: it carries part of the iWAT-vs-gWAT discordance claim.
say("")
say("=========================================================")
say(" SUPPRESSED-BUT-PRESENT check (paralogue pairs)")
say("=========================================================")
say("Large ratios below are expected for genuinely different genes; they are")
say("listed so an implausible one (e.g. Actb >> Actg1) is easy to spot.")
pairs <- list(c("Actb", "Actg1"), c("Aldoa", "Aldoc"), c("Ldha", "Ldhb"),
              c("Pgk1", "Gapdh"), c("Eno1", "Eno3"))
prow <- lapply(pairs, function(p) {
  g <- function(x) { i <- which(tolower(sym) == tolower(x)); if (!length(i)) NA_real_ else max(tot[i]) }
  a <- g(p[1]); b <- g(p[2])
  data.frame(gene_a = p[1], count_a = a, gene_b = p[2], count_b = b,
             ratio_a_over_b = if (is.na(a) || is.na(b) || b == 0) NA_real_ else round(a / b, 2),
             stringsAsFactors = FALSE)
})
pr <- do.call(rbind, prow)
print(pr, row.names = FALSE)
report[[length(report) + 1]] <- paste(capture.output(print(pr, row.names = FALSE)), collapse = "\n")

## 5) Write the report --------------------------------------
f <- file.path(OUTDIR, paste0("missing_gene_diagnosis_", STAMP, ".txt"))
writeLines(unlist(report), f)
write.csv(aud, file.path(OUTDIR, paste0("canonical_gene_audit_", STAMP, ".csv")),
          row.names = FALSE)
cat("\n[OK] Wrote ", f, "\n", sep = "")
cat("[OK] Wrote ", file.path(OUTDIR, paste0("canonical_gene_audit_", STAMP, ".csv")), "\n", sep = "")
cat("\nWhat to do with the verdict:\n")
cat(" - Case A (present, near-zero counts): a quantification setting, and it\n")
cat("   affects every multimapping-prone gene, not just Aldoa.\n")
cat(" - Case B (absent from counts.txt): an annotation/GTF problem.\n")
cat(" - Passes prefilter but missing from the DE CSVs: the bug is in part1's\n")
cat("   symbol handling, and that is the one worth fixing immediately.\n")
