#!/usr/bin/env Rscript
## =========================================================
## TASK 7 - Harmonise DE_cells column names
## =========================================================
## DE_cells carries DESeq2 and limma-style names for the same quantities:
## pvalue / P.Value, and log2FoldChange / logFC. Worse, cells use adj.P.Val
## where both tissue tables use padj, so a script written across all three
## silently returns nothing for the cells instead of erroring.
##
## The duplication is verified BEFORE anything is dropped. If the columns are
## not identical they are not duplicates, and the script stops rather than
## deleting a column that carries different numbers.
##
## The NaN values in adj.P.Val are DESeq2 independent filtering -- genes the
## filter set aside to gain power. They are kept as NA. Imputing them, or
## dropping those rows, would silently change what the table means.
##
## Writes a NEW timestamped file. The original is not modified.
##
## USAGE
##   setwd("path/to/KAT8KD_RNAseq")
##   source("followups/task7_harmonise_de_cells.R")
## =========================================================

RUN    <- "savepoints/RUN_20260727_172121"
OUTDIR <- "followups/output"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

f <- list.files(RUN, pattern = "^DE_cells_.*_[0-9]{8}_[0-9]{6}\\.csv$",
                recursive = TRUE, full.names = TRUE)
if (!length(f)) stop("STOP: DE_cells table not found under ", RUN, call. = FALSE)
src <- f[order(basename(f))][length(f)]
d <- read.csv(src, row.names = 1, stringsAsFactors = FALSE)
cat("[INFO] source: ", basename(src), "  (", nrow(d), " rows)\n", sep = "")

need <- c("pvalue", "P.Value", "log2FoldChange", "logFC", "adj.P.Val")
miss <- setdiff(need, names(d))
if (length(miss))
  stop("STOP: expected columns absent: ", paste(miss, collapse = ", "),
       "\n      This file may already have been harmonised.", call. = FALSE)

## ---- verify duplication before removing anything ----
same <- function(a, b) {
  A <- d[[a]]; B <- d[[b]]
  eq <- (A == B) | (is.na(A) & is.na(B))
  n_bad <- sum(!eq, na.rm = TRUE) + sum(is.na(eq))
  cat(sprintf("[%s] %s == %s   mismatches: %d\n",
              ifelse(n_bad == 0, "OK  ", "FAIL"), a, b, n_bad))
  n_bad == 0
}
ok_p <- same("pvalue", "P.Value")
ok_l <- same("log2FoldChange", "logFC")
if (!(ok_p && ok_l))
  stop("STOP: columns are NOT identical, so they are not duplicates. ",
       "Nothing dropped - investigate before proceeding.", call. = FALSE)

n_na <- sum(is.na(d$adj.P.Val))
cat("[INFO] adj.P.Val NA (independent filtering, retained as NA): ", n_na, "\n", sep = "")

out <- d
out$P.Value <- NULL
out$logFC   <- NULL
names(out)[names(out) == "adj.P.Val"] <- "padj"

stopifnot(sum(is.na(out$padj)) == n_na, nrow(out) == nrow(d))

dest <- file.path(OUTDIR, paste0("DE_cells_KAT8KD_vs_CTL_harmonised_", STAMP, ".csv"))
write.csv(out, dest, row.names = TRUE)
cat("[OK] written: ", dest, "\n", sep = "")
cat("[OK] original untouched: ", src, "\n", sep = "")

md <- c(
  "# Task 7 - Harmonise DE_cells column names", "",
  paste0("Generated ", Sys.time()), "",
  paste0("Source: `", basename(src), "` (", nrow(d), " rows) - **not modified**"),
  paste0("Output: `", basename(dest), "`"), "",
  "## Verification (run before any column was dropped)", "",
  "| Check | Result |", "|---|---|",
  paste0("| `pvalue` == `P.Value` | identical, 0 mismatches |"),
  paste0("| `log2FoldChange` == `logFC` | identical, 0 mismatches |"),
  paste0("| `adj.P.Val` NA count | ", n_na, ", retained as NA |"), "",
  "## Changes", "",
  "- dropped `P.Value` (byte-identical to `pvalue`)",
  "- dropped `logFC` (byte-identical to `log2FoldChange`)",
  "- renamed `adj.P.Val` -> `padj`, matching both tissue tables", "",
  paste0("The ", n_na, " NA values in `padj` are DESeq2 independent filtering and are ",
         "retained as NA. They were not imputed, dropped, or filled."), "",
  "Columns now: `" , paste(names(out), collapse = "`, `"), "`", "",
  "## sessionInfo", "```", capture.output(sessionInfo()), "```")
writeLines(md, file.path(OUTDIR, paste0("task7_report_", STAMP, ".md")))
cat("[DONE] -> ", OUTDIR, "/task7_report_", STAMP, ".md\n", sep = "")
