#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 1.5: DEPOT OVERLAP ANALYSIS
## =========================================================
## Compares iWAT vs gWAT DEGs to identify:
##   - Shared vs. depot-specific responses
##   - Concordance of effect directions
##   - logFC correlation between depots
## =========================================================

## 1) Load packages & parameters ----------------------------
required_pkgs <- c("dplyr", "ggplot2", "ggrepel", "RColorBrewer", "VennDiagram", "grid", "gridExtra")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
  library(gridExtra)
})

## Source central parameters for thresholds and colors
source(file.path(getwd(), "parameters.R"))

## 2) Find most recent Part 1 run ---------------------------
cat("\n=== FINDING PART 1 RESULTS ===\n")
savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run part1_main_analysis.R first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No Part 1 run directories found in savepoints/")
}

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using Part 1 results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

## 3) Parameters (from parameters.R) ------------------------
logFC_cut <- default_logFC_cut
fdr_cut   <- default_fdr_cut

## RdYlBu palette colors for consistent styling
col_iwat     <- "#4575B4"   ## blue (matches iWAT in QC plots)
col_gwat     <- "#D73027"   ## red  (matches gWAT in QC plots)
col_both     <- "#FC8D59"   ## orange (shared)
col_neither  <- "gray80"

## 4) Load DE tables ----------------------------------------
cat("\n=== LOADING DE TABLES ===\n")
tables_dir <- file.path(outdir, "tables")

de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("OVERALL", de_files)]  # Exclude any OVERALL files if present

if (length(de_files) != 2) {
  stop("Expected 2 DE tables (iWAT and gWAT), found: ", length(de_files))
}

## Load both tables
de_list <- lapply(de_files, function(f) {
  read.csv(f, row.names = 1, stringsAsFactors = FALSE)
})

## Extract contrast names
contrast_names <- sapply(de_files, function(f) {
  gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
})
names(de_list) <- contrast_names

## Identify which is iWAT and which is gWAT
iwat_idx <- grep("iWAT", contrast_names)
gwat_idx <- grep("gWAT", contrast_names)

if (length(iwat_idx) != 1 || length(gwat_idx) != 1) {
  stop("Could not identify iWAT and gWAT contrasts")
}

de_iwat <- de_list[[iwat_idx]]
de_gwat <- de_list[[gwat_idx]]

cat("[OK] Loaded iWAT contrast: ", contrast_names[iwat_idx], "\n", sep = "")
cat("[OK] Loaded gWAT contrast: ", contrast_names[gwat_idx], "\n", sep = "")

## 5) Identify DEG sets -------------------------------------
cat("\n=== IDENTIFYING DEG SETS ===\n")

## Use depot-specific thresholds from parameters.R
iwat_thresh <- get_tissue_thresholds("iWAT")
gwat_thresh <- get_tissue_thresholds("gWAT")

## iWAT DEGs
iwat_up <- rownames(de_iwat %>%
  dplyr::filter(!is.na(padj), padj < iwat_thresh$fdr_cut, logFC > iwat_thresh$logFC_cut))
iwat_down <- rownames(de_iwat %>%
  dplyr::filter(!is.na(padj), padj < iwat_thresh$fdr_cut, logFC < -iwat_thresh$logFC_cut))
iwat_all_deg <- c(iwat_up, iwat_down)

## gWAT DEGs
gwat_up <- rownames(de_gwat %>%
  dplyr::filter(!is.na(padj), padj < gwat_thresh$fdr_cut, logFC > gwat_thresh$logFC_cut))
gwat_down <- rownames(de_gwat %>%
  dplyr::filter(!is.na(padj), padj < gwat_thresh$fdr_cut, logFC < -gwat_thresh$logFC_cut))
gwat_all_deg <- c(gwat_up, gwat_down)

cat("[INFO] iWAT DEGs: ", length(iwat_all_deg), " (", length(iwat_up), " up, ", length(iwat_down), " down)\n", sep = "")
cat("[INFO] gWAT DEGs: ", length(gwat_all_deg), " (", length(gwat_up), " up, ", length(gwat_down), " down)\n", sep = "")

## 6) Calculate overlaps ------------------------------------
cat("\n=== CALCULATING OVERLAPS ===\n")

## Overall overlap
overlap_all <- intersect(iwat_all_deg, gwat_all_deg)
iwat_specific <- setdiff(iwat_all_deg, gwat_all_deg)
gwat_specific <- setdiff(gwat_all_deg, iwat_all_deg)

cat("[INFO] Shared DEGs (any direction): ", length(overlap_all), "\n", sep = "")
cat("[INFO] iWAT-specific DEGs: ", length(iwat_specific), "\n", sep = "")
cat("[INFO] gWAT-specific DEGs: ", length(gwat_specific), "\n", sep = "")

## Direction-specific overlaps
overlap_up_up <- intersect(iwat_up, gwat_up)
overlap_down_down <- intersect(iwat_down, gwat_down)
overlap_opposite <- union(
  intersect(iwat_up, gwat_down),
  intersect(iwat_down, gwat_up)
)

cat("\n--- Direction Concordance ---\n")
cat("[INFO] Both UP: ", length(overlap_up_up), " genes\n", sep = "")
cat("[INFO] Both DOWN: ", length(overlap_down_down), " genes\n", sep = "")
cat("[INFO] Opposite directions: ", length(overlap_opposite), " genes\n", sep = "")
concordance_pct <- if (length(overlap_all) > 0) {
  round(100 * (length(overlap_up_up) + length(overlap_down_down)) / length(overlap_all), 1)
} else {
  NA_real_
}
cat("[INFO] Concordance rate: ", concordance_pct, "%\n", sep = "")

## 7) Create overlap summary table --------------------------
overlap_summary <- data.frame(
  Category = c(
    "Total DEGs (iWAT)",
    "Total DEGs (gWAT)",
    "Shared DEGs (any direction)",
    "iWAT-specific DEGs",
    "gWAT-specific DEGs",
    "Shared: both UP",
    "Shared: both DOWN",
    "Shared: opposite directions",
    "Concordance rate (%)"
  ),
  Count = c(
    length(iwat_all_deg),
    length(gwat_all_deg),
    length(overlap_all),
    length(iwat_specific),
    length(gwat_specific),
    length(overlap_up_up),
    length(overlap_down_down),
    length(overlap_opposite),
    concordance_pct
  ),
  stringsAsFactors = FALSE
)

overlap_file <- paste0("Overlap_summary_", run_tag, ".csv")
write.csv(overlap_summary, file = file.path(outdir, "tables", overlap_file), row.names = FALSE)
cat("\n[OK] Saved overlap summary: ", overlap_file, "\n", sep = "")

## 8) Venn diagrams -----------------------------------------
cat("\n=== CREATING VENN DIAGRAMS ===\n")

## Overall Venn (all DEGs)
venn_all_file <- paste0("Venn_all_DEGs_", run_tag, ".png")
png(file.path(outdir, "plots", venn_all_file), width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_all <- venn.diagram(
  x = list(
    iWAT = iwat_all_deg,
    gWAT = gwat_all_deg
  ),
  category.names = c("iWAT", "gWAT"),
  filename = NULL,
  output = TRUE,
  fill = c(col_iwat, col_gwat),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "All DEGs (FDR<0.05, |logFC|>1)",
  main.cex = 1.5,
  main.pos = c(0.5, 1.05)
)
grid.draw(venn_all)
popViewport()
dev.off()
cat("[OK] Saved Venn diagram (all DEGs): ", venn_all_file, "\n", sep = "")

## UP-regulated Venn
venn_up_file <- paste0("Venn_UP_DEGs_", run_tag, ".png")
png(file.path(outdir, "plots", venn_up_file), width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_up <- venn.diagram(
  x = list(
    iWAT_UP = iwat_up,
    gWAT_UP = gwat_up
  ),
  category.names = c("iWAT UP", "gWAT UP"),
  filename = NULL,
  output = TRUE,
  fill = c(col_iwat, col_gwat),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "UP-regulated DEGs",
  main.cex = 1.5,
  main.pos = c(0.5, 1.05)
)
grid.draw(venn_up)
popViewport()
dev.off()
cat("[OK] Saved Venn diagram (UP DEGs): ", venn_up_file, "\n", sep = "")

## DOWN-regulated Venn
venn_down_file <- paste0("Venn_DOWN_DEGs_", run_tag, ".png")
png(file.path(outdir, "plots", venn_down_file), width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_down <- venn.diagram(
  x = list(
    iWAT_DOWN = iwat_down,
    gWAT_DOWN = gwat_down
  ),
  category.names = c("iWAT DOWN", "gWAT DOWN"),
  filename = NULL,
  output = TRUE,
  fill = c(col_iwat, col_gwat),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "DOWN-regulated DEGs",
  main.cex = 1.5,
  main.pos = c(0.5, 1.05)
)
grid.draw(venn_down)
popViewport()
dev.off()
cat("[OK] Saved Venn diagram (DOWN DEGs): ", venn_down_file, "\n", sep = "")

## 9) LogFC correlation scatter plot ------------------------
cat("\n=== CREATING LOGFC CORRELATION PLOT ===\n")

## Merge DE tables on gene names
de_merged <- merge(
  de_iwat[, c("logFC", "padj")],
  de_gwat[, c("logFC", "padj")],
  by = "row.names",
  suffixes = c("_iWAT", "_gWAT")
)
rownames(de_merged) <- de_merged$Row.names
de_merged$Row.names <- NULL

## Add significance categories (using depot-specific thresholds)
de_merged <- de_merged %>%
  mutate(
    sig_cat = case_when(
      padj_iWAT < iwat_thresh$fdr_cut & abs(logFC_iWAT) > iwat_thresh$logFC_cut &
        padj_gWAT < gwat_thresh$fdr_cut & abs(logFC_gWAT) > gwat_thresh$logFC_cut ~ "Both DEG",
      padj_iWAT < iwat_thresh$fdr_cut & abs(logFC_iWAT) > iwat_thresh$logFC_cut ~ "iWAT only",
      padj_gWAT < gwat_thresh$fdr_cut & abs(logFC_gWAT) > gwat_thresh$logFC_cut ~ "gWAT only",
      TRUE ~ "Neither"
    ),
    gene_name = rownames(de_merged)
  )

## Calculate correlation
cor_pearson <- cor(de_merged$logFC_iWAT, de_merged$logFC_gWAT, use = "complete.obs")
cor_spearman <- cor(de_merged$logFC_iWAT, de_merged$logFC_gWAT, method = "spearman", use = "complete.obs")

cat("[INFO] Pearson correlation: ", round(cor_pearson, 3), "\n", sep = "")
cat("[INFO] Spearman correlation: ", round(cor_spearman, 3), "\n", sep = "")

## Label interesting genes (extreme in both or opposite)
label_genes <- de_merged %>%
  dplyr::filter(
    sig_cat == "Both DEG" &
    (abs(logFC_iWAT) > 2 | abs(logFC_gWAT) > 2 | sign(logFC_iWAT) != sign(logFC_gWAT))
  ) %>%
  dplyr::slice_max(order_by = abs(logFC_iWAT) + abs(logFC_gWAT), n = 20)

## Create scatter plot
p_scatter <- ggplot(de_merged, aes(x = logFC_iWAT, y = logFC_gWAT, color = sig_cat)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.3) +
  geom_hline(yintercept = c(-logFC_cut, logFC_cut), linetype = "dashed", color = "gray40", linewidth = 0.3) +
  geom_vline(xintercept = c(-logFC_cut, logFC_cut), linetype = "dashed", color = "gray40", linewidth = 0.3) +
  scale_color_manual(
    values = c(
      "Both DEG" = col_both,
      "iWAT only" = col_iwat,
      "gWAT only" = col_gwat,
      "Neither" = col_neither
    ),
    name = "DEG Status"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_name),
    color = "black",
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "gray40",
    segment.size = 0.3
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste0("Pearson r = ", round(cor_pearson, 3), "\nSpearman \u03c1 = ", round(cor_spearman, 3)),
    hjust = -0.1, vjust = 1.5,
    size = 4,
    fontface = "bold",
    color = "black"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "logFC Correlation: iWAT vs gWAT",
    x = "log2 Fold Change (iWAT KD vs CTL)",
    y = "log2 Fold Change (gWAT KD vs CTL)"
  ) +
  coord_fixed(ratio = 1)

scatter_file <- paste0("Scatter_logFC_correlation_", run_tag, ".png")
ggsave(file.path(outdir, "plots", scatter_file), plot = p_scatter, width = 10, height = 8, dpi = 300)
cat("[OK] Saved logFC correlation plot: ", scatter_file, "\n", sep = "")

## 10) Export gene lists for follow-up ---------------------
cat("\n=== EXPORTING GENE LISTS ===\n")

## Shared genes (concordant)
shared_up <- data.frame(
  gene = overlap_up_up,
  logFC_iWAT = de_iwat[overlap_up_up, "logFC"],
  logFC_gWAT = de_gwat[overlap_up_up, "logFC"],
  padj_iWAT = de_iwat[overlap_up_up, "padj"],
  padj_gWAT = de_gwat[overlap_up_up, "padj"],
  stringsAsFactors = FALSE
)
shared_up <- shared_up[order(shared_up$padj_iWAT), ]

shared_down <- data.frame(
  gene = overlap_down_down,
  logFC_iWAT = de_iwat[overlap_down_down, "logFC"],
  logFC_gWAT = de_gwat[overlap_down_down, "logFC"],
  padj_iWAT = de_iwat[overlap_down_down, "padj"],
  padj_gWAT = de_gwat[overlap_down_down, "padj"],
  stringsAsFactors = FALSE
)
shared_down <- shared_down[order(shared_down$padj_iWAT), ]

## Depot-specific genes
iwat_specific_df <- data.frame(
  gene = iwat_specific,
  logFC = de_iwat[iwat_specific, "logFC"],
  padj = de_iwat[iwat_specific, "padj"],
  direction = ifelse(de_iwat[iwat_specific, "logFC"] > 0, "UP", "DOWN"),
  stringsAsFactors = FALSE
)
iwat_specific_df <- iwat_specific_df[order(iwat_specific_df$padj), ]

gwat_specific_df <- data.frame(
  gene = gwat_specific,
  logFC = de_gwat[gwat_specific, "logFC"],
  padj = de_gwat[gwat_specific, "padj"],
  direction = ifelse(de_gwat[gwat_specific, "logFC"] > 0, "UP", "DOWN"),
  stringsAsFactors = FALSE
)
gwat_specific_df <- gwat_specific_df[order(gwat_specific_df$padj), ]

## Save gene lists
write.csv(shared_up, file.path(outdir, "tables", paste0("Shared_UP_genes_", run_tag, ".csv")), row.names = FALSE)
write.csv(shared_down, file.path(outdir, "tables", paste0("Shared_DOWN_genes_", run_tag, ".csv")), row.names = FALSE)
write.csv(iwat_specific_df, file.path(outdir, "tables", paste0("iWAT_specific_genes_", run_tag, ".csv")), row.names = FALSE)
write.csv(gwat_specific_df, file.path(outdir, "tables", paste0("gWAT_specific_genes_", run_tag, ".csv")), row.names = FALSE)

cat("[OK] Saved shared UP genes: ", nrow(shared_up), " genes\n", sep = "")
cat("[OK] Saved shared DOWN genes: ", nrow(shared_down), " genes\n", sep = "")
cat("[OK] Saved iWAT-specific genes: ", nrow(iwat_specific_df), " genes\n", sep = "")
cat("[OK] Saved gWAT-specific genes: ", nrow(gwat_specific_df), " genes\n", sep = "")

## 11) Summary output --------------------------------------
cat("\n==========================================================\n")
cat("=== OVERLAP ANALYSIS SUMMARY ===\n")
cat("==========================================================\n\n")
cat("Total DEGs:\n")
cat("  iWAT: ", length(iwat_all_deg), " (", length(iwat_up), " up, ", length(iwat_down), " down)\n", sep = "")
cat("  gWAT: ", length(gwat_all_deg), " (", length(gwat_up), " up, ", length(gwat_down), " down)\n", sep = "")
cat("\nOverlap:\n")
cat("  Shared (any direction): ", length(overlap_all), " (",
    round(100 * length(overlap_all) / length(iwat_all_deg), 1), "% of iWAT, ",
    round(100 * length(overlap_all) / length(gwat_all_deg), 1), "% of gWAT)\n", sep = "")
cat("  Concordant UP: ", length(overlap_up_up), "\n", sep = "")
cat("  Concordant DOWN: ", length(overlap_down_down), "\n", sep = "")
cat("  Opposite directions: ", length(overlap_opposite), "\n", sep = "")
cat("\nDepot-specific:\n")
cat("  iWAT-only: ", length(iwat_specific), " (",
    round(100 * length(iwat_specific) / length(iwat_all_deg), 1), "% of iWAT DEGs)\n", sep = "")
cat("  gWAT-only: ", length(gwat_specific), " (",
    round(100 * length(gwat_specific) / length(gwat_all_deg), 1), "% of gWAT DEGs)\n", sep = "")
cat("\nCorrelation:\n")
cat("  Pearson r: ", round(cor_pearson, 3), "\n", sep = "")
cat("  Spearman \u03c1: ", round(cor_spearman, 3), "\n", sep = "")
cat("\n==========================================================\n\n")

cat("======================================\n")
cat("=== PART 1.5 (OVERLAP) COMPLETE ===\n")
cat("======================================\n")
cat("Files created:\n")
cat("  - Overlap summary table\n")
cat("  - 3 Venn diagrams (all, up, down)\n")
cat("  - logFC correlation scatter plot\n")
cat("  - 4 gene lists (shared up/down, depot-specific)\n\n")
