#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 7: DEPOT CONCORDANCE
## =========================================================
## Compares iWAT vs gWAT KAT8-KD responses to distinguish
## shared (depot-general) vs depot-specific transcriptional effects
##
## KEY ANALYSES:
## 1. DEG Overlap: Direction-specific Venn/classification
## 2. logFC Correlation: Scatter of iWAT vs gWAT FCs
## 3. Concordance Classification: Same-direction, opposite, unique
## 4. Shared KAT8 Response: Genes DEG in BOTH depots (depot-general)
## 5. Depot-Specific Genes: DEG in one depot only
##
## REQUIRES:
##   - Part 1 tissue results (savepoints/RUN_*/tables/DE_tissue_*.csv)
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "scales",
                   "RColorBrewer", "VennDiagram", "grid")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
})

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  use_save_core <- TRUE
} else {
  cat("[WARN] save_core not found; proceeding without it\n")
  use_save_core <- FALSE
}

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## ============================================================
## 3) FIND PART 1 RESULTS
## ============================================================
cat("\n=== FINDING TISSUE RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run part1_main_analysis.R first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No run directories found in savepoints/")
}

## Find tissue DE files across all run directories
tissue_de_files <- unlist(lapply(run_dirs, function(d) {
  list.files(file.path(d, "tables"), pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
}))

## Exclude OVERALL contrasts
tissue_de_files <- tissue_de_files[!grepl("OVERALL", tissue_de_files)]

if (length(tissue_de_files) == 0) stop("No tissue DE tables found. Run Part 1 first.")

## Use most recent files
tissue_de_files <- tissue_de_files[order(file.info(tissue_de_files)$mtime, decreasing = TRUE)]

## Load tissue results, deduplicate by contrast name
tissue_tables <- list()
for (f in tissue_de_files) {
  contrast <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  if (!(contrast %in% names(tissue_tables))) {
    tt <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
    ## Skip tables missing required columns
    if (!all(c("log2FoldChange", "padj") %in% colnames(tt)) &&
        !all(c("logFC", "padj") %in% colnames(tt))) {
      cat("[INFO] Skipping '", contrast, "': missing required columns\n", sep = "")
      next
    }
    tissue_tables[[contrast]] <- tt
    cat("[OK] Tissue DE: ", contrast, " (", nrow(tt), " genes)\n", sep = "")
  }
}

## Identify iWAT and gWAT contrasts
iwat_name <- grep("iWAT", names(tissue_tables), value = TRUE)
gwat_name <- grep("gWAT", names(tissue_tables), value = TRUE)

if (length(iwat_name) != 1 || length(gwat_name) != 1) {
  stop("Expected exactly 1 iWAT and 1 gWAT contrast. Found: ",
       paste(names(tissue_tables), collapse = ", "))
}

de_iwat <- tissue_tables[[iwat_name]]
de_gwat <- tissue_tables[[gwat_name]]

cat("\n[OK] iWAT contrast: ", iwat_name, " (", nrow(de_iwat), " genes)\n", sep = "")
cat("[OK] gWAT contrast: ", gwat_name, " (", nrow(de_gwat), " genes)\n", sep = "")

## Detect logFC column name (log2FoldChange from DESeq2 or logFC from earlier processing)
lfc_col_iwat <- if ("log2FoldChange" %in% colnames(de_iwat)) "log2FoldChange" else "logFC"
lfc_col_gwat <- if ("log2FoldChange" %in% colnames(de_gwat)) "log2FoldChange" else "logFC"

## Set up output directory
tissue_run_dir <- dirname(dirname(tissue_de_files[1]))
outdir <- tissue_run_dir
run_tag <- gsub("^RUN_", "", basename(outdir))
tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
concordance_dir <- file.path(plots_dir, "depot_concordance")
if (!dir.exists(concordance_dir)) dir.create(concordance_dir, recursive = TRUE)

cat("[INFO] Output directory: ", outdir, "\n", sep = "")

## ============================================================
## 4) ANALYSIS 1: DEG OVERLAP (iWAT vs gWAT)
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 1: DEG OVERLAP (iWAT vs gWAT) ===\n")
cat("============================================================\n\n")

## Depot-specific thresholds from parameters.R
iwat_thresh <- get_tissue_thresholds("iWAT")
gwat_thresh <- get_tissue_thresholds("gWAT")

## iWAT DEGs
iwat_up <- rownames(de_iwat)[!is.na(de_iwat$padj) &
  de_iwat$padj < iwat_thresh$fdr_cut &
  de_iwat[[lfc_col_iwat]] > iwat_thresh$logFC_cut]
iwat_down <- rownames(de_iwat)[!is.na(de_iwat$padj) &
  de_iwat$padj < iwat_thresh$fdr_cut &
  de_iwat[[lfc_col_iwat]] < -iwat_thresh$logFC_cut]
iwat_all_deg <- union(iwat_up, iwat_down)

## gWAT DEGs
gwat_up <- rownames(de_gwat)[!is.na(de_gwat$padj) &
  de_gwat$padj < gwat_thresh$fdr_cut &
  de_gwat[[lfc_col_gwat]] > gwat_thresh$logFC_cut]
gwat_down <- rownames(de_gwat)[!is.na(de_gwat$padj) &
  de_gwat$padj < gwat_thresh$fdr_cut &
  de_gwat[[lfc_col_gwat]] < -gwat_thresh$logFC_cut]
gwat_all_deg <- union(gwat_up, gwat_down)

cat("iWAT DEGs: ", length(iwat_all_deg), " (", length(iwat_up), " up, ",
    length(iwat_down), " down)\n", sep = "")
cat("gWAT DEGs: ", length(gwat_all_deg), " (", length(gwat_up), " up, ",
    length(gwat_down), " down)\n", sep = "")

## Overall overlap
shared_all <- intersect(iwat_all_deg, gwat_all_deg)
iwat_only  <- setdiff(iwat_all_deg, gwat_all_deg)
gwat_only  <- setdiff(gwat_all_deg, iwat_all_deg)

cat("\nShared DEGs (any direction): ", length(shared_all), "\n", sep = "")

## Direction-specific overlap
shared_up_up     <- intersect(iwat_up, gwat_up)
shared_down_down <- intersect(iwat_down, gwat_down)
shared_concordant <- union(shared_up_up, shared_down_down)

shared_up_down   <- intersect(iwat_up, gwat_down)  # UP in iWAT, DOWN in gWAT
shared_down_up   <- intersect(iwat_down, gwat_up)   # DOWN in iWAT, UP in gWAT
shared_discordant <- union(shared_up_down, shared_down_up)

cat("  Concordant (same direction): ", length(shared_concordant), "\n", sep = "")
cat("    Both up:   ", length(shared_up_up), "\n", sep = "")
cat("    Both down: ", length(shared_down_down), "\n", sep = "")
cat("  Discordant (opposite direction): ", length(shared_discordant), "\n", sep = "")
cat("    UP iWAT / DOWN gWAT: ", length(shared_up_down), "\n", sep = "")
cat("    DOWN iWAT / UP gWAT: ", length(shared_down_up), "\n", sep = "")
cat("  iWAT-only DEGs: ", length(iwat_only), "\n", sep = "")
cat("  gWAT-only DEGs: ", length(gwat_only), "\n", sep = "")

concordance_pct <- if (length(shared_all) > 0) {
  round(100 * length(shared_concordant) / length(shared_all), 1)
} else {
  NA_real_
}
cat("  Concordance rate: ", concordance_pct, "%\n", sep = "")

## Fisher's exact test for overlap enrichment
common_genes <- intersect(rownames(de_iwat), rownames(de_gwat))
a  <- length(shared_all)
b  <- length(setdiff(gwat_all_deg, iwat_all_deg))
cc <- length(setdiff(iwat_all_deg, gwat_all_deg))
d  <- length(common_genes) - a - b - cc

fisher_res <- fisher.test(matrix(c(a, b, cc, d), nrow = 2))
cat("  Fisher's exact test: OR = ", round(fisher_res$estimate, 2),
    ", p = ", format(fisher_res$p.value, digits = 3), "\n", sep = "")

## ============================================================
## 5) VENN DIAGRAMS
## ============================================================
cat("\n=== CREATING VENN DIAGRAMS ===\n")

col_iwat <- "#4575B4"
col_gwat <- "#D73027"

## All DEGs Venn
venn_all_file <- file.path(concordance_dir, paste0("Venn_depot_all_DEGs_", run_tag, ".png"))
png(venn_all_file, width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_all <- venn.diagram(
  x = list(iWAT = iwat_all_deg, gWAT = gwat_all_deg),
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
  main = paste0("All DEGs (iWAT: FDR<", iwat_thresh$fdr_cut,
                " |logFC|>", iwat_thresh$logFC_cut,
                "; gWAT: FDR<", gwat_thresh$fdr_cut,
                " |logFC|>", gwat_thresh$logFC_cut, ")"),
  main.cex = 1.2,
  main.pos = c(0.5, 1.05)
)
grid.draw(venn_all)
popViewport()
dev.off()
cat("[OK] Saved Venn (all DEGs)\n")

## UP DEGs Venn
venn_up_file <- file.path(concordance_dir, paste0("Venn_depot_UP_DEGs_", run_tag, ".png"))
png(venn_up_file, width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_up <- venn.diagram(
  x = list(iWAT_UP = iwat_up, gWAT_UP = gwat_up),
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
cat("[OK] Saved Venn (UP DEGs)\n")

## DOWN DEGs Venn
venn_down_file <- file.path(concordance_dir, paste0("Venn_depot_DOWN_DEGs_", run_tag, ".png"))
png(venn_down_file, width = 2400, height = 2400, res = 300)
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
venn_down <- venn.diagram(
  x = list(iWAT_DOWN = iwat_down, gWAT_DOWN = gwat_down),
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
cat("[OK] Saved Venn (DOWN DEGs)\n")

## ============================================================
## 6) ANALYSIS 2: logFC CORRELATION
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 2: logFC CORRELATION ===\n")
cat("============================================================\n\n")

## Build correlation data frame from common genes
cor_df <- data.frame(
  gene       = common_genes,
  iwat_lfc   = de_iwat[common_genes, lfc_col_iwat],
  gwat_lfc   = de_gwat[common_genes, lfc_col_gwat],
  iwat_padj  = de_iwat[common_genes, "padj"],
  gwat_padj  = de_gwat[common_genes, "padj"],
  stringsAsFactors = FALSE
)

## Remove NAs
cor_df <- cor_df %>% filter(is.finite(iwat_lfc), is.finite(gwat_lfc))

## Classify each gene
cor_df <- cor_df %>%
  mutate(
    iwat_sig = !is.na(iwat_padj) & iwat_padj < iwat_thresh$fdr_cut &
      abs(iwat_lfc) > iwat_thresh$logFC_cut,
    gwat_sig = !is.na(gwat_padj) & gwat_padj < gwat_thresh$fdr_cut &
      abs(gwat_lfc) > gwat_thresh$logFC_cut,
    category = case_when(
      iwat_sig & gwat_sig & sign(iwat_lfc) == sign(gwat_lfc) ~ "Concordant",
      iwat_sig & gwat_sig & sign(iwat_lfc) != sign(gwat_lfc) ~ "Discordant",
      iwat_sig & !gwat_sig ~ "iWAT only",
      !iwat_sig & gwat_sig ~ "gWAT only",
      TRUE ~ "NS"
    )
  )

## Correlation (all genes)
cor_method <- depot_concordance_params$cor_method
cor_all <- cor(cor_df$iwat_lfc, cor_df$gwat_lfc, method = cor_method)
cor_test_all <- cor.test(cor_df$iwat_lfc, cor_df$gwat_lfc, method = cor_method)

## Correlation (DEGs only)
sig_df <- cor_df %>% filter(iwat_sig | gwat_sig)
if (nrow(sig_df) > 10) {
  cor_sig <- cor(sig_df$iwat_lfc, sig_df$gwat_lfc, method = cor_method)
} else {
  cor_sig <- NA
}

## Spearman as secondary metric
cor_spearman <- cor(cor_df$iwat_lfc, cor_df$gwat_lfc, method = "spearman")

cat("Correlation (all genes):  r = ", round(cor_all, 3),
    " (p = ", format(cor_test_all$p.value, digits = 3), ")\n", sep = "")
cat("Correlation (DEGs only):  r = ", round(cor_sig, 3), "\n", sep = "")
cat("Spearman (all genes):     rho = ", round(cor_spearman, 3), "\n", sep = "")

## Save concordance gene table (non-NS genes)
cor_output <- cor_df %>%
  filter(category != "NS") %>%
  arrange(desc(abs(iwat_lfc) + abs(gwat_lfc)))

cor_file <- paste0("depot_concordance_genes_", run_tag, ".csv")
write.csv(cor_output, file = file.path(tables_dir, cor_file), row.names = FALSE)
cat("[OK] Saved: ", cor_file, "\n", sep = "")

## --- Scatter plot ---
top_n_label <- depot_concordance_params$top_n_label

label_df <- cor_df %>%
  filter(category %in% c("Concordant", "Discordant")) %>%
  arrange(desc(abs(iwat_lfc) + abs(gwat_lfc))) %>%
  slice_head(n = top_n_label)

category_colors <- c(
  "Concordant" = "#4CAF50",
  "Discordant" = "#F44336",
  "iWAT only"  = col_iwat,
  "gWAT only"  = col_gwat,
  "NS"         = "grey85"
)

## Determine axis limits: 99.5th percentile to avoid outlier squishing
axis_max <- max(
  quantile(abs(cor_df$iwat_lfc), 0.995, na.rm = TRUE),
  quantile(abs(cor_df$gwat_lfc), 0.995, na.rm = TRUE)
) * 1.05
n_clipped <- sum(abs(cor_df$iwat_lfc) > axis_max | abs(cor_df$gwat_lfc) > axis_max,
                 na.rm = TRUE)
if (n_clipped > 0) cat("[INFO] Clipping ", n_clipped, " extreme points for readability\n", sep = "")

p_scatter <- ggplot(cor_df, aes(x = iwat_lfc, y = gwat_lfc, color = category)) +
  geom_point(data = cor_df %>% filter(category == "NS"), size = 0.5, alpha = 0.3) +
  geom_point(data = cor_df %>% filter(category != "NS"), size = 1.5, alpha = 0.7) +
  geom_text_repel(
    data = label_df,
    aes(label = gene),
    color = "black", size = 3, max.overlaps = 20,
    segment.color = "grey50", segment.size = 0.3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = category_colors, name = "Category") +
  labs(
    title = "logFC Concordance: iWAT vs gWAT (KAT8 KD vs CTL)",
    subtitle = paste0(
      "Pearson r = ", round(cor_all, 3), " (all genes); r = ",
      round(cor_sig, 3), " (DEGs); Spearman \u03c1 = ", round(cor_spearman, 3)
    ),
    x = expression(log[2]~FC~"(iWAT KD vs CTL)"),
    y = expression(log[2]~FC~"(gWAT KD vs CTL)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right"
  ) +
  coord_cartesian(xlim = c(-axis_max, axis_max), ylim = c(-axis_max, axis_max))

scatter_file <- paste0("Concordance_scatter_depot_", run_tag, ".png")
ggsave(file.path(concordance_dir, scatter_file), plot = p_scatter,
       width = depot_concordance_params$plot_width,
       height = depot_concordance_params$plot_height,
       dpi = 300, bg = "white")
cat("[OK] Saved: ", scatter_file, "\n", sep = "")

## ============================================================
## 7) ANALYSIS 3: SHARED KAT8 RESPONSE (Depot-General Signature)
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 3: SHARED KAT8 RESPONSE (DEPOT-GENERAL) ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  Genes DEG in BOTH iWAT and gWAT with the same direction are\n")
cat("  likely core KAT8-dependent targets conserved across depots,\n")
cat("  rather than depot-specific microenvironment effects.\n\n")

if (length(shared_concordant) == 0) {
  cat("  No concordant DEGs found between depots\n\n")
} else {
  ## Build shared response gene table
  shared_df <- data.frame(
    gene      = shared_concordant,
    iwat_lfc  = de_iwat[shared_concordant, lfc_col_iwat],
    iwat_padj = de_iwat[shared_concordant, "padj"],
    gwat_lfc  = de_gwat[shared_concordant, lfc_col_gwat],
    gwat_padj = de_gwat[shared_concordant, "padj"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      direction = ifelse(iwat_lfc > 0, "Up", "Down"),
      mean_lfc  = (iwat_lfc + gwat_lfc) / 2,
      lfc_diff  = abs(iwat_lfc - gwat_lfc)
    ) %>%
    arrange(desc(abs(mean_lfc)))

  cat("Shared depot-general genes: ", nrow(shared_df), "\n", sep = "")
  cat("  Up in both:   ", sum(shared_df$direction == "Up"), "\n", sep = "")
  cat("  Down in both: ", sum(shared_df$direction == "Down"), "\n", sep = "")
  cat("  Mean |logFC| difference between depots: ", round(mean(shared_df$lfc_diff), 2), "\n", sep = "")

  ## Show top hits
  cat("\nTop shared KAT8-response genes (by mean |logFC|):\n")
  for (i in seq_len(min(15, nrow(shared_df)))) {
    row <- shared_df[i, ]
    cat("  ", row$gene, ": iWAT=", round(row$iwat_lfc, 2),
        " gWAT=", round(row$gwat_lfc, 2), " (", row$direction, ")\n", sep = "")
  }

  ## Save
  shared_file <- paste0("depot_shared_response_genes_", run_tag, ".csv")
  write.csv(shared_df, file = file.path(tables_dir, shared_file), row.names = FALSE)
  cat("\n[OK] Saved: ", shared_file, "\n", sep = "")
}

## ============================================================
## 8) ANALYSIS 4: DISCORDANT GENES (Opposite Direction)
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 4: DISCORDANT GENES (OPPOSITE DIRECTION) ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  Genes that are DEG in both depots but in OPPOSITE directions\n")
cat("  may reflect depot-divergent regulatory mechanisms downstream\n")
cat("  of KAT8 or different cell composition effects.\n\n")

if (length(shared_discordant) == 0) {
  cat("  No discordant DEGs found\n\n")
} else {
  discord_df <- data.frame(
    gene      = shared_discordant,
    iwat_lfc  = de_iwat[shared_discordant, lfc_col_iwat],
    iwat_padj = de_iwat[shared_discordant, "padj"],
    gwat_lfc  = de_gwat[shared_discordant, lfc_col_gwat],
    gwat_padj = de_gwat[shared_discordant, "padj"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      iwat_direction = ifelse(iwat_lfc > 0, "Up", "Down"),
      gwat_direction = ifelse(gwat_lfc > 0, "Up", "Down"),
      pattern = paste0("iWAT_", iwat_direction, "__gWAT_", gwat_direction)
    ) %>%
    arrange(desc(abs(iwat_lfc) + abs(gwat_lfc)))

  cat("Discordant genes: ", nrow(discord_df), "\n", sep = "")
  cat("  UP in iWAT / DOWN in gWAT: ", sum(discord_df$pattern == "iWAT_Up__gWAT_Down"), "\n", sep = "")
  cat("  DOWN in iWAT / UP in gWAT: ", sum(discord_df$pattern == "iWAT_Down__gWAT_Up"), "\n", sep = "")

  ## Show all or top hits
  cat("\nDiscordant genes (by combined |logFC|):\n")
  for (i in seq_len(min(15, nrow(discord_df)))) {
    row <- discord_df[i, ]
    cat("  ", row$gene, ": iWAT=", round(row$iwat_lfc, 2),
        " gWAT=", round(row$gwat_lfc, 2), "\n", sep = "")
  }

  discord_file <- paste0("depot_discordant_genes_", run_tag, ".csv")
  write.csv(discord_df, file = file.path(tables_dir, discord_file), row.names = FALSE)
  cat("\n[OK] Saved: ", discord_file, "\n", sep = "")
}

## ============================================================
## 9) ANALYSIS 5: DEPOT-SPECIFIC GENES
## ============================================================
cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 5: DEPOT-SPECIFIC GENES ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  Genes DEG in only one depot may reflect:\n")
cat("  - Depot-specific adipocyte biology (e.g., browning in iWAT)\n")
cat("  - Differences in immune/stromal cell composition\n")
cat("  - Depot-specific KAT8 target accessibility\n\n")

## --- iWAT-specific ---
cat("--- iWAT-SPECIFIC GENES ---\n")
cat("  iWAT-only DEGs: ", length(iwat_only), "\n", sep = "")

if (length(iwat_only) > 0) {
  iwat_only_df <- data.frame(
    gene      = iwat_only,
    iwat_lfc  = de_iwat[iwat_only, lfc_col_iwat],
    iwat_padj = de_iwat[iwat_only, "padj"],
    direction = ifelse(de_iwat[iwat_only, lfc_col_iwat] > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )

  ## Cross-reference with gWAT data (what does gWAT show for these genes?)
  common_iwat_only <- intersect(iwat_only, rownames(de_gwat))
  iwat_only_df$gwat_lfc  <- NA_real_
  iwat_only_df$gwat_padj <- NA_real_
  if (length(common_iwat_only) > 0) {
    iwat_only_df[iwat_only_df$gene %in% common_iwat_only, "gwat_lfc"] <-
      de_gwat[common_iwat_only, lfc_col_gwat]
    iwat_only_df[iwat_only_df$gene %in% common_iwat_only, "gwat_padj"] <-
      de_gwat[common_iwat_only, "padj"]
  }

  iwat_only_df <- iwat_only_df %>% arrange(desc(abs(iwat_lfc)))

  cat("    Up:   ", sum(iwat_only_df$direction == "Up"), "\n", sep = "")
  cat("    Down: ", sum(iwat_only_df$direction == "Down"), "\n", sep = "")

  iwat_only_file <- paste0("depot_iWAT_specific_genes_", run_tag, ".csv")
  write.csv(iwat_only_df, file = file.path(tables_dir, iwat_only_file), row.names = FALSE)
  cat("  [OK] Saved: ", iwat_only_file, "\n\n", sep = "")
}

## --- gWAT-specific ---
cat("--- gWAT-SPECIFIC GENES ---\n")
cat("  gWAT-only DEGs: ", length(gwat_only), "\n", sep = "")

if (length(gwat_only) > 0) {
  gwat_only_df <- data.frame(
    gene      = gwat_only,
    gwat_lfc  = de_gwat[gwat_only, lfc_col_gwat],
    gwat_padj = de_gwat[gwat_only, "padj"],
    direction = ifelse(de_gwat[gwat_only, lfc_col_gwat] > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )

  ## Cross-reference with iWAT data
  common_gwat_only <- intersect(gwat_only, rownames(de_iwat))
  gwat_only_df$iwat_lfc  <- NA_real_
  gwat_only_df$iwat_padj <- NA_real_
  if (length(common_gwat_only) > 0) {
    gwat_only_df[gwat_only_df$gene %in% common_gwat_only, "iwat_lfc"] <-
      de_iwat[common_gwat_only, lfc_col_iwat]
    gwat_only_df[gwat_only_df$gene %in% common_gwat_only, "iwat_padj"] <-
      de_iwat[common_gwat_only, "padj"]
  }

  gwat_only_df <- gwat_only_df %>% arrange(desc(abs(gwat_lfc)))

  cat("    Up:   ", sum(gwat_only_df$direction == "Up"), "\n", sep = "")
  cat("    Down: ", sum(gwat_only_df$direction == "Down"), "\n", sep = "")

  gwat_only_file <- paste0("depot_gWAT_specific_genes_", run_tag, ".csv")
  write.csv(gwat_only_df, file = file.path(tables_dir, gwat_only_file), row.names = FALSE)
  cat("  [OK] Saved: ", gwat_only_file, "\n\n", sep = "")
}

## ============================================================
## 10) SUMMARY TABLE
## ============================================================
cat("\n=== GENERATING SUMMARY ===\n")

summary_df <- data.frame(
  Metric = c(
    "iWAT total DEGs",
    "iWAT UP",
    "iWAT DOWN",
    "gWAT total DEGs",
    "gWAT UP",
    "gWAT DOWN",
    "Shared DEGs (any direction)",
    "Shared concordant (same direction)",
    "  Both UP",
    "  Both DOWN",
    "Shared discordant (opposite direction)",
    "  UP iWAT / DOWN gWAT",
    "  DOWN iWAT / UP gWAT",
    "iWAT-only DEGs",
    "gWAT-only DEGs",
    "Concordance rate (%)",
    "Pct iWAT DEGs shared",
    "Pct gWAT DEGs shared",
    paste0("Pearson r (all genes, ", cor_method, ")"),
    "Pearson r (DEGs only)",
    "Spearman rho (all genes)",
    "Fisher OR (overlap enrichment)",
    "Fisher p-value"
  ),
  Value = c(
    length(iwat_all_deg),
    length(iwat_up),
    length(iwat_down),
    length(gwat_all_deg),
    length(gwat_up),
    length(gwat_down),
    length(shared_all),
    length(shared_concordant),
    length(shared_up_up),
    length(shared_down_down),
    length(shared_discordant),
    length(shared_up_down),
    length(shared_down_up),
    length(iwat_only),
    length(gwat_only),
    concordance_pct,
    round(100 * length(shared_all) / max(1, length(iwat_all_deg)), 1),
    round(100 * length(shared_all) / max(1, length(gwat_all_deg)), 1),
    round(cor_all, 4),
    round(cor_sig, 4),
    round(cor_spearman, 4),
    round(fisher_res$estimate, 2),
    format(fisher_res$p.value, digits = 4)
  ),
  stringsAsFactors = FALSE
)

summary_file <- paste0("depot_concordance_summary_", run_tag, ".csv")
write.csv(summary_df, file = file.path(tables_dir, summary_file), row.names = FALSE)
cat("[OK] Saved summary: ", summary_file, "\n", sep = "")

## Print summary
cat("\n")
cat("============================================================\n")
cat("=== DEPOT CONCORDANCE SUMMARY ===\n")
cat("============================================================\n")
print(summary_df, row.names = FALSE)
cat("\n")

## ============================================================
## 11) SESSION INFO
## ============================================================
log_dir <- file.path(outdir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
session_file <- file.path(log_dir, paste0("sessionInfo_depot_concordance_", run_tag, ".txt"))
sink(session_file)
cat("=== R SESSION INFORMATION - Part 7 (Depot Concordance) ===\n\n")
print(sessionInfo())
sink()
cat("[INFO] Session info saved\n")

cat("\n=== Part 7 (Depot Concordance) COMPLETE ===\n")
cat("Output directory: ", concordance_dir, "\n", sep = "")
cat("\nFiles created:\n")
cat("  - 3 Venn diagrams (all, UP, DOWN)\n")
cat("  - logFC correlation scatter plot\n")
cat("  - Concordance gene table (all non-NS genes)\n")
cat("  - Shared depot-general response genes\n")
cat("  - Discordant genes (opposite direction)\n")
cat("  - iWAT-specific genes (with gWAT cross-reference)\n")
cat("  - gWAT-specific genes (with iWAT cross-reference)\n")
cat("  - Summary table\n\n")
