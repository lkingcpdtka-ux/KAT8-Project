#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 7: IN VITRO / IN VIVO CONCORDANCE
## =========================================================
## Compares cell culture (3T3-L1) vs tissue (iWAT/gWAT) results
## to distinguish cell-autonomous vs microenvironment-driven effects
##
## KEY ANALYSES:
## 1. DEG Overlap: Venn/UpSet of shared vs unique DEGs
## 2. logFC Correlation: Scatter of tissue vs cell FCs
## 3. Concordance Classification: Same-direction, opposite, unique
## 4. Pathway-Level Comparison: Shared enriched pathways
## 5. Cell-Autonomous Signature: Genes changed in BOTH systems
##
## REQUIRES:
##   - Part 1 tissue results (savepoints/RUN_*/tables/DE_tissue_*.csv)
##   - Cell culture results (savepoints/RUN_*/tables/DE_cells_*.csv)
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "scales", "RColorBrewer", "VennDiagram", "grid")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
  library(UpSetR)
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
## 3) FIND RESULTS FROM BOTH ANALYSES
## ============================================================

cat("\n=== FINDING TISSUE AND CELL CULTURE RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run Part 1 and cell culture analysis first.")
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

## Find cell culture DE files
cell_de_files <- unlist(lapply(run_dirs, function(d) {
  list.files(file.path(d, "tables"), pattern = "^DE_cells_.*\\.csv$", full.names = TRUE)
}))

if (length(tissue_de_files) == 0) stop("No tissue DE tables found. Run Part 1 first.")
if (length(cell_de_files) == 0) stop("No cell culture DE tables found. Run KAT8_bulk_cells_comprehensive.R first.")

## Use most recent files
tissue_de_files <- tissue_de_files[order(file.info(tissue_de_files)$mtime, decreasing = TRUE)]
cell_de_files <- cell_de_files[order(file.info(cell_de_files)$mtime, decreasing = TRUE)]

## Load cell culture results (most recent)
cell_de <- read.csv(cell_de_files[1], row.names = 1, stringsAsFactors = FALSE)
cat("[OK] Cell culture DE: ", nrow(cell_de), " genes from ", basename(cell_de_files[1]), "\n", sep = "")

## Load tissue results (separate iWAT and gWAT)
tissue_tables <- list()
for (f in tissue_de_files) {
  contrast <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  if (!(contrast %in% names(tissue_tables))) {
    tissue_tables[[contrast]] <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
    cat("[OK] Tissue DE: ", contrast, " (", nrow(tissue_tables[[contrast]]), " genes)\n", sep = "")
  }
}

## Set up output directory (use most recent tissue run)
tissue_run_dir <- dirname(dirname(tissue_de_files[1]))
outdir <- tissue_run_dir
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")
concordance_dir <- file.path(plots_dir, "concordance")
if (!dir.exists(concordance_dir)) dir.create(concordance_dir, recursive = TRUE)

cat("[INFO] Output directory: ", outdir, "\n", sep = "")

## ============================================================
## 4) ANALYSIS 1: GENE-LEVEL OVERLAP
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 1: DEG OVERLAP (CELL CULTURE vs TISSUE) ===\n")
cat("============================================================\n\n")

## Extract thresholds from parameters
lfc_tissue <- concordance_params$logFC_cut_tissue
fdr_tissue <- concordance_params$fdr_cut_tissue
lfc_cells  <- concordance_params$logFC_cut_cells
fdr_cells  <- concordance_params$fdr_cut_cells

## Get cell culture DEGs
cell_deg_up <- rownames(cell_de)[!is.na(cell_de$padj) & cell_de$padj < fdr_cells &
                                   cell_de$log2FoldChange > lfc_cells]
cell_deg_down <- rownames(cell_de)[!is.na(cell_de$padj) & cell_de$padj < fdr_cells &
                                     cell_de$log2FoldChange < -lfc_cells]
cell_deg_all <- union(cell_deg_up, cell_deg_down)

cat("Cell culture DEGs: ", length(cell_deg_all), " (", length(cell_deg_up), " up, ",
    length(cell_deg_down), " down)\n", sep = "")

## Get tissue DEGs for each contrast
overlap_results <- list()

for (contrast in names(tissue_tables)) {
  tt <- tissue_tables[[contrast]]

  tissue_deg_up <- rownames(tt)[!is.na(tt$padj) & tt$padj < fdr_tissue &
                                  tt$log2FoldChange > lfc_tissue]
  tissue_deg_down <- rownames(tt)[!is.na(tt$padj) & tt$padj < fdr_tissue &
                                    tt$log2FoldChange < -lfc_tissue]
  tissue_deg_all <- union(tissue_deg_up, tissue_deg_down)

  cat("\n--- ", contrast, " ---\n", sep = "")
  cat("Tissue DEGs: ", length(tissue_deg_all), " (", length(tissue_deg_up), " up, ",
      length(tissue_deg_down), " down)\n", sep = "")

  ## Overall overlap
  shared_all <- intersect(cell_deg_all, tissue_deg_all)
  cat("Shared DEGs (any direction): ", length(shared_all), "\n", sep = "")

  ## Direction-specific overlap
  shared_up_up <- intersect(cell_deg_up, tissue_deg_up)
  shared_down_down <- intersect(cell_deg_down, tissue_deg_down)
  shared_concordant <- union(shared_up_up, shared_down_down)

  shared_up_down <- intersect(cell_deg_up, tissue_deg_down)
  shared_down_up <- intersect(cell_deg_down, tissue_deg_up)
  shared_discordant <- union(shared_up_down, shared_down_up)

  cat("  Concordant (same direction): ", length(shared_concordant), "\n", sep = "")
  cat("    Both up:   ", length(shared_up_up), "\n", sep = "")
  cat("    Both down: ", length(shared_down_down), "\n", sep = "")
  cat("  Discordant (opposite direction): ", length(shared_discordant), "\n", sep = "")

  cell_only <- setdiff(cell_deg_all, tissue_deg_all)
  tissue_only <- setdiff(tissue_deg_all, cell_deg_all)
  cat("  Cell-only DEGs: ", length(cell_only), "\n", sep = "")
  cat("  Tissue-only DEGs: ", length(tissue_only), "\n", sep = "")

  overlap_results[[contrast]] <- list(
    shared_all = shared_all,
    shared_concordant = shared_concordant,
    shared_discordant = shared_discordant,
    shared_up_up = shared_up_up,
    shared_down_down = shared_down_down,
    cell_only = cell_only,
    tissue_only = tissue_only,
    tissue_deg_all = tissue_deg_all,
    tissue_deg_up = tissue_deg_up,
    tissue_deg_down = tissue_deg_down
  )

  ## Fisher's exact test for enrichment
  common_genes <- intersect(rownames(cell_de), rownames(tt))
  a <- length(intersect(cell_deg_all, tissue_deg_all))
  b <- length(setdiff(tissue_deg_all, cell_deg_all))
  cc <- length(setdiff(cell_deg_all, tissue_deg_all))
  d <- length(common_genes) - a - b - cc

  fisher_res <- fisher.test(matrix(c(a, b, cc, d), nrow = 2))
  cat("  Fisher's exact test: OR = ", round(fisher_res$estimate, 2),
      ", p = ", format(fisher_res$p.value, digits = 3), "\n", sep = "")
}

## ============================================================
## 5) ANALYSIS 2: logFC CORRELATION
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 2: logFC CORRELATION ===\n")
cat("============================================================\n\n")

for (contrast in names(tissue_tables)) {
  tt <- tissue_tables[[contrast]]

  ## Find common genes
  common_genes <- intersect(rownames(cell_de), rownames(tt))
  cat("--- ", contrast, " ---\n", sep = "")
  cat("Common genes between cell and tissue: ", length(common_genes), "\n", sep = "")

  ## Build correlation data frame
  cor_df <- data.frame(
    gene = common_genes,
    cell_lfc = cell_de[common_genes, "log2FoldChange"],
    tissue_lfc = tt[common_genes, "log2FoldChange"],
    cell_padj = cell_de[common_genes, "padj"],
    tissue_padj = tt[common_genes, "padj"],
    stringsAsFactors = FALSE
  )

  ## Remove NAs
  cor_df <- cor_df %>% filter(is.finite(cell_lfc), is.finite(tissue_lfc))

  ## Classify each gene
  cor_df <- cor_df %>%
    mutate(
      cell_sig = !is.na(cell_padj) & cell_padj < fdr_cells & abs(cell_lfc) > lfc_cells,
      tissue_sig = !is.na(tissue_padj) & tissue_padj < fdr_tissue & abs(tissue_lfc) > lfc_tissue,
      category = case_when(
        cell_sig & tissue_sig & sign(cell_lfc) == sign(tissue_lfc) ~ "Concordant",
        cell_sig & tissue_sig & sign(cell_lfc) != sign(tissue_lfc) ~ "Discordant",
        cell_sig & !tissue_sig ~ "Cell only",
        !cell_sig & tissue_sig ~ "Tissue only",
        TRUE ~ "NS"
      )
    )

  ## Correlation (all genes)
  cor_all <- cor(cor_df$cell_lfc, cor_df$tissue_lfc, method = concordance_params$cor_method)
  cor_test_all <- cor.test(cor_df$cell_lfc, cor_df$tissue_lfc, method = concordance_params$cor_method)

  ## Correlation (significant genes only)
  sig_df <- cor_df %>% filter(cell_sig | tissue_sig)
  if (nrow(sig_df) > 10) {
    cor_sig <- cor(sig_df$cell_lfc, sig_df$tissue_lfc, method = concordance_params$cor_method)
  } else {
    cor_sig <- NA
  }

  cat("  Correlation (all genes):     r = ", round(cor_all, 3),
      " (p = ", format(cor_test_all$p.value, digits = 3), ")\n", sep = "")
  cat("  Correlation (DEGs only):     r = ", round(cor_sig, 3), "\n", sep = "")

  ## Save correlation table
  cor_output <- cor_df %>%
    filter(category != "NS") %>%
    arrange(desc(abs(cell_lfc) + abs(tissue_lfc)))

  cor_file <- paste0("concordance_genes_", contrast, "_", run_tag, ".csv")
  write.csv(cor_output, file = file.path(tables_dir, cor_file), row.names = FALSE)
  cat("  [OK] Saved: ", cor_file, "\n", sep = "")

  ## --- Scatter plot ---
  label_df <- cor_df %>%
    filter(category %in% c("Concordant", "Discordant")) %>%
    arrange(desc(abs(cell_lfc) + abs(tissue_lfc))) %>%
    slice_head(n = concordance_params$top_n_label)

  category_colors <- c(
    "Concordant" = "#4CAF50",
    "Discordant" = "#F44336",
    "Cell only"  = "#FF9800",
    "Tissue only" = "#2196F3",
    "NS"         = "grey85"
  )

  p_scatter <- ggplot(cor_df, aes(x = cell_lfc, y = tissue_lfc, color = category)) +
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
      title = paste0("logFC Concordance: Cells vs ", contrast),
      subtitle = paste0("r = ", round(cor_all, 3), " (all genes); r = ", round(cor_sig, 3), " (DEGs)"),
      x = expression(log[2]~FC~"(3T3-L1 Cell Culture)"),
      y = expression(log[2]~FC~"("*Tissue*")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)

  scatter_file <- paste0("Concordance_scatter_", contrast, "_", run_tag, ".png")
  ggsave(file.path(concordance_dir, scatter_file), plot = p_scatter,
         width = concordance_params$plot_width, height = concordance_params$plot_height,
         dpi = 300, bg = "white")
  cat("  [OK] Saved: ", scatter_file, "\n", sep = "")

  ## --- Venn diagram ---
  venn_file <- file.path(concordance_dir, paste0("Venn_DEGs_cells_vs_", contrast, "_", run_tag, ".png"))
  png(venn_file, width = 8, height = 6, units = "in", res = 300)
  venn.plot <- draw.pairwise.venn(
    area1 = length(cell_deg_all),
    area2 = length(overlap_results[[contrast]]$tissue_deg_all),
    cross.area = length(overlap_results[[contrast]]$shared_all),
    category = c("Cell Culture", contrast),
    fill = c("#FF9800", "#2196F3"),
    alpha = 0.5,
    cat.cex = 1.2,
    cex = 1.5,
    fontfamily = "sans",
    cat.fontfamily = "sans"
  )
  grid.draw(venn.plot)
  dev.off()
  cat("  [OK] Saved Venn diagram\n")
}

## ============================================================
## 6) ANALYSIS 3: CELL-AUTONOMOUS SIGNATURE
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 3: CELL-AUTONOMOUS SIGNATURE ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  Genes DEG in BOTH cell culture AND tissue are likely cell-autonomous\n")
cat("  (direct KAT8 targets, not secondary microenvironment effects)\n\n")

## For each tissue contrast, identify cell-autonomous genes
for (contrast in names(tissue_tables)) {
  tt <- tissue_tables[[contrast]]
  res <- overlap_results[[contrast]]

  cat("--- ", contrast, " ---\n", sep = "")

  if (length(res$shared_concordant) == 0) {
    cat("  No concordant DEGs found\n\n")
    next
  }

  ## Build cell-autonomous gene table
  ca_genes <- res$shared_concordant
  ca_df <- data.frame(
    gene = ca_genes,
    cell_lfc = cell_de[ca_genes, "log2FoldChange"],
    cell_padj = cell_de[ca_genes, "padj"],
    tissue_lfc = tt[ca_genes, "log2FoldChange"],
    tissue_padj = tt[ca_genes, "padj"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      direction = ifelse(cell_lfc > 0, "Up", "Down"),
      mean_lfc = (cell_lfc + tissue_lfc) / 2
    ) %>%
    arrange(desc(abs(mean_lfc)))

  cat("  Cell-autonomous genes: ", nrow(ca_df), "\n", sep = "")
  cat("    Up in both:   ", sum(ca_df$direction == "Up"), "\n", sep = "")
  cat("    Down in both: ", sum(ca_df$direction == "Down"), "\n", sep = "")

  ## Show top hits
  cat("  Top cell-autonomous genes (by mean |logFC|):\n")
  for (i in seq_len(min(10, nrow(ca_df)))) {
    row <- ca_df[i, ]
    cat("    ", row$gene, ": cell=", round(row$cell_lfc, 2),
        " tissue=", round(row$tissue_lfc, 2), " (", row$direction, ")\n", sep = "")
  }

  ## Save
  ca_file <- paste0("cell_autonomous_genes_", contrast, "_", run_tag, ".csv")
  write.csv(ca_df, file = file.path(tables_dir, ca_file), row.names = FALSE)
  cat("  [OK] Saved: ", ca_file, "\n\n", sep = "")
}

## ============================================================
## 7) ANALYSIS 4: MICROENVIRONMENT-ONLY GENES
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 4: MICROENVIRONMENT-DRIVEN GENES ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  Genes DEG in tissue but NOT in cells likely reflect:\n")
cat("  - Non-adipocyte cell type changes (immune cells, endothelial, etc.)\n")
cat("  - Microenvironment-dependent adipocyte responses\n")
cat("  - Paracrine/systemic effects absent in culture\n\n")

for (contrast in names(tissue_tables)) {
  res <- overlap_results[[contrast]]
  tt <- tissue_tables[[contrast]]

  tissue_only_genes <- res$tissue_only
  cat("--- ", contrast, " ---\n", sep = "")
  cat("  Tissue-only DEGs: ", length(tissue_only_genes), "\n", sep = "")

  if (length(tissue_only_genes) == 0) next

  ## Get details
  tissue_only_df <- data.frame(
    gene = tissue_only_genes,
    tissue_lfc = tt[tissue_only_genes, "log2FoldChange"],
    tissue_padj = tt[tissue_only_genes, "padj"],
    stringsAsFactors = FALSE
  )

  ## Add cell culture data if gene exists there
  common_with_cells <- intersect(tissue_only_genes, rownames(cell_de))
  if (length(common_with_cells) > 0) {
    tissue_only_df$cell_lfc <- NA
    tissue_only_df$cell_padj <- NA
    tissue_only_df[tissue_only_df$gene %in% common_with_cells, "cell_lfc"] <-
      cell_de[common_with_cells, "log2FoldChange"]
    tissue_only_df[tissue_only_df$gene %in% common_with_cells, "cell_padj"] <-
      cell_de[common_with_cells, "padj"]
  }

  tissue_only_df <- tissue_only_df %>% arrange(desc(abs(tissue_lfc)))

  env_file <- paste0("microenvironment_genes_", contrast, "_", run_tag, ".csv")
  write.csv(tissue_only_df, file = file.path(tables_dir, env_file), row.names = FALSE)
  cat("  [OK] Saved: ", env_file, "\n\n", sep = "")
}

## ============================================================
## 8) SUMMARY TABLE
## ============================================================

cat("\n=== GENERATING SUMMARY ===\n")

summary_rows <- list()
for (contrast in names(overlap_results)) {
  res <- overlap_results[[contrast]]
  summary_rows[[contrast]] <- data.frame(
    Tissue_Contrast = contrast,
    N_Cell_DEGs = length(cell_deg_all),
    N_Tissue_DEGs = length(res$tissue_deg_all),
    N_Shared = length(res$shared_all),
    N_Concordant = length(res$shared_concordant),
    N_Discordant = length(res$shared_discordant),
    N_Cell_Only = length(res$cell_only),
    N_Tissue_Only = length(res$tissue_only),
    Pct_Cell_Shared = round(100 * length(res$shared_all) / max(1, length(cell_deg_all)), 1),
    Pct_Tissue_Shared = round(100 * length(res$shared_all) / max(1, length(res$tissue_deg_all)), 1),
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
summary_file <- paste0("concordance_summary_", run_tag, ".csv")
write.csv(summary_df, file = file.path(tables_dir, summary_file), row.names = FALSE)
cat("[OK] Saved summary: ", summary_file, "\n", sep = "")

## Print summary
cat("\n")
cat("============================================================\n")
cat("=== CONCORDANCE SUMMARY ===\n")
cat("============================================================\n")
print(summary_df)
cat("\n")

## ============================================================
## 9) SESSION INFO
## ============================================================

session_file <- file.path(outdir, "logs", paste0("sessionInfo_concordance_", run_tag, ".txt"))
if (dir.exists(dirname(session_file))) {
  sink(session_file)
  cat("=== R SESSION INFORMATION - Part 7 (In Vitro / In Vivo Concordance) ===\n\n")
  print(sessionInfo())
  sink()
  cat("[INFO] Session info saved\n")
}

cat("\n=== Part 7 (In Vitro / In Vivo Concordance) COMPLETE ===\n")
cat("Output files in: ", concordance_dir, "\n", sep = "")
