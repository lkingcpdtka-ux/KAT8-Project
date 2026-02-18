#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 7.5: SEX CONCORDANCE WITHIN DEPOTS
## =========================================================
## Compares Male vs Female KAT8-KD responses WITHIN each depot
## to justify merging sexes in the main analysis.
##
## KEY QUESTION FOR REVIEWERS:
##   "Do males and females show the same transcriptional response
##    to KAT8 knockdown within each depot?"
##
## KEY ANALYSES (per depot):
## 1. logFC Correlation: M vs F scatter (high r = justify merge)
## 2. DEG Overlap: Venn diagrams (M vs F DEGs)
## 3. Concordance Classification: Same-direction, opposite, sex-specific
## 4. Summary table with stats across both depots
##
## SELF-CONTAINED: Computes per-sex DESeq2 from counts.txt
## (no dependency on Part 1 sex-specific DE tables)
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "scales",
                   "RColorBrewer", "VennDiagram", "grid")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c("DESeq2")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
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

## Prefilter defaults
if (!exists("prefilter_min_count"))   prefilter_min_count   <- 10
if (!exists("prefilter_min_samples")) prefilter_min_samples <- 4

## Helper: safe DESeq2 results to data.frame
deseq_res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  std_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  if (ncol(df) == length(std_cols)) colnames(df) <- std_cols
  df
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

## ============================================================
## 3) LOAD RAW COUNTS + SAMPLE ANNOTATION
## ============================================================
cat("\n=== LOADING DATA FOR SEX CONCORDANCE ===\n")

counts_file <- file.path(getwd(), "counts.txt")
if (!file.exists(counts_file)) stop("counts.txt not found")
counts_raw <- read.delim(counts_file, header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)

## Sample annotation (same as part1)
sample_ids <- paste0("JS_", sprintf("%02d", 1:40))
sample_annot <- data.frame(
  Sample   = sample_ids,
  Depot    = c(rep("iWAT", 10), rep("iWAT", 10),
               rep("gWAT", 10), rep("gWAT", 10)),
  Sex      = c(rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5),
               rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5)),
  Genotype = c(rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5)),
  stringsAsFactors = FALSE
)

## Remove outliers
outlier_samples <- c("JS_08", "JS_28")
sample_annot <- sample_annot[!sample_annot$Sample %in% outlier_samples, ]
tissue_samples <- intersect(sample_annot$Sample, colnames(counts_raw))
sample_annot <- sample_annot[sample_annot$Sample %in% tissue_samples, ]

## Factorize
sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
sample_annot$Depot    <- factor(sample_annot$Depot, levels = c("iWAT", "gWAT"))
sample_annot$Sex      <- factor(sample_annot$Sex, levels = c("F", "M"))
rownames(sample_annot) <- sample_annot$Sample

## Build count matrix
if ("gene" %in% colnames(counts_raw) || "gene_name" %in% colnames(counts_raw)) {
  counts_tissue <- counts_raw
  if (!"gene_name" %in% colnames(counts_tissue) && "gene" %in% colnames(counts_tissue))
    counts_tissue$gene_name <- counts_tissue$gene
  if (!"ensembl_gene_id" %in% colnames(counts_tissue) && "gene_id" %in% colnames(counts_tissue))
    counts_tissue$ensembl_gene_id <- counts_tissue$gene_id
  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  if ("ensembl_gene_id" %in% colnames(counts_tissue))
    counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]
  if (any(duplicated(counts_tissue$gene_name))) {
    dup_flag <- duplicated(counts_tissue$gene_name) | duplicated(counts_tissue$gene_name, fromLast = TRUE)
    if ("ensembl_gene_id" %in% colnames(counts_tissue))
      counts_tissue$gene_name[dup_flag] <- paste0(counts_tissue$ensembl_gene_id[dup_flag], "_",
                                                   counts_tissue$gene_name[dup_flag])
    counts_tissue$gene_name <- make.unique(counts_tissue$gene_name)
  }
  rownames(counts_tissue) <- counts_tissue$gene_name
  count_matrix <- as.matrix(counts_tissue[, tissue_samples])
  rownames(count_matrix) <- counts_tissue$gene_name
} else {
  count_matrix <- as.matrix(counts_raw[, tissue_samples])
}

## Prefilter
keep <- rowSums(count_matrix >= prefilter_min_count) >= prefilter_min_samples
count_matrix <- count_matrix[keep, , drop = FALSE]
cat("[INFO] Genes after prefiltering: ", nrow(count_matrix), "\n")

## ============================================================
## 3b) OUTPUT SETUP
## ============================================================
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) > 0) {
  outdir <- run_dirs[which.max(file.info(run_dirs)$mtime)]
} else {
  outdir <- file.path(savepoint_dir,
                      paste0("RUN_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
sex_conc_dir <- file.path(plots_dir, "sex_concordance")
if (!dir.exists(sex_conc_dir)) dir.create(sex_conc_dir, recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

## Colors
col_male   <- "#2166AC"   ## dark blue
col_female <- "#B2182B"   ## dark red

## Category colors for scatter
category_colors <- c(
  "Concordant"  = "#4CAF50",
  "Discordant"  = "#F44336",
  "Male only"   = col_male,
  "Female only" = col_female,
  "NS"          = "grey85"
)

## ============================================================
## 4) ANALYSIS: SEX CONCORDANCE WITHIN EACH DEPOT
## ============================================================

## Collect summary rows for the combined table
all_summary_rows <- list()

for (depot in c("iWAT", "gWAT")) {
  cat("\n")
  cat("============================================================\n")
  cat("=== SEX CONCORDANCE: ", depot, " (Male vs Female) ===\n", sep = "")
  cat("============================================================\n\n")

  ## ---------------------------------------------------------
  ## 4.0) Subset to depot and run per-sex DESeq2
  ## ---------------------------------------------------------
  depot_samples <- sample_annot$Sample[sample_annot$Depot == depot]
  depot_annot   <- sample_annot[depot_samples, ]
  depot_counts  <- count_matrix[, depot_samples]

  ## Re-prefilter within depot
  depot_keep <- rowSums(depot_counts >= prefilter_min_count) >= prefilter_min_samples
  depot_counts <- depot_counts[depot_keep, , drop = FALSE]
  depot_gene_names <- rownames(depot_counts)

  cat("[INFO] ", depot, ": ", nrow(depot_counts), " genes, ",
      ncol(depot_counts), " samples\n", sep = "")

  ## Female-only DESeq2
  female_idx <- depot_annot$Sample[depot_annot$Sex == "F"]
  cat("[INFO] Running female-only DESeq2 (n = ", length(female_idx), ")...\n", sep = "")
  dds_f <- DESeqDataSetFromMatrix(
    countData = round(depot_counts[, female_idx]),
    colData   = depot_annot[female_idx, ],
    design    = ~ Genotype
  )
  dds_f <- DESeq(dds_f, quiet = TRUE)
  de_f <- deseq_res_to_df(results(dds_f))
  de_f$gene <- depot_gene_names
  rownames(de_f) <- depot_gene_names

  ## Male-only DESeq2
  male_idx <- depot_annot$Sample[depot_annot$Sex == "M"]
  cat("[INFO] Running male-only DESeq2 (n = ", length(male_idx), ")...\n", sep = "")
  dds_m <- DESeqDataSetFromMatrix(
    countData = round(depot_counts[, male_idx]),
    colData   = depot_annot[male_idx, ],
    design    = ~ Genotype
  )
  dds_m <- DESeq(dds_m, quiet = TRUE)
  de_m <- deseq_res_to_df(results(dds_m))
  de_m$gene <- depot_gene_names
  rownames(de_m) <- depot_gene_names

  ## Detect logFC column
  lfc_col_m <- "log2FoldChange"
  lfc_col_f <- "log2FoldChange"

  ## Get thresholds (same for M and F within a depot)
  thresh <- get_tissue_thresholds(depot)

  ## -----------------------------------------------------------
  ## 4a) DEG sets
  ## -----------------------------------------------------------
  m_up <- rownames(de_m)[!is.na(de_m$padj) & de_m$padj < thresh$fdr_cut &
                           de_m[[lfc_col_m]] > thresh$logFC_cut]
  m_down <- rownames(de_m)[!is.na(de_m$padj) & de_m$padj < thresh$fdr_cut &
                             de_m[[lfc_col_m]] < -thresh$logFC_cut]
  m_all <- union(m_up, m_down)

  f_up <- rownames(de_f)[!is.na(de_f$padj) & de_f$padj < thresh$fdr_cut &
                           de_f[[lfc_col_f]] > thresh$logFC_cut]
  f_down <- rownames(de_f)[!is.na(de_f$padj) & de_f$padj < thresh$fdr_cut &
                             de_f[[lfc_col_f]] < -thresh$logFC_cut]
  f_all <- union(f_up, f_down)

  cat("Male DEGs:   ", length(m_all), " (", length(m_up), " up, ", length(m_down), " down)\n", sep = "")
  cat("Female DEGs: ", length(f_all), " (", length(f_up), " up, ", length(f_down), " down)\n", sep = "")

  ## -----------------------------------------------------------
  ## 4b) Overlap
  ## -----------------------------------------------------------
  shared_all   <- intersect(m_all, f_all)
  m_only       <- setdiff(m_all, f_all)
  f_only       <- setdiff(f_all, m_all)

  shared_up_up     <- intersect(m_up, f_up)
  shared_down_down <- intersect(m_down, f_down)
  shared_concordant <- union(shared_up_up, shared_down_down)

  shared_discordant <- union(intersect(m_up, f_down), intersect(m_down, f_up))

  concordance_pct <- if (length(shared_all) > 0) {
    round(100 * length(shared_concordant) / length(shared_all), 1)
  } else {
    NA_real_
  }

  cat("\nShared DEGs: ", length(shared_all), "\n", sep = "")
  cat("  Concordant (same direction): ", length(shared_concordant), "\n", sep = "")
  cat("    Both up:   ", length(shared_up_up), "\n", sep = "")
  cat("    Both down: ", length(shared_down_down), "\n", sep = "")
  cat("  Discordant (opposite): ", length(shared_discordant), "\n", sep = "")
  cat("  Male-only:   ", length(m_only), "\n", sep = "")
  cat("  Female-only: ", length(f_only), "\n", sep = "")
  cat("  Concordance rate: ", concordance_pct, "%\n", sep = "")

  ## Fisher's exact test
  common_genes <- intersect(rownames(de_m), rownames(de_f))
  a  <- length(shared_all)
  b  <- length(setdiff(f_all, m_all))
  cc <- length(setdiff(m_all, f_all))
  d  <- length(common_genes) - a - b - cc
  fisher_res <- fisher.test(matrix(c(a, b, cc, d), nrow = 2))
  cat("  Fisher's exact test: OR = ", round(fisher_res$estimate, 2),
      ", p = ", format(fisher_res$p.value, digits = 3), "\n", sep = "")

  ## -----------------------------------------------------------
  ## 4c) logFC correlation
  ## -----------------------------------------------------------
  cor_df <- data.frame(
    gene     = common_genes,
    m_lfc    = de_m[common_genes, lfc_col_m],
    f_lfc    = de_f[common_genes, lfc_col_f],
    m_padj   = de_m[common_genes, "padj"],
    f_padj   = de_f[common_genes, "padj"],
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(m_lfc), is.finite(f_lfc))

  ## Classify
  cor_df <- cor_df %>%
    mutate(
      m_sig = !is.na(m_padj) & m_padj < thresh$fdr_cut & abs(m_lfc) > thresh$logFC_cut,
      f_sig = !is.na(f_padj) & f_padj < thresh$fdr_cut & abs(f_lfc) > thresh$logFC_cut,
      category = case_when(
        m_sig & f_sig & sign(m_lfc) == sign(f_lfc) ~ "Concordant",
        m_sig & f_sig & sign(m_lfc) != sign(f_lfc) ~ "Discordant",
        m_sig & !f_sig ~ "Male only",
        !m_sig & f_sig ~ "Female only",
        TRUE ~ "NS"
      )
    )

  ## Correlations
  cor_pearson  <- cor(cor_df$m_lfc, cor_df$f_lfc, method = "pearson")
  cor_spearman <- cor(cor_df$m_lfc, cor_df$f_lfc, method = "spearman")
  cor_test     <- cor.test(cor_df$m_lfc, cor_df$f_lfc, method = "pearson")

  sig_df <- cor_df %>% filter(m_sig | f_sig)
  cor_sig <- if (nrow(sig_df) > 10) {
    cor(sig_df$m_lfc, sig_df$f_lfc, method = "pearson")
  } else {
    NA
  }

  cat("\nlogFC Correlation:\n")
  cat("  Pearson r (all genes):  ", round(cor_pearson, 3),
      " (p = ", format(cor_test$p.value, digits = 3), ")\n", sep = "")
  cat("  Pearson r (DEGs only):  ", round(cor_sig, 3), "\n", sep = "")
  cat("  Spearman rho:           ", round(cor_spearman, 3), "\n", sep = "")

  ## -----------------------------------------------------------
  ## 4d) Scatter plot (M vs F)
  ## -----------------------------------------------------------
  top_n_label <- depot_concordance_params$top_n_label

  label_df <- cor_df %>%
    filter(category %in% c("Concordant", "Discordant")) %>%
    arrange(desc(abs(m_lfc) + abs(f_lfc))) %>%
    slice_head(n = top_n_label)

  ## Axis limits: 99th pct but guarantee labeled genes fit
  pct_max <- max(
    quantile(abs(cor_df$m_lfc), 0.99, na.rm = TRUE),
    quantile(abs(cor_df$f_lfc), 0.99, na.rm = TRUE)
  )
  label_max <- if (nrow(label_df) > 0) {
    max(abs(label_df$m_lfc), abs(label_df$f_lfc), na.rm = TRUE)
  } else {
    pct_max
  }
  axis_max <- max(pct_max, label_max) * 1.10

  p_scatter <- ggplot(cor_df, aes(x = m_lfc, y = f_lfc, color = category)) +
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
      title = paste0(depot, ": Male vs Female KAT8-KD Response"),
      subtitle = paste0(
        "Pearson r = ", round(cor_pearson, 3),
        " (all); r = ", round(cor_sig, 3),
        " (DEGs); concordance = ", concordance_pct, "%"
      ),
      x = expression(log[2]~FC~"(Male KD vs CTL)"),
      y = expression(log[2]~FC~"(Female KD vs CTL)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(-axis_max, axis_max), ylim = c(-axis_max, axis_max))

  scatter_file <- paste0("Sex_concordance_scatter_", depot, "_", run_tag, ".png")
  ggsave(file.path(sex_conc_dir, scatter_file), plot = p_scatter,
         width = depot_concordance_params$plot_width,
         height = depot_concordance_params$plot_height,
         dpi = 300, bg = "white")
  cat("[OK] Saved: ", scatter_file, "\n", sep = "")

  ## -----------------------------------------------------------
  ## 4e) Venn diagrams (all, UP, DOWN)
  ## -----------------------------------------------------------
  ## All DEGs
  venn_file <- file.path(sex_conc_dir, paste0("Venn_sex_all_", depot, "_", run_tag, ".png"))
  png(venn_file, width = 2400, height = 2400, res = 300)
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
  v <- venn.diagram(
    x = list(Male = m_all, Female = f_all),
    category.names = c("Male", "Female"),
    filename = NULL, output = TRUE,
    fill = c(col_male, col_female), alpha = 0.5,
    cex = 1.5, cat.cex = 1.5, cat.fontface = "bold",
    cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05),
    main = paste0(depot, ": All DEGs (Male vs Female)"),
    main.cex = 1.5, main.pos = c(0.5, 1.05)
  )
  grid.draw(v)
  popViewport()
  dev.off()

  ## UP DEGs
  venn_file <- file.path(sex_conc_dir, paste0("Venn_sex_UP_", depot, "_", run_tag, ".png"))
  png(venn_file, width = 2400, height = 2400, res = 300)
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
  v <- venn.diagram(
    x = list(Male_UP = m_up, Female_UP = f_up),
    category.names = c("Male UP", "Female UP"),
    filename = NULL, output = TRUE,
    fill = c(col_male, col_female), alpha = 0.5,
    cex = 1.5, cat.cex = 1.5, cat.fontface = "bold",
    cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05),
    main = paste0(depot, ": UP-regulated DEGs"),
    main.cex = 1.5, main.pos = c(0.5, 1.05)
  )
  grid.draw(v)
  popViewport()
  dev.off()

  ## DOWN DEGs
  venn_file <- file.path(sex_conc_dir, paste0("Venn_sex_DOWN_", depot, "_", run_tag, ".png"))
  png(venn_file, width = 2400, height = 2400, res = 300)
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9))
  v <- venn.diagram(
    x = list(Male_DOWN = m_down, Female_DOWN = f_down),
    category.names = c("Male DOWN", "Female DOWN"),
    filename = NULL, output = TRUE,
    fill = c(col_male, col_female), alpha = 0.5,
    cex = 1.5, cat.cex = 1.5, cat.fontface = "bold",
    cat.pos = c(-20, 20), cat.dist = c(0.05, 0.05),
    main = paste0(depot, ": DOWN-regulated DEGs"),
    main.cex = 1.5, main.pos = c(0.5, 1.05)
  )
  grid.draw(v)
  popViewport()
  dev.off()
  cat("[OK] Saved 3 Venn diagrams for ", depot, "\n", sep = "")

  ## -----------------------------------------------------------
  ## 4f) Save gene lists
  ## -----------------------------------------------------------
  ## Concordant genes (shared M/F response)
  if (length(shared_concordant) > 0) {
    conc_df <- data.frame(
      gene     = shared_concordant,
      m_lfc    = de_m[shared_concordant, lfc_col_m],
      m_padj   = de_m[shared_concordant, "padj"],
      f_lfc    = de_f[shared_concordant, lfc_col_f],
      f_padj   = de_f[shared_concordant, "padj"],
      stringsAsFactors = FALSE
    ) %>%
      mutate(direction = ifelse(m_lfc > 0, "Up", "Down"),
             mean_lfc = (m_lfc + f_lfc) / 2) %>%
      arrange(desc(abs(mean_lfc)))

    write.csv(conc_df,
              file.path(tables_dir, paste0("sex_concordant_genes_", depot, "_", run_tag, ".csv")),
              row.names = FALSE)
    cat("[OK] Saved concordant gene list for ", depot, "\n", sep = "")
  }

  ## Sex-specific genes
  if (length(m_only) > 0) {
    m_only_df <- data.frame(
      gene   = m_only,
      m_lfc  = de_m[m_only, lfc_col_m],
      m_padj = de_m[m_only, "padj"],
      stringsAsFactors = FALSE
    )
    ## Cross-reference: what does the female analysis show?
    common_m_only <- intersect(m_only, rownames(de_f))
    m_only_df$f_lfc  <- NA_real_
    m_only_df$f_padj <- NA_real_
    if (length(common_m_only) > 0) {
      m_only_df[m_only_df$gene %in% common_m_only, "f_lfc"]  <- de_f[common_m_only, lfc_col_f]
      m_only_df[m_only_df$gene %in% common_m_only, "f_padj"] <- de_f[common_m_only, "padj"]
    }
    write.csv(m_only_df %>% arrange(desc(abs(m_lfc))),
              file.path(tables_dir, paste0("sex_male_only_genes_", depot, "_", run_tag, ".csv")),
              row.names = FALSE)
  }

  if (length(f_only) > 0) {
    f_only_df <- data.frame(
      gene   = f_only,
      f_lfc  = de_f[f_only, lfc_col_f],
      f_padj = de_f[f_only, "padj"],
      stringsAsFactors = FALSE
    )
    common_f_only <- intersect(f_only, rownames(de_m))
    f_only_df$m_lfc  <- NA_real_
    f_only_df$m_padj <- NA_real_
    if (length(common_f_only) > 0) {
      f_only_df[f_only_df$gene %in% common_f_only, "m_lfc"]  <- de_m[common_f_only, lfc_col_m]
      f_only_df[f_only_df$gene %in% common_f_only, "m_padj"] <- de_m[common_f_only, "padj"]
    }
    write.csv(f_only_df %>% arrange(desc(abs(f_lfc))),
              file.path(tables_dir, paste0("sex_female_only_genes_", depot, "_", run_tag, ".csv")),
              row.names = FALSE)
  }

  ## -----------------------------------------------------------
  ## 4g) Collect summary row
  ## -----------------------------------------------------------
  all_summary_rows[[depot]] <- data.frame(
    Depot = depot,
    Male_DEGs = length(m_all),
    Male_UP = length(m_up),
    Male_DOWN = length(m_down),
    Female_DEGs = length(f_all),
    Female_UP = length(f_up),
    Female_DOWN = length(f_down),
    Shared = length(shared_all),
    Concordant = length(shared_concordant),
    Discordant = length(shared_discordant),
    Male_only = length(m_only),
    Female_only = length(f_only),
    Concordance_pct = concordance_pct,
    Pct_Male_shared = round(100 * length(shared_all) / max(1, length(m_all)), 1),
    Pct_Female_shared = round(100 * length(shared_all) / max(1, length(f_all)), 1),
    Pearson_r_all = round(cor_pearson, 4),
    Pearson_r_DEGs = round(cor_sig, 4),
    Spearman_rho = round(cor_spearman, 4),
    Fisher_OR = round(fisher_res$estimate, 2),
    Fisher_p = format(fisher_res$p.value, digits = 4),
    stringsAsFactors = FALSE
  )
}

## ============================================================
## 5) COMBINED SUMMARY TABLE
## ============================================================
cat("\n=== GENERATING COMBINED SUMMARY ===\n")

summary_df <- bind_rows(all_summary_rows)

summary_file <- paste0("sex_concordance_summary_", run_tag, ".csv")
write.csv(summary_df, file = file.path(tables_dir, summary_file), row.names = FALSE)
cat("[OK] Saved: ", summary_file, "\n", sep = "")

## Print summary
cat("\n")
cat("============================================================\n")
cat("=== SEX CONCORDANCE SUMMARY ===\n")
cat("============================================================\n")
cat("INTERPRETATION FOR REVIEWERS:\n")
cat("  High Pearson r + high concordance % = sexes respond similarly\n")
cat("  This justifies pooling males and females within each depot.\n\n")
print(summary_df, row.names = FALSE)
cat("\n")

## ============================================================
## 6) SESSION INFO
## ============================================================
log_dir <- file.path(outdir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
session_file <- file.path(log_dir, paste0("sessionInfo_sex_concordance_", run_tag, ".txt"))
sink(session_file)
cat("=== R SESSION INFORMATION - Part 7.5 (Sex Concordance) ===\n\n")
print(sessionInfo())
sink()
cat("[INFO] Session info saved\n")

cat("\n=== Part 7.5 (Sex Concordance) COMPLETE ===\n")
cat("Output directory: ", sex_conc_dir, "\n", sep = "")
cat("\nFiles created (per depot):\n")
cat("  - logFC scatter plot (Male vs Female)\n")
cat("  - 3 Venn diagrams (all, UP, DOWN)\n")
cat("  - Concordant gene list\n")
cat("  - Male-only / Female-only gene lists (with cross-reference)\n")
cat("  - Combined summary table\n\n")
