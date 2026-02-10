#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 9: UPSTREAM REGULATOR / TF ANALYSIS
## =========================================================
## Identifies transcription factors whose targets are enriched
## among KAT8-regulated genes, connecting expression changes
## to the upstream epigenetic mechanism
##
## KEY ANALYSES:
## 1. MSigDB C3 TFT: TF target enrichment (motif-based)
## 2. TRRUST-style: Literature-curated TF-target relationships
## 3. Chromatin Modifier Focus: Highlight HATs/HDACs/chromatin TFs
## 4. Cross-Contrast Comparison: Shared upstream regulators
## 5. KAT8-Centric Integration: Connect to known KAT8/MOF biology
##
## REQUIRES:
##   - Part 1 tissue results (savepoints/RUN_*/tables/DE_tissue_*.csv)
##   - msigdbr package for gene set databases
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "tibble", "stringr", "msigdbr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "fgsea"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(stringr)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(fgsea)
})

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

set.seed(MASTER_SEED)

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
## 3) KNOWN CHROMATIN MODIFIERS / EPIGENETIC REGULATORS
## ============================================================
## Curated list for highlighting in results
## KAT8/MOF is part of the NSL and MSL complexes

chromatin_modifiers <- list(
  ## Histone acetyltransferases (HATs)
  HATs = c(
    "KAT8", "MOF",           ## The gene of interest (NSL/MSL complex)
    "KAT2A", "GCN5",         ## SAGA complex
    "KAT2B", "PCAF",         ## SAGA/ATAC
    "EP300", "P300",          ## CBP/p300
    "CREBBP", "CBP",          ## CBP
    "KAT5", "TIP60",         ## NuA4/TIP60
    "KAT6A", "MOZ",          ## MOZ/MORF
    "KAT6B", "MORF",
    "KAT7", "HBO1",          ## HBO1
    "TAF1"                    ## TFIID
  ),

  ## Histone deacetylases (HDACs)
  HDACs = c(
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5",
    "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10",
    "HDAC11", "SIRT1", "SIRT2", "SIRT3", "SIRT6", "SIRT7"
  ),

  ## KAT8/MOF complex members
  NSL_MSL_complex = c(
    "KAT8", "KANSL1", "KANSL2", "KANSL3",   ## NSL complex
    "MCRS1", "WDR5", "PHF20", "OGT", "HCF1",
    "MSL1", "MSL2", "MSL3"                    ## MSL complex
  ),

  ## Chromatin remodelers
  Remodelers = c(
    "SMARCA4", "BRG1", "SMARCA2", "BRM",    ## SWI/SNF
    "ARID1A", "ARID1B", "ARID2",
    "CHD1", "CHD2", "CHD4", "CHD7",          ## CHD family
    "EP400", "SRCAP"                           ## INO80/SRCAP
  ),

  ## Histone methyltransferases
  HMTs = c(
    "EZH2", "EZH1", "SUV39H1", "SUV39H2",
    "SETD1A", "SETD1B", "KMT2A", "MLL1",
    "DOT1L", "NSD1", "NSD2", "NSD3"
  )
)

## Flatten to a lookup vector
chromatin_all <- unique(unlist(chromatin_modifiers))
chromatin_category <- data.frame(
  TF = rep(names(chromatin_modifiers), sapply(chromatin_modifiers, length)),
  Gene = unlist(chromatin_modifiers),
  stringsAsFactors = FALSE
) %>% distinct(Gene, .keep_all = TRUE)

## ============================================================
## 4) LOAD TF TARGET GENE SETS FROM MSIGDB
## ============================================================

cat("\n=== LOADING TF TARGET GENE SETS ===\n")

## C3: Regulatory target gene sets (TFT = transcription factor targets)
if (tf_params$use_c3_tft) {
  cat("[INFO] Loading MSigDB C3:TFT (TF targets, motif-based)...\n")

  ## C3 TFT: TF binding motifs
  tryCatch({
    c3_tft <- msigdbr(species = "Mus musculus", collection = "C3", subcollection = "TFT:GTRD")
    if (nrow(c3_tft) == 0) {
      ## Fallback to broader C3 TFT
      c3_tft <- msigdbr(species = "Mus musculus", collection = "C3", subcollection = "TFT:TFT_Legacy")
    }
    if (nrow(c3_tft) == 0) {
      c3_tft <- msigdbr(species = "Mus musculus", collection = "C3")
    }

    tft_list <- split(c3_tft$gene_symbol, c3_tft$gs_name)
    tft_list <- tft_list[sapply(tft_list, length) >= tf_params$min_targets &
                           sapply(tft_list, length) <= 500]
    cat("[OK] TF target gene sets: ", length(tft_list), "\n", sep = "")
  }, error = function(e) {
    cat("[WARN] Failed to load C3:TFT: ", conditionMessage(e), "\n", sep = "")
    tft_list <<- NULL
  })
} else {
  tft_list <- NULL
}

## Also load Hallmark for context
tryCatch({
  hallmark_sets <- msigdbr(species = "Mus musculus", collection = "H")
  hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
  cat("[OK] Hallmark gene sets: ", length(hallmark_list), "\n", sep = "")
}, error = function(e) {
  hallmark_list <- NULL
})

## ============================================================
## 5) FIND PART 1 RESULTS
## ============================================================

cat("\n=== FINDING PART 1 RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run Part 1 first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")
tf_dir <- file.path(plots_dir, "tf_analysis")
if (!dir.exists(tf_dir)) dir.create(tf_dir, recursive = TRUE)

## Load DE tables
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
de_tables <- list()
for (f in de_files) {
  contrast <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(f))
  de_tables[[contrast]] <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  cat("[OK] Loaded: ", contrast, " (", nrow(de_tables[[contrast]]), " genes)\n", sep = "")
}

## ============================================================
## 6) ANALYSIS 1: TF TARGET ENRICHMENT (GSEA-based)
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 1: TF TARGET ENRICHMENT (fGSEA) ===\n")
cat("============================================================\n\n")
cat("RATIONALE:\n")
cat("  KAT8 is a histone acetyltransferase (H4K16ac).\n")
cat("  Loss of KAT8 -> loss of H4K16ac -> altered TF binding/activity.\n")
cat("  This analysis identifies TFs whose targets are enriched among DEGs.\n\n")

tf_enrichment_results <- list()

for (contrast in names(de_tables)) {
  cat("--- ", contrast, " ---\n", sep = "")

  tt <- de_tables[[contrast]]

  ## Build ranked gene list (Wald statistic)
  if (!"stat" %in% colnames(tt)) {
    cat("  [WARN] No 'stat' column found. Skipping.\n\n")
    next
  }

  ranked_df <- tt %>%
    filter(is.finite(stat)) %>%
    arrange(desc(stat))
  ranked_vec <- setNames(ranked_df$stat, rownames(ranked_df))

  ## Run fGSEA against TF target sets
  if (!is.null(tft_list) && length(tft_list) > 0) {
    tryCatch({
      fgsea_tf <- fgsea(
        pathways = tft_list,
        stats = ranked_vec,
        minSize = tf_params$min_targets,
        maxSize = 500,
        nPermSimple = 10000
      )

      ## Filter significant
      sig_tf <- fgsea_tf %>%
        filter(padj < tf_params$pvalue_cutoff) %>%
        arrange(padj)

      cat("  Significant TF target sets: ", nrow(sig_tf), "\n", sep = "")

      if (nrow(sig_tf) > 0) {
        ## Extract TF name from gene set name
        sig_tf <- sig_tf %>%
          mutate(
            TF_name = str_extract(pathway, "^[A-Z0-9]+"),
            Direction = ifelse(NES > 0, "Targets Up", "Targets Down"),
            Is_Chromatin_Modifier = toupper(TF_name) %in% toupper(chromatin_all)
          )

        ## Convert list columns for CSV
        for (col_name in colnames(sig_tf)) {
          if (is.list(sig_tf[[col_name]])) {
            sig_tf[[col_name]] <- sapply(sig_tf[[col_name]], function(x) {
              if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ";")
            })
          }
        }

        tf_enrichment_results[[contrast]] <- sig_tf

        ## Save
        tf_file <- paste0("tf_enrichment_", contrast, "_", run_tag, ".csv")
        write.csv(sig_tf, file = file.path(tables_dir, tf_file), row.names = FALSE)
        cat("  [OK] Saved: ", tf_file, "\n", sep = "")

        ## Report top TFs
        cat("  Top 10 enriched TFs:\n")
        top_tfs <- sig_tf %>% slice_head(n = 10)
        for (k in seq_len(nrow(top_tfs))) {
          row <- top_tfs[k, ]
          chrom_tag <- if (row$Is_Chromatin_Modifier) " *CHROMATIN*" else ""
          cat("    ", row$TF_name, ": NES = ", round(row$NES, 2),
              " (padj = ", format(row$padj, digits = 3), ")", chrom_tag, "\n", sep = "")
        }

        ## Flag chromatin modifiers
        chrom_tfs <- sig_tf %>% filter(Is_Chromatin_Modifier)
        if (nrow(chrom_tfs) > 0) {
          cat("\n  Chromatin modifier TFs enriched:\n")
          for (k in seq_len(nrow(chrom_tfs))) {
            cat("    ", chrom_tfs$TF_name[k], " (NES = ", round(chrom_tfs$NES[k], 2), ")\n", sep = "")
          }
        }
      }
    }, error = function(e) {
      cat("  [WARN] fGSEA TF analysis failed: ", conditionMessage(e), "\n", sep = "")
    })
  }
  cat("\n")
}

## ============================================================
## 7) ANALYSIS 2: ORA-BASED TF ENRICHMENT (UP/DOWN SEPARATELY)
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 2: TF TARGET ORA (DIRECTION-SPECIFIC) ===\n")
cat("============================================================\n\n")

tf_ora_results <- list()

for (contrast in names(de_tables)) {
  cat("--- ", contrast, " ---\n", sep = "")

  tt <- de_tables[[contrast]]
  thresholds <- get_tissue_thresholds(contrast)

  ## Get up and down DEGs
  deg_up <- rownames(tt)[!is.na(tt$padj) & tt$padj < thresholds$fdr_cut &
                           tt$log2FoldChange > thresholds$logFC_cut]
  deg_down <- rownames(tt)[!is.na(tt$padj) & tt$padj < thresholds$fdr_cut &
                             tt$log2FoldChange < -thresholds$logFC_cut]

  for (direction in c("Up", "Down")) {
    deg_genes <- if (direction == "Up") deg_up else deg_down
    if (length(deg_genes) < 10) {
      cat("  ", direction, ": Too few DEGs (", length(deg_genes), ")\n", sep = "")
      next
    }

    cat("  ", direction, "-regulated DEGs: ", length(deg_genes), "\n", sep = "")

    ## Test each TF target set for overlap
    if (!is.null(tft_list)) {
      background_genes <- rownames(tt)

      tf_tests <- lapply(names(tft_list), function(tf_set_name) {
        targets <- tft_list[[tf_set_name]]
        targets_in_bg <- intersect(targets, background_genes)
        deg_in_targets <- intersect(deg_genes, targets_in_bg)

        if (length(targets_in_bg) < tf_params$min_targets) return(NULL)

        ## Fisher's exact test
        a <- length(deg_in_targets)
        b <- length(targets_in_bg) - a
        cc <- length(deg_genes) - a
        d <- length(background_genes) - a - b - cc

        ft <- fisher.test(matrix(c(a, b, cc, d), nrow = 2), alternative = "greater")

        data.frame(
          TF_Set = tf_set_name,
          TF_Name = str_extract(tf_set_name, "^[A-Z0-9]+"),
          Direction = direction,
          N_Targets_In_BG = length(targets_in_bg),
          N_DEG_In_Targets = a,
          Odds_Ratio = ft$estimate,
          P_Value = ft$p.value,
          stringsAsFactors = FALSE
        )
      })

      tf_tests_df <- bind_rows(tf_tests)
      if (nrow(tf_tests_df) > 0) {
        tf_tests_df$FDR <- p.adjust(tf_tests_df$P_Value, method = "BH")
        tf_tests_df$Is_Chromatin <- toupper(tf_tests_df$TF_Name) %in% toupper(chromatin_all)

        sig_ora <- tf_tests_df %>%
          filter(FDR < tf_params$qvalue_cutoff) %>%
          arrange(FDR)

        key <- paste0(contrast, "_", direction)
        tf_ora_results[[key]] <- sig_ora

        if (nrow(sig_ora) > 0) {
          cat("    Significant TFs (ORA): ", nrow(sig_ora), "\n", sep = "")

          ## Top 5
          for (k in seq_len(min(5, nrow(sig_ora)))) {
            row <- sig_ora[k, ]
            chrom_tag <- if (row$Is_Chromatin) " *CHROM*" else ""
            cat("      ", row$TF_Name, ": OR=", round(row$Odds_Ratio, 2),
                " (FDR=", format(row$FDR, digits = 3), ")", chrom_tag, "\n", sep = "")
          }
        }
      }
    }
  }
  cat("\n")
}

## Save combined ORA results
if (length(tf_ora_results) > 0) {
  all_ora <- bind_rows(tf_ora_results, .id = "Contrast_Direction")
  ora_file <- paste0("tf_ora_all_contrasts_", run_tag, ".csv")
  write.csv(all_ora, file = file.path(tables_dir, ora_file), row.names = FALSE)
  cat("[OK] Saved TF ORA results: ", ora_file, "\n", sep = "")
}

## ============================================================
## 8) ANALYSIS 3: CROSS-CONTRAST TF COMPARISON
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== ANALYSIS 3: CROSS-CONTRAST TF COMPARISON ===\n")
cat("============================================================\n\n")

if (length(tf_enrichment_results) >= 2) {
  ## Find TFs significant in multiple contrasts
  all_tf_gsea <- bind_rows(tf_enrichment_results, .id = "Contrast")

  tf_counts <- all_tf_gsea %>%
    group_by(TF_name) %>%
    summarise(
      N_Contrasts = n_distinct(Contrast),
      Mean_NES = mean(NES, na.rm = TRUE),
      Contrasts = paste(unique(Contrast), collapse = "; "),
      Is_Chromatin = any(Is_Chromatin_Modifier),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Contrasts), desc(abs(Mean_NES)))

  shared_tfs <- tf_counts %>% filter(N_Contrasts >= 2)

  cat("TFs enriched in 2+ contrasts: ", nrow(shared_tfs), "\n", sep = "")
  if (nrow(shared_tfs) > 0) {
    for (k in seq_len(min(15, nrow(shared_tfs)))) {
      row <- shared_tfs[k, ]
      chrom_tag <- if (row$Is_Chromatin) " *CHROMATIN*" else ""
      cat("  ", row$TF_name, ": ", row$N_Contrasts, " contrasts, mean NES=",
          round(row$Mean_NES, 2), chrom_tag, "\n", sep = "")
    }
  }

  cross_file <- paste0("tf_cross_contrast_", run_tag, ".csv")
  write.csv(tf_counts, file = file.path(tables_dir, cross_file), row.names = FALSE)
  cat("[OK] Saved: ", cross_file, "\n", sep = "")
}

## ============================================================
## 9) VISUALIZATION: TF ENRICHMENT BAR PLOT
## ============================================================

cat("\n=== GENERATING TF ENRICHMENT PLOTS ===\n")

for (contrast in names(tf_enrichment_results)) {
  sig_tf <- tf_enrichment_results[[contrast]]

  if (nrow(sig_tf) == 0) next

  ## Top TFs by absolute NES
  top_up <- sig_tf %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = tf_params$top_n_tfs / 2)
  top_down <- sig_tf %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = tf_params$top_n_tfs / 2)
  plot_data <- bind_rows(top_down, top_up) %>%
    arrange(NES) %>%
    mutate(
      label = ifelse(!is.na(TF_name), TF_name, pathway),
      label = factor(label, levels = label),
      fill_col = case_when(
        Is_Chromatin_Modifier & NES > 0 ~ "Chromatin (Up)",
        Is_Chromatin_Modifier & NES < 0 ~ "Chromatin (Down)",
        NES > 0 ~ "Other TF (Up)",
        TRUE ~ "Other TF (Down)"
      )
    )

  if (nrow(plot_data) == 0) next

  p_tf <- ggplot(plot_data, aes(x = NES, y = label, fill = fill_col)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(
      values = c(
        "Chromatin (Up)" = "#E69F00",
        "Chromatin (Down)" = "#0072B2",
        "Other TF (Up)" = "#FDD49E",
        "Other TF (Down)" = "#9ECAE1"
      ),
      name = "Category"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      title = paste0("Upstream TF Enrichment: ", contrast),
      subtitle = "Chromatin modifiers highlighted (bold colors)",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.text.y = element_text(size = 9)
    )

  tf_plot_file <- paste0("TF_enrichment_", contrast, "_", run_tag, ".png")
  ggsave(file.path(tf_dir, tf_plot_file), plot = p_tf,
         width = tf_params$plot_width, height = max(6, nrow(plot_data) * 0.35),
         dpi = 300, bg = "white")
  cat("[OK] Saved: ", tf_plot_file, "\n", sep = "")
}

## ============================================================
## 10) KAT8-SPECIFIC CONTEXT
## ============================================================

cat("\n")
cat("============================================================\n")
cat("=== KAT8-SPECIFIC BIOLOGICAL CONTEXT ===\n")
cat("============================================================\n\n")

cat("KNOWN KAT8/MOF BIOLOGY:\n")
cat("  - KAT8 (MOF) acetylates H4K16 (primary substrate)\n")
cat("  - Part of NSL complex (broad gene activation) and MSL complex (dosage compensation)\n")
cat("  - H4K16ac is associated with active transcription and open chromatin\n")
cat("  - Loss of KAT8 -> reduced H4K16ac -> chromatin compaction -> gene silencing\n\n")

cat("EXPECTED DOWNSTREAM EFFECTS:\n")
cat("  - Down-regulated genes: Direct KAT8 targets (lost H4K16ac -> silencing)\n")
cat("  - Up-regulated genes: Indirect effects (compensatory, secondary signaling)\n")
cat("  - TFs whose targets are DOWN: May require KAT8-mediated chromatin opening\n")
cat("  - TFs whose targets are UP: May be repressed by KAT8 or activated by stress\n\n")

## Check if KAT8 complex members appear as enriched TFs
if (length(tf_enrichment_results) > 0) {
  all_tf <- bind_rows(tf_enrichment_results, .id = "Contrast")
  nsl_msl_hits <- all_tf %>%
    filter(toupper(TF_name) %in% toupper(chromatin_modifiers$NSL_MSL_complex))

  if (nrow(nsl_msl_hits) > 0) {
    cat("NSL/MSL COMPLEX MEMBERS DETECTED AS ENRICHED TFS:\n")
    for (k in seq_len(nrow(nsl_msl_hits))) {
      row <- nsl_msl_hits[k, ]
      cat("  ", row$TF_name, " in ", row$Contrast,
          ": NES = ", round(row$NES, 2), "\n", sep = "")
    }
    cat("\nThis suggests KAT8 complex targets are coordinately affected.\n")
  } else {
    cat("No NSL/MSL complex members detected as enriched TFs.\n")
    cat("This is not unexpected - the analysis may capture downstream effectors.\n")
  }
}

## ============================================================
## 11) SESSION INFO
## ============================================================

session_file <- file.path(outdir, "logs", paste0("sessionInfo_tf_analysis_", run_tag, ".txt"))
if (dir.exists(dirname(session_file))) {
  sink(session_file)
  cat("=== R SESSION INFORMATION - Part 9 (TF Analysis) ===\n\n")
  print(sessionInfo())
  sink()
}

cat("\n=== Part 9 (Upstream Regulator / TF Analysis) COMPLETE ===\n")
cat("Output directory: ", tf_dir, "\n", sep = "")
