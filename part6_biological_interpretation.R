#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 6: BIOLOGICAL INTERPRETATION
## =========================================================
## This script integrates ORA, fGSEA, and cell type marker analysis
## to provide biological interpretation of bulk RNA-seq data
##
## KEY ANALYSES:
## 1. Pathway Convergence: Compare ORA vs fGSEA results
## 2. Cell Type Marker Overlay: Score adipose cell signatures
## 3. Leading Edge Integration: Which cells drive pathway signals?
## 4. Publication Figures: Integrated summary visualizations
##
## INTERPRETATION PHILOSOPHY:
## Bulk RNA-seq captures composite signals from multiple cell types.
## We CANNOT claim:
##   - "More macrophages" (only: "macrophage markers are elevated")
##   - "Pathway X is activated" (only: "genes in pathway X are upregulated")
##   - Causality (only: "associated with" or "consistent with")
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "purrr", "tibble", "RColorBrewer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "org.Mm.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "ComplexHeatmap"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(tibble)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(circlize)
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

## ============================================================
## CELL TYPE MARKER DEFINITIONS - ADIPOSE TISSUE
## ============================================================
## These are curated marker gene sets for adipose-resident cell types
## Sources: Tabula Muris, literature consensus, single-cell atlases
##
## IMPORTANT: These are SIGNATURE genes, not absolute markers
## Bulk RNA-seq enrichment ≠ cell abundance (could be cell state change)
## ============================================================

adipose_cell_markers <- list(

  ## Mature adipocytes - lipid metabolism, insulin signaling
  Adipocyte = c(
    "Adipoq", "Lep", "Fabp4", "Plin1", "Pparg",    ## Core adipocyte markers
    "Cidec", "Cidea", "Retn", "Rbp4", "Aqp7",      ## Lipid droplet/metabolism
    "Slc2a4", "Irs1", "Cebpa", "Lpl", "Acaca",     ## Insulin/glucose/lipogenesis
    "Fasn", "Scd1", "Dgat1", "Dgat2", "Pck1"       ## Lipid synthesis
  ),

  ## Adipocyte precursors / preadipocytes - fibroblast-like
  Preadipocyte = c(
    "Pdgfra", "Pdgfrb", "Cd34", "Ly6a", "Dcn",     ## Stem/progenitor markers
    "Col1a1", "Col1a2", "Col3a1", "Fn1", "Vim",    ## ECM/mesenchymal
    "Dlk1", "Pref1", "Sox9", "Zfp423", "Pparg"     ## Adipogenic commitment
  ),

  ## Macrophages - M1/M2 continuum
  Macrophage = c(
    "Adgre1", "Cd68", "Csf1r", "Itgam", "Mrc1",    ## Pan-macrophage (F4/80, CD68)
    "Cd163", "Lyve1", "Mertk", "Timd4", "Fcgr1",   ## Tissue-resident/M2-like
    "Ccr2", "Ly6c2", "Itgax", "Cd86", "Nos2",      ## Inflammatory/M1-like
    "Arg1", "Chil3", "Il10", "Tgfb1", "Cd206"      ## Anti-inflammatory/M2-like
  ),

  ## Pro-inflammatory macrophage subset (M1-like)
  Macrophage_M1 = c(
    "Nos2", "Tnf", "Il1b", "Il6", "Ccl2",          ## Pro-inflammatory cytokines
    "Cxcl9", "Cxcl10", "Cd86", "Cd80", "Fcgr1",    ## M1 polarization
    "Stat1", "Irf5", "Nfkb1", "Socs3", "Ptgs2"     ## M1 signaling
  ),

  ## Anti-inflammatory macrophage subset (M2-like)
  Macrophage_M2 = c(
    "Arg1", "Mrc1", "Cd163", "Chil3", "Retnla",    ## M2 markers
    "Il10", "Tgfb1", "Mertk", "Pparg", "Clec10a",  ## Anti-inflammatory
    "Stat6", "Irf4", "Ccl22", "Ccl17", "Fcer2a"    ## M2 signaling
  ),

  ## Endothelial cells - vascular
  Endothelial = c(
    "Pecam1", "Cdh5", "Vwf", "Kdr", "Flt1",        ## Core EC markers (CD31, VE-cad)
    "Tek", "Emcn", "Icam1", "Vcam1", "Sele",       ## Vascular/adhesion
    "Cldn5", "Plvap", "Aqp1", "Esam", "Egfl7"      ## Capillary EC
  ),

  ## T cells - adaptive immunity
  T_cell = c(
    "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a",         ## Pan-T and subsets
    "Lck", "Zap70", "Lat", "Trac", "Trbc1",        ## TCR signaling
    "Foxp3", "Il2ra", "Ctla4", "Il7r", "Ccr7"      ## Treg/naive markers
  ),

  ## B cells - adaptive immunity (minor in adipose)
  B_cell = c(
    "Cd19", "Cd79a", "Cd79b", "Ms4a1", "Pax5",     ## Pan-B markers
    "Ighm", "Cd22", "Blnk", "Bank1", "Ebf1"        ## B cell signaling
  ),

  ## Natural Killer cells
  NK_cell = c(
    "Ncr1", "Klrb1c", "Klrd1", "Nkg7", "Prf1",     ## NK markers (NKp46)
    "Gzma", "Gzmb", "Ifng", "Klrk1", "Fcgr3"       ## Cytotoxicity
  ),

  ## Fibroblasts / stromal cells
  Fibroblast = c(
    "Dcn", "Lum", "Col1a1", "Col1a2", "Col3a1",    ## ECM production
    "Pdgfra", "Fap", "Thy1", "Acta2", "S100a4",    ## Fibroblast markers
    "Vim", "Fn1", "Postn", "Tnc", "Sparc"          ## Matrix proteins
  ),

  ## Smooth muscle / pericytes
  Smooth_muscle = c(
    "Acta2", "Myh11", "Tagln", "Cnn1", "Des",      ## Smooth muscle actin/myosin
    "Pdgfrb", "Rgs5", "Notch3", "Mcam", "Kcnj8"    ## Pericyte/mural
  ),

  ## Dendritic cells
  Dendritic = c(
    "Itgax", "Cd74", "H2-Ab1", "H2-Eb1", "Flt3",   ## DC markers (CD11c, MHCII)
    "Batf3", "Irf8", "Xcr1", "Clec9a", "Cd209a"    ## DC subsets
  ),

  ## Mast cells
  Mast_cell = c(
    "Kit", "Cpa3", "Tpsb2", "Tpsab1", "Ms4a2",     ## Mast cell markers
    "Fcer1a", "Mcpt1", "Mcpt4", "Mcpt8", "Hdc"     ## Tryptase/chymase
  ),

  ## Neutrophils (typically transient in adipose)
  Neutrophil = c(
    "S100a8", "S100a9", "Ly6g", "Cxcr2", "Mmp9",   ## Neutrophil markers
    "Ltf", "Ngp", "Camp", "Elane", "Mpo"           ## Granule proteins
  ),

  ## Crown-like structure / dying adipocytes
  ## Marker of adipocyte stress/death
  Crown_like_structure = c(
    "Trem2", "Gpnmb", "Cd9", "Spp1", "Lpl",        ## CLS macrophage markers
    "Mmp12", "Ctsl", "Fabp5", "Lgals3", "Cd36"     ## Lipid-laden macrophages
  )
)

## ============================================================
## BIOLOGICAL PROCESS CATEGORIES
## ============================================================
## For grouping pathways into interpretable themes

pathway_themes <- list(
  Inflammation = c(
    "inflammatory", "cytokine", "chemokine", "interferon",
    "nf-kappa", "nfkb", "interleukin", "tnf", "il-1", "il-6",
    "innate immune", "acute phase"
  ),

  Lipid_Metabolism = c(
    "lipid", "fatty acid", "triglyceride", "cholesterol",
    "lipoprotein", "adipogenesis", "lipogenesis", "lipolysis",
    "beta-oxidation", "peroxisome"
  ),

  ECM_Remodeling = c(
    "extracellular matrix", "collagen", "matrix metalloproteinase",
    "ecm", "basement membrane", "integrin", "focal adhesion",
    "cell-matrix", "proteoglycan"
  ),

  Mitochondria = c(
    "mitochondr", "oxidative phosphorylation", "electron transport",
    "respiratory chain", "atp synthesis", "nad", "fadh"
  ),

  Fibrosis = c(
    "fibrosis", "fibrotic", "collagen", "tgf-beta", "tgfb",
    "wound healing", "tissue remodeling"
  ),

  Angiogenesis = c(
    "angiogenesis", "blood vessel", "vascular", "vegf",
    "endothelial", "vasculature"
  ),

  Immune_Cell_Migration = c(
    "chemotaxis", "leukocyte migration", "cell migration",
    "extravasation", "diapedesis", "homing"
  ),

  Apoptosis = c(
    "apoptosis", "apoptotic", "cell death", "programmed cell death",
    "caspase", "bcl-2"
  ),

  Insulin_Signaling = c(
    "insulin", "glucose", "glycolysis", "gluconeogenesis",
    "pi3k", "akt", "mtor"
  ),

  Oxidative_Stress = c(
    "oxidative stress", "reactive oxygen", "antioxidant",
    "glutathione", "nrf2", "ros"
  )
)

## ============================================================
## UTILITY FUNCTIONS
## ============================================================

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## Calculate signature score for a gene set
## Method: Mean log2FC of marker genes present in DEG data
calculate_signature_score <- function(de_table, marker_genes, method = "mean_logfc") {

  ## Find marker genes present in data
  present_markers <- intersect(rownames(de_table), marker_genes)

  if (length(present_markers) < 3) {
    return(list(
      score = NA,
      n_markers = length(present_markers),
      n_up = NA,
      n_down = NA,
      markers_found = present_markers,
      warning = "Too few markers found"
    ))
  }

  ## Get logFC values for present markers
  marker_lfc <- de_table[present_markers, "log2FoldChange"]
  marker_padj <- de_table[present_markers, "padj"]

  ## Calculate score based on method
  if (method == "mean_logfc") {
    score <- mean(marker_lfc, na.rm = TRUE)
  } else if (method == "median_logfc") {
    score <- median(marker_lfc, na.rm = TRUE)
  } else if (method == "ssgsea") {
    ## Simplified ssGSEA-like: sum of signed ranks
    all_lfc <- de_table$log2FoldChange
    ranks <- rank(all_lfc)
    names(ranks) <- rownames(de_table)
    score <- mean(ranks[present_markers], na.rm = TRUE) / length(all_lfc)
  }

  ## Count significant up/down
  sig_markers <- present_markers[!is.na(marker_padj) & marker_padj < 0.05]
  n_up <- sum(marker_lfc[sig_markers] > 0, na.rm = TRUE)
  n_down <- sum(marker_lfc[sig_markers] < 0, na.rm = TRUE)

  return(list(
    score = score,
    n_markers = length(present_markers),
    n_up = n_up,
    n_down = n_down,
    markers_found = present_markers,
    warning = NULL
  ))
}

## Categorize pathway by theme
categorize_pathway <- function(pathway_name, theme_keywords) {
  pathway_lower <- tolower(pathway_name)

  for (theme in names(theme_keywords)) {
    keywords <- theme_keywords[[theme]]
    if (any(sapply(keywords, function(kw) grepl(kw, pathway_lower, ignore.case = TRUE)))) {
      return(theme)
    }
  }
  return("Other")
}

## Fisher's exact test for marker overlap with DEGs
test_marker_enrichment <- function(deg_genes, marker_genes, background_genes) {

  ## Create contingency table
  deg_set <- intersect(deg_genes, background_genes)
  marker_set <- intersect(marker_genes, background_genes)

  a <- length(intersect(deg_set, marker_set))  ## DEG AND marker
  b <- length(setdiff(marker_set, deg_set))     ## marker but NOT DEG
  c <- length(setdiff(deg_set, marker_set))     ## DEG but NOT marker
  d <- length(background_genes) - a - b - c     ## neither

  ## Fisher's exact test
  mat <- matrix(c(a, b, c, d), nrow = 2)
  test <- fisher.test(mat, alternative = "greater")

  return(list(
    overlap = a,
    odds_ratio = test$estimate,
    p_value = test$p.value,
    marker_in_deg = a,
    total_markers = length(marker_set),
    fraction = if (length(marker_set) > 0) a / length(marker_set) else 0
  ))
}

## ============================================================
## MAIN ANALYSIS
## ============================================================

tryCatch({

  ## 3) Find Part 1 results -----------------------------------
  cat("\n=== FINDING PART 1/2/3 RESULTS ===\n")

  savepoint_dir <- file.path(getwd(), "savepoints")
  if (!dir.exists(savepoint_dir)) {
    stop("savepoints directory not found. Please run Part 1-3 first.")
  }

  run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
  if (length(run_dirs) == 0) {
    stop("No Part 1 run directories found in savepoints/")
  }

  run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
  outdir <- run_dirs_sorted[1]
  run_tag <- gsub("^RUN_", "", basename(outdir))

  cat("[INFO] Using results from: ", outdir, "\n", sep = "")
  cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

  tables_dir <- file.path(outdir, "tables")
  plots_dir <- file.path(outdir, "plots")

  ## Create interpretation subdirectory
  interp_dir <- file.path(plots_dir, "interpretation")
  if (!dir.exists(interp_dir)) dir.create(interp_dir, recursive = TRUE)

  ## ============================================================
  ## 4) LOAD ALL RESULTS
  ## ============================================================

  cat("\n=== LOADING DE, ORA, AND GSEA RESULTS ===\n")

  ## 4.1) Load DE tables
  de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
  de_tables <- list()
  for (de_file in de_files) {
    contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", basename(de_file))
    de_tables[[contrast_name]] <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)
    cat("[OK] Loaded DE table: ", contrast_name, " (", nrow(de_tables[[contrast_name]]), " genes)\n", sep = "")
  }

  ## 4.2) Load ORA results
  ora_files <- list.files(tables_dir, pattern = "^ORA_gobp_.*\\.csv$", full.names = TRUE)
  ora_results <- list()
  for (ora_file in ora_files) {
    bn <- basename(ora_file)
    ## Pattern: ORA_gobp_[contrast]_[direction]_[timestamp].csv
    parts <- strsplit(gsub("\\.csv$", "", bn), "_")[[1]]
    ## Find the Up or Down direction
    if ("Up" %in% parts) {
      direction <- "Up"
      idx <- which(parts == "Up")
      contrast_name <- paste(parts[3:(idx-1)], collapse = "_")
    } else if ("Down" %in% parts) {
      direction <- "Down"
      idx <- which(parts == "Down")
      contrast_name <- paste(parts[3:(idx-1)], collapse = "_")
    } else {
      next
    }
    key <- paste0(contrast_name, "_", direction)
    ora_results[[key]] <- read.csv(ora_file, stringsAsFactors = FALSE)
    cat("[OK] Loaded ORA: ", key, " (", nrow(ora_results[[key]]), " pathways)\n", sep = "")
  }

  ## 4.3) Load fGSEA results
  gsea_files <- list.files(tables_dir, pattern = "^fgsea_gobp_.*\\.csv$", full.names = TRUE)
  gsea_results <- list()
  for (gsea_file in gsea_files) {
    bn <- basename(gsea_file)
    contrast_name <- gsub("^fgsea_gobp_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", bn)
    gsea_results[[contrast_name]] <- read.csv(gsea_file, stringsAsFactors = FALSE)
    cat("[OK] Loaded fGSEA: ", contrast_name, " (", nrow(gsea_results[[contrast_name]]), " pathways)\n", sep = "")
  }

  ## ============================================================
  ## 5) ANALYSIS 1: PATHWAY CONVERGENCE (ORA vs fGSEA)
  ## ============================================================
  ## High-confidence pathways: significant in BOTH methods

  cat("\n")
  cat("============================================================\n")
  cat("=== ANALYSIS 1: PATHWAY CONVERGENCE (ORA vs fGSEA) ===\n")
  cat("============================================================\n")
  cat("\n")
  cat("RATIONALE:\n")
  cat("  - ORA uses discrete DEG lists (conservative, depends on cutoffs)\n")
  cat("  - fGSEA uses ranked gene lists (more sensitive to subtle changes)\n")
  cat("  - Pathways significant in BOTH methods = highest confidence\n")
  cat("  - Pathways in only one method may still be biologically relevant\n")
  cat("\n")

  convergence_results <- list()

  for (contrast in names(de_tables)) {
    cat("--- ", contrast, " ---\n", sep = "")

    ## Get ORA results for this contrast
    ora_up_key <- paste0(contrast, "_Up")
    ora_down_key <- paste0(contrast, "_Down")

    ora_up <- if (ora_up_key %in% names(ora_results)) ora_results[[ora_up_key]] else NULL
    ora_down <- if (ora_down_key %in% names(ora_results)) ora_results[[ora_down_key]] else NULL

    ## Get fGSEA results
    gsea_df <- if (contrast %in% names(gsea_results)) gsea_results[[contrast]] else NULL

    if (is.null(gsea_df)) {
      cat("  [WARN] No fGSEA results found\n")
      next
    }

    ## Split fGSEA by direction (NES sign)
    gsea_up <- gsea_df %>% filter(NES > 0, padj < 0.05)
    gsea_down <- gsea_df %>% filter(NES < 0, padj < 0.05)

    ## Find convergent pathways (in both ORA and fGSEA)
    ## Use pathway ID for matching
    ora_up_ids <- if (!is.null(ora_up)) ora_up$ID else character(0)
    ora_down_ids <- if (!is.null(ora_down)) ora_down$ID else character(0)
    gsea_up_ids <- gsea_up$pathway_id
    gsea_down_ids <- gsea_down$pathway_id

    convergent_up <- intersect(ora_up_ids, gsea_up_ids)
    convergent_down <- intersect(ora_down_ids, gsea_down_ids)

    cat("  Up-regulated pathways:\n")
    cat("    ORA only: ", length(setdiff(ora_up_ids, gsea_up_ids)), "\n", sep = "")
    cat("    fGSEA only: ", length(setdiff(gsea_up_ids, ora_up_ids)), "\n", sep = "")
    cat("    CONVERGENT (both): ", length(convergent_up), " <-- HIGH CONFIDENCE\n", sep = "")

    cat("  Down-regulated pathways:\n")
    cat("    ORA only: ", length(setdiff(ora_down_ids, gsea_down_ids)), "\n", sep = "")
    cat("    fGSEA only: ", length(setdiff(gsea_down_ids, ora_down_ids)), "\n", sep = "")
    cat("    CONVERGENT (both): ", length(convergent_down), " <-- HIGH CONFIDENCE\n", sep = "")

    ## Store convergent results
    convergence_results[[contrast]] <- list(
      convergent_up = convergent_up,
      convergent_down = convergent_down,
      ora_only_up = setdiff(ora_up_ids, gsea_up_ids),
      ora_only_down = setdiff(ora_down_ids, gsea_down_ids),
      gsea_only_up = setdiff(gsea_up_ids, ora_up_ids),
      gsea_only_down = setdiff(gsea_down_ids, ora_down_ids)
    )

    ## Save convergent pathway table
    if (length(convergent_up) > 0 || length(convergent_down) > 0) {

      ## Build convergent table with details from both methods
      convergent_table <- data.frame()

      if (length(convergent_up) > 0 && !is.null(ora_up)) {
        up_df <- ora_up %>%
          dplyr::filter(ID %in% convergent_up) %>%
          dplyr::mutate(Direction = "Up") %>%
          dplyr::left_join(
            gsea_up %>% dplyr::select(pathway_id, NES, padj) %>% dplyr::rename(gsea_NES = NES, gsea_padj = padj),
            by = c("ID" = "pathway_id")
          ) %>%
          dplyr::rename(ora_padj = p.adjust)
        convergent_table <- dplyr::bind_rows(convergent_table, up_df)
      }

      if (length(convergent_down) > 0 && !is.null(ora_down)) {
        down_df <- ora_down %>%
          dplyr::filter(ID %in% convergent_down) %>%
          dplyr::mutate(Direction = "Down") %>%
          dplyr::left_join(
            gsea_down %>% dplyr::select(pathway_id, NES, padj) %>% dplyr::rename(gsea_NES = NES, gsea_padj = padj),
            by = c("ID" = "pathway_id")
          ) %>%
          dplyr::rename(ora_padj = p.adjust)
        convergent_table <- dplyr::bind_rows(convergent_table, down_df)
      }

      if (nrow(convergent_table) > 0) {
        ## Add theme categorization
        convergent_table$Theme <- sapply(convergent_table$Description, categorize_pathway, theme_keywords = pathway_themes)

        conv_file <- paste0("convergent_pathways_", contrast, "_", run_tag, ".csv")
        write.csv(convergent_table, file = file.path(tables_dir, conv_file), row.names = FALSE)
        cat("  [OK] Saved convergent pathways: ", conv_file, "\n", sep = "")
      }
    }
    cat("\n")
  }

  ## ============================================================
  ## 6) ANALYSIS 2: CELL TYPE MARKER OVERLAY
  ## ============================================================
  ## Score each cell type's marker enrichment in DEGs

  cat("\n")
  cat("============================================================\n")
  cat("=== ANALYSIS 2: CELL TYPE MARKER OVERLAY ===\n")
  cat("============================================================\n")
  cat("\n")
  cat("RATIONALE:\n")
  cat("  - Bulk RNA-seq = composite signal from multiple cell types\n")
  cat("  - If macrophage markers are upregulated, this could indicate:\n")
  cat("    (a) More macrophages infiltrating the tissue\n")
  cat("    (b) Resident macrophages changing their expression state\n")
  cat("    (c) Other cells expressing macrophage-associated genes\n")
  cat("  - We CANNOT distinguish these with bulk data alone\n")
  cat("  - Language: 'macrophage signature is enriched' NOT 'more macrophages'\n")
  cat("\n")

  cell_type_scores <- list()

  for (contrast in names(de_tables)) {
    cat("--- ", contrast, " ---\n", sep = "")

    de_table <- de_tables[[contrast]]
    background_genes <- rownames(de_table)

    ## Get DEG thresholds for this contrast
    thresholds <- get_tissue_thresholds(contrast)

    ## Define DEGs
    deg_up <- rownames(de_table)[
      !is.na(de_table$padj) &
      de_table$padj < thresholds$fdr_cut &
      de_table$log2FoldChange > thresholds$logFC_cut
    ]
    deg_down <- rownames(de_table)[
      !is.na(de_table$padj) &
      de_table$padj < thresholds$fdr_cut &
      de_table$log2FoldChange < -thresholds$logFC_cut
    ]

    scores_list <- list()

    for (cell_type in names(adipose_cell_markers)) {
      markers <- adipose_cell_markers[[cell_type]]

      ## Calculate signature score (mean logFC of markers)
      sig_result <- calculate_signature_score(de_table, markers)

      ## Test enrichment in up-regulated DEGs
      enrich_up <- test_marker_enrichment(deg_up, markers, background_genes)

      ## Test enrichment in down-regulated DEGs
      enrich_down <- test_marker_enrichment(deg_down, markers, background_genes)

      scores_list[[cell_type]] <- data.frame(
        CellType = cell_type,
        Contrast = contrast,
        N_Markers_Defined = length(markers),
        N_Markers_Found = sig_result$n_markers,
        Mean_LogFC = sig_result$score,
        N_Sig_Up = sig_result$n_up,
        N_Sig_Down = sig_result$n_down,
        Enrich_Up_OR = enrich_up$odds_ratio,
        Enrich_Up_P = enrich_up$p_value,
        Enrich_Up_N = enrich_up$overlap,
        Enrich_Down_OR = enrich_down$odds_ratio,
        Enrich_Down_P = enrich_down$p_value,
        Enrich_Down_N = enrich_down$overlap,
        stringsAsFactors = FALSE
      )
    }

    cell_type_scores[[contrast]] <- dplyr::bind_rows(scores_list)

    ## Print summary
    ct_summary <- cell_type_scores[[contrast]] %>%
      dplyr::arrange(dplyr::desc(abs(Mean_LogFC)))

    cat("  Cell type signature scores (mean logFC of markers):\n")
    for (i in seq_len(min(5, nrow(ct_summary)))) {
      row <- ct_summary[i, ]
      dir_label <- if (!is.na(row$Mean_LogFC) && row$Mean_LogFC > 0) "UP" else "DOWN"
      cat("    ", row$CellType, ": ", round(row$Mean_LogFC, 3), " (", dir_label, ")",
          " [", row$N_Markers_Found, "/", row$N_Markers_Defined, " markers]\n", sep = "")
    }
    cat("\n")
  }

  ## Combine all cell type scores
  all_ct_scores <- bind_rows(cell_type_scores)
  ct_scores_file <- paste0("cell_type_signatures_", run_tag, ".csv")
  write.csv(all_ct_scores, file = file.path(tables_dir, ct_scores_file), row.names = FALSE)
  cat("[OK] Saved cell type signature scores: ", ct_scores_file, "\n", sep = "")

  ## ============================================================
  ## 7) ANALYSIS 3: LEADING EDGE × CELL TYPE INTEGRATION
  ## ============================================================
  ## Which cell types' markers appear in pathway leading edges?

  cat("\n")
  cat("============================================================\n")
  cat("=== ANALYSIS 3: LEADING EDGE × CELL TYPE INTEGRATION ===\n")
  cat("============================================================\n")
  cat("\n")
  cat("RATIONALE:\n")
  cat("  - Leading edge = core genes driving pathway enrichment\n")
  cat("  - If inflammatory pathway's leading edge contains macrophage markers:\n")
  cat("    -> The signal may be driven by macrophage-derived expression\n")
  cat("  - If it contains adipocyte markers:\n")
  cat("    -> Adipocytes themselves may be source of inflammation\n")
  cat("  - This helps interpret pathway enrichment mechanistically\n")
  cat("\n")

  le_celltype_results <- list()

  for (contrast in names(gsea_results)) {
    cat("--- ", contrast, " ---\n", sep = "")

    gsea_df <- gsea_results[[contrast]]

    ## Get top pathways
    top_pathways <- gsea_df %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = 20)

    if (nrow(top_pathways) == 0) {
      cat("  No significant pathways\n\n")
      next
    }

    pathway_celltype_list <- list()

    for (i in seq_len(nrow(top_pathways))) {
      pathway <- top_pathways[i, ]

      ## Parse leading edge genes
      le_genes <- character(0)
      if ("leadingEdge_symbols" %in% colnames(pathway) && !is.na(pathway$leadingEdge_symbols)) {
        le_genes <- unlist(strsplit(pathway$leadingEdge_symbols, ","))
      }

      if (length(le_genes) == 0) next

      ## Check overlap with each cell type
      celltype_overlap <- data.frame()
      for (cell_type in names(adipose_cell_markers)) {
        markers <- adipose_cell_markers[[cell_type]]
        overlap <- intersect(le_genes, markers)

        if (length(overlap) > 0) {
          celltype_overlap <- rbind(celltype_overlap, data.frame(
            Pathway = pathway$pathway_name,
            PathwayID = pathway$pathway_id,
            NES = pathway$NES,
            padj = pathway$padj,
            CellType = cell_type,
            N_Overlap = length(overlap),
            N_LeadingEdge = length(le_genes),
            Fraction = length(overlap) / length(le_genes),
            Overlap_Genes = paste(overlap, collapse = ", "),
            stringsAsFactors = FALSE
          ))
        }
      }

      if (nrow(celltype_overlap) > 0) {
        pathway_celltype_list[[pathway$pathway_id]] <- celltype_overlap
      }
    }

    if (length(pathway_celltype_list) > 0) {
      le_celltype_df <- dplyr::bind_rows(pathway_celltype_list)
      le_celltype_results[[contrast]] <- le_celltype_df

      ## Print notable findings
      notable <- le_celltype_df %>%
        dplyr::filter(N_Overlap >= 3) %>%
        dplyr::arrange(dplyr::desc(N_Overlap))

      if (nrow(notable) > 0) {
        cat("  Notable pathway-cell type overlaps (>=3 genes):\n")
        for (j in seq_len(min(5, nrow(notable)))) {
          row <- notable[j, ]
          cat("    ", substr(row$Pathway, 1, 40), "...\n", sep = "")
          cat("      ", row$CellType, ": ", row$N_Overlap, "/", row$N_LeadingEdge,
              " (", round(100*row$Fraction, 1), "%)\n", sep = "")
        }
      }
    }
    cat("\n")
  }

  ## Save leading edge × cell type results
  if (length(le_celltype_results) > 0) {
    all_le_ct <- bind_rows(le_celltype_results, .id = "Contrast")
    le_ct_file <- paste0("leading_edge_celltype_", run_tag, ".csv")
    write.csv(all_le_ct, file = file.path(tables_dir, le_ct_file), row.names = FALSE)
    cat("[OK] Saved leading edge × cell type analysis: ", le_ct_file, "\n", sep = "")
  }

  ## ============================================================
  ## 8) PUBLICATION FIGURE 1: CELL TYPE SIGNATURE HEATMAP
  ## ============================================================

  cat("\n=== GENERATING PUBLICATION FIGURES ===\n\n")

  ## Figure 1: Heatmap of cell type signature scores across contrasts
  if (nrow(all_ct_scores) > 0) {

    cat("[INFO] Creating cell type signature heatmap...\n")

    ## Pivot to matrix
    ct_matrix <- all_ct_scores %>%
      dplyr::select(CellType, Contrast, Mean_LogFC) %>%
      tidyr::pivot_wider(names_from = Contrast, values_from = Mean_LogFC) %>%
      tibble::column_to_rownames("CellType")

    ## Handle NAs
    ct_matrix[is.na(ct_matrix)] <- 0

    ## Order cell types by biological grouping
    celltype_order <- c(
      "Adipocyte", "Preadipocyte",
      "Macrophage", "Macrophage_M1", "Macrophage_M2", "Crown_like_structure",
      "T_cell", "B_cell", "NK_cell", "Dendritic", "Neutrophil", "Mast_cell",
      "Endothelial", "Fibroblast", "Smooth_muscle"
    )
    celltype_order <- intersect(celltype_order, rownames(ct_matrix))
    ct_matrix <- ct_matrix[celltype_order, , drop = FALSE]

    ## Create color scale (diverging, centered at 0)
    ## Using viridis blue and teal (avoiding extreme purple/yellow)
    max_abs <- max(abs(ct_matrix), na.rm = TRUE)
    col_fun <- colorRamp2(
      c(-max_abs, 0, max_abs),
      c("#3b528b", "white", "#21918c")  ## Viridis blue -> white -> viridis teal
    )

    ## Create cell type group annotation
    celltype_groups <- c(
      "Adipocyte" = "Adipocyte",
      "Preadipocyte" = "Adipocyte",
      "Macrophage" = "Immune",
      "Macrophage_M1" = "Immune",
      "Macrophage_M2" = "Immune",
      "Crown_like_structure" = "Immune",
      "T_cell" = "Immune",
      "B_cell" = "Immune",
      "NK_cell" = "Immune",
      "Dendritic" = "Immune",
      "Neutrophil" = "Immune",
      "Mast_cell" = "Immune",
      "Endothelial" = "Stromal",
      "Fibroblast" = "Stromal",
      "Smooth_muscle" = "Stromal"
    )

    row_groups <- celltype_groups[rownames(ct_matrix)]

    row_ha <- rowAnnotation(
      Category = row_groups,
      col = list(Category = c(
        "Adipocyte" = "#FFC107",
        "Immune" = "#E91E63",
        "Stromal" = "#4CAF50"
      )),
      show_legend = TRUE,
      annotation_name_side = "top"
    )

    ## Create heatmap
    ht <- Heatmap(
      as.matrix(ct_matrix),
      name = "Mean logFC\nof markers",
      col = col_fun,
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "left",
      column_names_rot = 45,
      left_annotation = row_ha,
      row_names_gp = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 11),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9)
      ),
      column_title = "Cell Type Signature Enrichment Across Contrasts",
      column_title_gp = gpar(fontsize = 14, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        ## Add asterisks for significant enrichment
        val <- ct_matrix[i, j]
        if (!is.na(val) && abs(val) > 0.3) {
          grid.text("*", x, y, gp = gpar(fontsize = 10, col = "black"))
        }
      }
    )

    ## Save
    png(file.path(interp_dir, paste0("Fig_celltype_signatures_", run_tag, ".png")),
        width = 10, height = 10, units = "in", res = 300)
    draw(ht, padding = unit(c(2, 2, 2, 15), "mm"))
    dev.off()
    cat("[OK] Saved: Fig_celltype_signatures_", run_tag, ".png\n", sep = "")
  }

  ## ============================================================
  ## 9) PUBLICATION FIGURE 2: PATHWAY CONVERGENCE SUMMARY
  ## ============================================================

  ## Figure 2: Venn-like summary of ORA vs fGSEA convergence
  if (length(convergence_results) > 0) {

    cat("[INFO] Creating pathway convergence summary...\n")

    ## Build summary data
    conv_summary <- data.frame()
    for (contrast in names(convergence_results)) {
      res <- convergence_results[[contrast]]
      conv_summary <- rbind(conv_summary, data.frame(
        Contrast = contrast,
        Direction = c("Up", "Up", "Up", "Down", "Down", "Down"),
        Category = rep(c("Convergent", "ORA only", "GSEA only"), 2),
        N_Pathways = c(
          length(res$convergent_up),
          length(res$ora_only_up),
          length(res$gsea_only_up),
          length(res$convergent_down),
          length(res$ora_only_down),
          length(res$gsea_only_down)
        ),
        stringsAsFactors = FALSE
      ))
    }

    ## Create stacked bar plot
    conv_summary$Category <- factor(conv_summary$Category,
                                    levels = c("GSEA only", "Convergent", "ORA only"))

    p_conv <- ggplot(conv_summary, aes(x = Contrast, y = N_Pathways, fill = Category)) +
      geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.3) +
      facet_wrap(~ Direction, scales = "free_y") +
      scale_fill_manual(
        values = c(
          "Convergent" = "#4CAF50",     ## Green = high confidence
          "ORA only" = "#FF9800",        ## Orange = ORA-specific
          "GSEA only" = "#2196F3"        ## Blue = GSEA-specific
        ),
        name = "Method Agreement"
      ) +
      labs(
        title = "Pathway Analysis Convergence: ORA vs fGSEA",
        subtitle = "Convergent pathways (green) have highest confidence",
        x = NULL,
        y = "Number of Significant Pathways (GO:BP, FDR < 0.05)"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 11)
      )

    ggsave(file.path(interp_dir, paste0("Fig_pathway_convergence_", run_tag, ".png")),
           plot = p_conv, width = 10, height = 6, dpi = 300, bg = "white")
    cat("[OK] Saved: Fig_pathway_convergence_", run_tag, ".png\n", sep = "")
  }

  ## ============================================================
  ## 10) PUBLICATION FIGURE 3: INTEGRATED PATHWAY × CELL TYPE
  ## ============================================================

  ## Figure 3: For top pathways, show which cell types contribute
  if (length(le_celltype_results) > 0) {

    cat("[INFO] Creating pathway × cell type integration plot...\n")

    all_le_ct <- dplyr::bind_rows(le_celltype_results, .id = "Contrast")

    ## Filter to notable overlaps
    notable_le <- all_le_ct %>%
      dplyr::filter(N_Overlap >= 2) %>%
      dplyr::mutate(
        Pathway_Short = substr(Pathway, 1, 50),
        Direction = ifelse(NES > 0, "Up", "Down")
      )

    if (nrow(notable_le) > 10) {
      ## Create tile plot
      ## Summarize: for each pathway-celltype, average fraction across contrasts
      tile_data <- notable_le %>%
        dplyr::group_by(Pathway_Short, CellType) %>%
        dplyr::summarise(
          Mean_Fraction = mean(Fraction, na.rm = TRUE),
          Total_Overlap = sum(N_Overlap),
          .groups = "drop"
        ) %>%
        dplyr::filter(Mean_Fraction > 0.05)  ## At least 5% of leading edge

      if (nrow(tile_data) > 5) {
        ## Get top pathways by total cell type overlap
        top_pathways <- tile_data %>%
          dplyr::group_by(Pathway_Short) %>%
          dplyr::summarise(Total = sum(Total_Overlap)) %>%
          dplyr::arrange(dplyr::desc(Total)) %>%
          dplyr::slice_head(n = 15) %>%
          dplyr::pull(Pathway_Short)

        tile_data_filtered <- tile_data %>%
          dplyr::filter(Pathway_Short %in% top_pathways)

        p_tile <- ggplot(tile_data_filtered, aes(x = CellType, y = Pathway_Short, fill = Mean_Fraction)) +
          geom_tile(color = "white", linewidth = 0.5) +
          geom_text(aes(label = Total_Overlap), size = 3, color = "black") +
          scale_fill_gradient(
            low = "white",
            high = "#E91E63",
            name = "Fraction of\nLeading Edge",
            labels = scales::percent
          ) +
          labs(
            title = "Cell Type Markers in Pathway Leading Edges",
            subtitle = "Numbers show count of overlapping genes",
            x = "Cell Type",
            y = "Pathway"
          ) +
          theme_minimal(base_size = 11) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = element_text(size = 8),
            panel.grid = element_blank(),
            legend.position = "right"
          )

        ggsave(file.path(interp_dir, paste0("Fig_pathway_celltype_", run_tag, ".png")),
               plot = p_tile, width = 12, height = 8, dpi = 300, bg = "white")
        cat("[OK] Saved: Fig_pathway_celltype_", run_tag, ".png\n", sep = "")
      }
    }
  }

  ## ============================================================
  ## 11) INTERPRETATION GUIDELINES
  ## ============================================================

  cat("\n")
  cat("============================================================\n")
  cat("=== INTERPRETATION GUIDELINES FOR MANUSCRIPT ===\n")
  cat("============================================================\n")
  cat("\n")
  cat("LANGUAGE TO USE:\n")
  cat("  GOOD: 'Macrophage marker genes were enriched among upregulated genes'\n")
  cat("  BAD:  'Macrophage abundance was increased'\n")
  cat("\n")
  cat("  GOOD: 'Genes associated with inflammatory processes were upregulated'\n")
  cat("  BAD:  'The inflammatory pathway was activated'\n")
  cat("\n")
  cat("  GOOD: 'These findings are consistent with increased immune infiltration'\
")
  cat("  BAD:  'KAT8 knockdown causes macrophage infiltration'\n")
  cat("\n")
  cat("KEY CAVEATS TO ACKNOWLEDGE:\n")
  cat("  1. Bulk RNA-seq cannot distinguish cell composition from cell state\n")
  cat("  2. Pathway enrichment reflects coordinated gene expression, not activity\n")
  cat("  3. Causality cannot be inferred from transcriptomic associations\n")
  cat("  4. Validation with orthogonal methods (histology, flow cytometry,\n")
  cat("     single-cell RNA-seq) would strengthen compositional claims\n")
  cat("\n")
  cat("SUGGESTED MANUSCRIPT FRAMING:\n")
  cat("  Results section:\n")
  cat("    'To characterize the transcriptomic changes, we performed pathway\n")
  cat("    enrichment analysis using both ORA and GSEA. Pathways showing\n")
  cat("    enrichment in both methods (convergent pathways) were considered\n")
  cat("    high confidence. To assess whether observed changes might reflect\n")
  cat("    altered cellular composition, we scored cell type marker signatures.'\n")
  cat("\n")
  cat("  Discussion section:\n")
  cat("    'The enrichment of macrophage marker genes among upregulated genes\n")
  cat("    is consistent with either increased macrophage abundance or altered\n")
  cat("    macrophage activation state. Distinguishing these possibilities would\n")
  cat("    require single-cell resolution or histological quantification.'\n")
  cat("\n")
  cat("============================================================\n")

  ## ============================================================
  ## 12) SESSION INFO & COMPLETION
  ## ============================================================

  cat("\n=== PART 6 (BIOLOGICAL INTERPRETATION) COMPLETE ===\n\n")
  cat("Output files generated:\n")
  cat("  Tables:\n")
  cat("    - convergent_pathways_[contrast]_*.csv (high-confidence pathways)\n")
  cat("    - cell_type_signatures_*.csv (marker enrichment scores)\n")
  cat("    - leading_edge_celltype_*.csv (pathway-cell type integration)\n")
  cat("\n")
  cat("  Figures (in plots/interpretation/):\n")
  cat("    - Fig_celltype_signatures_*.png (heatmap of cell type scores)\n")
  cat("    - Fig_pathway_convergence_*.png (ORA vs fGSEA agreement)\n")
  cat("    - Fig_pathway_celltype_*.png (pathway × cell type tile plot)\n")
  cat("\n")

  ## Save session info
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_interpretation_", run_tag, ".txt"))
  sink(session_file)
  cat("=== R SESSION INFORMATION - Part 6 (Biological Interpretation) ===\n\n")
  print(sessionInfo())
  sink()
  cat("[INFO] Session info saved to: ", session_file, "\n", sep = "")

}, error = function(e) {

  cat("\n======================================\n")
  cat("=== ERROR IN INTERPRETATION PIPELINE ===\n")
  cat("======================================\n")
  cat(conditionMessage(e), "\n\n")
  print(traceback())

  stop(e)
})

cat("\n=== Part 6 finished ===\n")
