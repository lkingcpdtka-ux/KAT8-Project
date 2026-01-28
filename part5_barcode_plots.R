#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 5: BARCODE PLOTS FOR KEY PATHWAYS
## =========================================================
## Creates GSEA-style barcode (enrichment) plots for selected
## pathways of interest identified from ORA/GSEA analysis
##
## Run AFTER parts 1-3 have completed
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "tidyr", "stringr", "scales")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "fgsea",
  "AnnotationDbi",
  "msigdbr"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(scales)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(fgsea)
  library(AnnotationDbi)
  library(msigdbr)
})

## 2) Find most recent run directory -----------------------
cat("\n=== FINDING ANALYSIS RESULTS ===\n")

savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run parts 1-3 first.")
}

run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No run directories found in savepoints/")
}

run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")

## Create barcode plots subdirectory
barcode_dir <- file.path(plots_dir, "barcode_plots")
if (!dir.exists(barcode_dir)) {
  dir.create(barcode_dir, recursive = TRUE)
}

## ============================================================
## SECTION 1: DEFINE PATHWAYS OF INTEREST
## ============================================================
## Specify pathways to create barcode plots for
## These are matched using partial/fuzzy matching against
## GO:BP and KEGG pathway names
## ============================================================

cat("\n=== DEFINING PATHWAYS OF INTEREST ===\n")

## User-specified pathways of interest
## Format: list with display name and search patterns
pathways_of_interest <- list(

  "Mitochondrial Gene Expression" = list(
    patterns = c(
      "mitochondrial gene expression",
      "mitochondrial translation",
      "mitochondrial transcription",
      "mitochondrial RNA"
    ),
    databases = c("GO:BP", "KEGG")
  ),

  "Cytokine Receptor" = list(
    patterns = c(
      "cytokine receptor",
      "cytokine-cytokine receptor",
      "cytokine mediated signaling",
      "cytokine signaling"
    ),
    databases = c("GO:BP", "KEGG")
  ),

  "Apoptosis" = list(
    patterns = c(
      "apoptosis",
      "apoptotic process",
      "programmed cell death",
      "apoptotic signaling"
    ),
    databases = c("GO:BP", "KEGG")
  ),

  "ATP-dependent Chromatin Remodeling" = list(
    patterns = c(
      "atp-dependent chromatin remodeling",
      "chromatin remodeling",
      "atp dependent chromatin",
      "nucleosome remodeling",
      "swi/snf"
    ),
    databases = c("GO:BP", "KEGG")
  ),

  "ECM Receptor Interaction" = list(
    patterns = c(
      "ecm-receptor interaction",
      "ecm receptor interaction",
      "extracellular matrix receptor",
      "cell-matrix adhesion",
      "integrin"
    ),
    databases = c("GO:BP", "KEGG")
  ),

  "Chemical/ROS Response" = list(
    patterns = c(
      "reactive oxygen species",
      "response to oxidative stress",
      "ros",
      "oxidative stress",
      "cellular response to chemical",
      "chemical carcinogenesis - reactive oxygen"
    ),
    databases = c("GO:BP", "KEGG")
  )

)

cat("[INFO] Pathways to analyze:\n")
for (pw_name in names(pathways_of_interest)) {
  cat("  - ", pw_name, "\n", sep = "")
}

## ============================================================
## SECTION 2: LOAD GENE SETS AND DE RESULTS
## ============================================================

cat("\n=== LOADING GENE SETS ===\n")

## Load GO:BP gene sets
go_bp_data <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
go_bp_list <- split(go_bp_data$gene_symbol, go_bp_data$gs_name)
cat("[OK] Loaded GO:BP gene sets: ", length(go_bp_list), " terms\n", sep = "")

## Load KEGG gene sets
kegg_data <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
kegg_list <- split(kegg_data$gene_symbol, kegg_data$gs_name)
cat("[OK] Loaded KEGG gene sets: ", length(kegg_list), " pathways\n", sep = "")

## Combine for searching
all_gene_sets <- c(go_bp_list, kegg_list)

## Load DE tables
de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
cat("[OK] Found ", length(de_files), " DE tables\n", sep = "")

## ============================================================
## SECTION 3: BARCODE PLOT FUNCTION
## ============================================================

create_barcode_plot <- function(ranked_genes, gene_set, pathway_name,
                                 contrast_name, direction_color = "auto",
                                 output_dir = barcode_dir) {

  ## ranked_genes: named numeric vector (gene names -> ranking metric)
  ## gene_set: character vector of gene symbols in the pathway
  ## pathway_name: string for plot title
  ## contrast_name: string for subtitle
  ## direction_color: "up" (orange), "down" (blue), or "auto" (based on ES)

  if (length(gene_set) == 0) {
    cat("[WARN] Empty gene set for: ", pathway_name, "\n")
    return(NULL)
  }

  ## Find pathway genes in ranked list
  pathway_genes_in_data <- intersect(gene_set, names(ranked_genes))

  if (length(pathway_genes_in_data) < 3) {
    cat("[WARN] Too few pathway genes found (", length(pathway_genes_in_data),
        ") for: ", pathway_name, "\n", sep = "")
    return(NULL)
  }

  ## Sort ranked genes
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  n_genes <- length(ranked_genes)

  ## Get positions of pathway genes in ranked list
  gene_positions <- match(pathway_genes_in_data, names(ranked_genes))
  gene_positions <- gene_positions[!is.na(gene_positions)]

  if (length(gene_positions) == 0) {
    cat("[WARN] No pathway genes found in ranked list for: ", pathway_name, "\n")
    return(NULL)
  }

  ## Calculate running enrichment score (simplified GSEA)
  ## Uses the Kolmogorov-Smirnov-like statistic
  n_hit <- length(gene_positions)
  n_miss <- n_genes - n_hit

  ## Calculate step sizes
  hit_step <- 1 / n_hit
  miss_step <- 1 / n_miss

  ## Calculate running score
  running_score <- numeric(n_genes)
  current_score <- 0

  for (i in seq_len(n_genes)) {
    if (i %in% gene_positions) {
      current_score <- current_score + hit_step
    } else {
      current_score <- current_score - miss_step
    }
    running_score[i] <- current_score
  }

  ## Find enrichment score (max deviation from 0)
  es_max <- max(running_score)
  es_min <- min(running_score)
  enrichment_score <- ifelse(abs(es_max) > abs(es_min), es_max, es_min)

  ## Determine direction and color
  if (direction_color == "auto") {
    direction_color <- ifelse(enrichment_score > 0, "up", "down")
  }

  if (direction_color == "up") {
    line_color <- "#E69F00"  ## Orange
    fill_color <- "#FDD49E"  ## Light orange
    direction_label <- "Enriched in Up-regulated"
  } else {
    line_color <- "#0072B2"  ## Blue
    fill_color <- "#9ECAE1"  ## Light blue
    direction_label <- "Enriched in Down-regulated"
  }

  ## Create data frame for plotting
  plot_df <- data.frame(
    rank = seq_len(n_genes),
    score = running_score
  )

  ## Create barcode data
  barcode_df <- data.frame(
    rank = gene_positions,
    y = 0
  )

  ## Create the plot
  p <- ggplot() +
    ## Running enrichment score line
    geom_line(data = plot_df, aes(x = rank, y = score),
              color = line_color, linewidth = 1.2) +
    ## Fill area under curve (optional)
    geom_area(data = plot_df, aes(x = rank, y = score),
              fill = fill_color, alpha = 0.3) +
    ## Zero line
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ## Barcode (gene hits)
    geom_segment(data = barcode_df,
                 aes(x = rank, xend = rank, y = -0.05, yend = 0.05),
                 color = line_color, linewidth = 0.3, alpha = 0.7) +
    ## Add enrichment score annotation
    annotate("text", x = n_genes * 0.98, y = enrichment_score * 0.9,
             label = paste0("ES = ", round(enrichment_score, 3)),
             hjust = 1, size = 3.5, fontface = "bold", color = line_color) +
    ## Add gene count annotation
    annotate("text", x = n_genes * 0.98, y = enrichment_score * 0.75,
             label = paste0("Genes: ", length(gene_positions), "/", length(gene_set)),
             hjust = 1, size = 3, color = "grey40") +
    ## Labels
    labs(
      title = str_wrap(pathway_name, width = 60),
      subtitle = paste0(contrast_name, " | ", direction_label),
      x = "Gene Rank",
      y = "Enrichment Score"
    ) +
    ## Theme
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, lineheight = 1.1),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      plot.margin = margin(10, 15, 10, 10)
    ) +
    ## Scale
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.02)),
      labels = comma
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.15, 0.15))
    )

  ## Add indicator for high/low rank regions
  p <- p +
    annotate("text", x = n_genes * 0.02, y = min(running_score) - 0.1,
             label = "Up-regulated", hjust = 0, size = 2.5, color = "#E69F00") +
    annotate("text", x = n_genes * 0.98, y = min(running_score) - 0.1,
             label = "Down-regulated", hjust = 1, size = 2.5, color = "#0072B2")

  ## Save plot
  safe_pathway_name <- gsub("[^A-Za-z0-9_-]", "_", pathway_name)
  safe_pathway_name <- gsub("_+", "_", safe_pathway_name)
  safe_pathway_name <- substr(safe_pathway_name, 1, 50)

  plot_file <- paste0("Barcode_", safe_pathway_name, "_", contrast_name, "_", run_tag, ".png")
  ggsave(
    file.path(output_dir, plot_file),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  cat("[OK] Saved barcode plot: ", plot_file, "\n", sep = "")

  return(list(
    plot = p,
    enrichment_score = enrichment_score,
    n_genes_found = length(gene_positions),
    n_genes_total = length(gene_set)
  ))
}

## ============================================================
## SECTION 4: MATCH PATHWAYS AND CREATE PLOTS
## ============================================================

cat("\n=== CREATING BARCODE PLOTS ===\n")

## Function to find matching gene sets
find_matching_gene_sets <- function(patterns, gene_set_list) {
  matches <- list()

  for (gs_name in names(gene_set_list)) {
    gs_name_lower <- tolower(gs_name)

    for (pattern in patterns) {
      pattern_lower <- tolower(pattern)

      if (grepl(pattern_lower, gs_name_lower, fixed = TRUE)) {
        matches[[gs_name]] <- gene_set_list[[gs_name]]
        break
      }
    }
  }

  return(matches)
}

## Summary tracker
barcode_summary <- data.frame(
  Pathway_Interest = character(),
  Matched_Pathway = character(),
  Contrast = character(),
  ES = numeric(),
  N_Genes_Found = integer(),
  N_Genes_Total = integer(),
  stringsAsFactors = FALSE
)

## Process each DE file
for (de_file in de_files) {

  de_filename <- basename(de_file)
  contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

  cat("\n--- Processing contrast: ", contrast_name, " ---\n", sep = "")

  ## Load DE table
  de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

  ## Ensure we have the stat column (Wald statistic for ranking)
  if (!"stat" %in% colnames(de_table)) {
    if ("log2FoldChange" %in% colnames(de_table)) {
      cat("[INFO] Using log2FoldChange for ranking (stat column not found)\n")
      de_table$stat <- de_table$log2FoldChange
    } else {
      cat("[WARN] Cannot find ranking metric in ", de_file, "\n")
      next
    }
  }

  ## Create ranked gene list
  ranked_genes <- de_table %>%
    dplyr::filter(is.finite(stat)) %>%
    dplyr::arrange(dplyr::desc(stat))

  ranked_vec <- setNames(ranked_genes$stat, rownames(ranked_genes))

  cat("[INFO] Ranked gene list: ", length(ranked_vec), " genes\n", sep = "")

  ## Process each pathway of interest
  for (pw_name in names(pathways_of_interest)) {

    pw_info <- pathways_of_interest[[pw_name]]
    patterns <- pw_info$patterns

    cat("\n  Searching for: ", pw_name, "\n", sep = "")

    ## Search in all gene sets
    matched_sets <- find_matching_gene_sets(patterns, all_gene_sets)

    if (length(matched_sets) == 0) {
      cat("    [WARN] No matching gene sets found\n")
      next
    }

    cat("    [OK] Found ", length(matched_sets), " matching gene sets\n", sep = "")

    ## Create barcode plot for each matched set (limit to top 3 by gene overlap)
    ## Calculate gene overlap for each matched set
    overlap_scores <- sapply(matched_sets, function(gs) {
      length(intersect(gs, names(ranked_vec)))
    })

    ## Sort by overlap and take top matches
    top_matches <- names(sort(overlap_scores, decreasing = TRUE))[1:min(3, length(overlap_scores))]

    for (gs_name in top_matches) {
      gene_set <- matched_sets[[gs_name]]

      ## Clean up pathway name for display
      display_name <- gsub("^GOBP_|^KEGG_|^GO_", "", gs_name)
      display_name <- gsub("_", " ", display_name)
      display_name <- str_to_title(display_name)

      result <- create_barcode_plot(
        ranked_genes = ranked_vec,
        gene_set = gene_set,
        pathway_name = display_name,
        contrast_name = contrast_name,
        direction_color = "auto",
        output_dir = barcode_dir
      )

      if (!is.null(result)) {
        barcode_summary <- rbind(
          barcode_summary,
          data.frame(
            Pathway_Interest = pw_name,
            Matched_Pathway = gs_name,
            Contrast = contrast_name,
            ES = result$enrichment_score,
            N_Genes_Found = result$n_genes_found,
            N_Genes_Total = result$n_genes_total,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

## ============================================================
## SECTION 5: SAVE SUMMARY AND FINISH
## ============================================================

cat("\n=== SAVING SUMMARY ===\n")

## Save summary table
summary_file <- paste0("barcode_plots_summary_", run_tag, ".csv")
write.csv(barcode_summary, file = file.path(tables_dir, summary_file), row.names = FALSE)
cat("[OK] Saved summary table: ", summary_file, "\n", sep = "")

## Print summary
cat("\n==========================================================\n")
cat("=== BARCODE PLOTS SUMMARY ===\n")
cat("==========================================================\n")
cat("Total plots created: ", nrow(barcode_summary), "\n\n", sep = "")

for (pw_name in unique(barcode_summary$Pathway_Interest)) {
  subset_data <- barcode_summary[barcode_summary$Pathway_Interest == pw_name, ]
  cat("--- ", pw_name, " ---\n", sep = "")
  for (i in seq_len(nrow(subset_data))) {
    row <- subset_data[i, ]
    direction <- ifelse(row$ES > 0, "Up", "Down")
    cat("  ", row$Contrast, ": ES=", round(row$ES, 3), " (", direction,
        "), Genes=", row$N_Genes_Found, "/", row$N_Genes_Total, "\n", sep = "")
  }
  cat("\n")
}
cat("==========================================================\n")

cat("\n======================================\n")
cat("=== PART 5 (BARCODE PLOTS) COMPLETE ===\n")
cat("======================================\n")
cat("Created:\n")
cat("  - GSEA-style barcode plots for pathways of interest\n")
cat("  - Summary table: ", summary_file, "\n", sep = "")
cat("\nAll plots saved to: ", barcode_dir, "\n\n", sep = "")
