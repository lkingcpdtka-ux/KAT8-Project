## =========================================================
## KAT8 bulk RNA-seq - CENTRAL PARAMETERS FILE
## =========================================================
## All analysis cutoffs and thresholds are defined HERE ONLY
## Other scripts (part1-6) source this file to use these values
## =========================================================

## ============================================================
## REPRODUCIBILITY: RANDOM SEED CONFIGURATION
## ============================================================
## Single seed for all stochastic operations to ensure reproducibility
## Used by: all parts (DESeq2, GSEA, clustering, PCA, etc.)

MASTER_SEED <- 12345
set.seed(MASTER_SEED)

## ============================================================
## DEG THRESHOLDS (Differential Expression)
## ============================================================
## Used by: part1_main_analysis.R, part2_ora.R, part3_fgsea.R

## iWAT cutoffs (stringent - strong signal)
iWAT_logFC_cut <- 1.0
iWAT_fdr_cut   <- 0.05

## gWAT cutoffs (can be relaxed if needed)
gWAT_logFC_cut <- 1.0
gWAT_fdr_cut   <- 0.05

## Default cutoffs (for any other contrasts)
default_logFC_cut <- 1.0
default_fdr_cut   <- 0.05

## ============================================================
## ORA PARAMETERS (Over-Representation Analysis)
## ============================================================
## Used by: part2_ora.R, part4_publication_plots.R

ora_params <- list(
  ## Enrichment analysis cutoffs (initial filtering)
  pvalue_cutoff   = 0.1,
  qvalue_cutoff   = 0.2,

  ## Gene set size constraints
  min_gs_size     = 5,
  max_gs_size     = 500,

  ## GO:BP simplification (removes redundant terms)
  simplify_cutoff = 0.7,

  ## Plotting parameters
  top_n           = 15,      ## Top pathways per direction for bar plots
  order_by        = "adj. p-value"
)

## ============================================================
## GSEA/fGSEA PARAMETERS (Gene Set Enrichment Analysis)
## ============================================================
## Used by: part3_fgsea.R, part5_barcode_plots.R

gsea_params <- list(
  ## Gene set size constraints
  min_gs_size     = 5,
  max_gs_size     = 500,

  ## Ranking metric (from DESeq2)
  rank_metric     = "DESeq2 Wald statistic",

  ## GO:BP simplification - reduces redundant GO terms using semantic similarity
  ## WARNING: simplify() is O(n²) complexity - can hang with many terms!
  ## Strategy: Only simplify the TOP N most significant terms (by p-value)
  simplify_go     = TRUE,    ## Enable GO:BP simplification
  simplify_cutoff = 0.5,     ## Semantic similarity threshold (lower = more aggressive merging)
  max_terms_for_simplify = 300,  ## Simplify top N terms by p-value (keeps computation manageable)

  ## Plotting parameters
  top_n_per_direction = 10,  ## Top pathways per direction for plots

  ## FDR cutoff for significance
  fdr_cutoff      = 0.05,

  ## Random seed for reproducibility (used by gseGO/gseKEGG)
  seed            = 12345
)

## ============================================================
## GSEA PLOT PARAMETERS
## ============================================================
## Used by: part4_visualizations.R

gsea_plot_params <- list(
  plot_width  = 7.5,
  plot_height = 4.8
)

## ============================================================
## PUBLICATION PLOT PARAMETERS
## ============================================================
## Used by: part4_publication_plots.R

dotplot_params <- list(
  ## Filtering
  fdr_cutoff      = 0.05,
  min_count       = 3,       ## Minimum genes per pathway

  ## Top pathways to show
  top_n_gobp      = 15,
  top_n_kegg      = 15,

  ## Plot sizing (fixed to keep dot plots consistent)
  plot_width      = 4.6,
  plot_height     = 4.6,

  ## GO:BP simplification cutoff (should match ora_params)
  gobp_simplify_cutoff = 0.7
)

## ============================================================
## ENRICHMENT BAR PLOT PARAMETERS
## ============================================================
## Used by: part4_visualizations.R

barplot_params <- list(
  plot_width  = 7.5,
  plot_height = 5.2
)

## ============================================================
## HEATMAP PARAMETERS
## ============================================================
## Used by: part4_visualizations.R

heatmap_params <- list(
  lfc_clip    = 2.0,         ## Clip z-scores at +/- this value

  ## Custom gradient: Green (low) -> Orange (high)
  ## For z-scores: negative = green, positive = orange
  ## For -log10(FDR): low significance = green, high significance = orange
  custom_gradient = c("#196B23", "#408F1C", "#8EB21E", "#D9B11C", "#E97132"),

  ## Display mode: "zscore" is standard for heatmaps
  display_mode = "zscore",

  ## Legacy labels
  lfc_label   = "log2FC (DESeq2)",
  reference   = "CTL=0 reference"
)

## ============================================================
## ENRICHMENT PLOT COLORS
## ============================================================
## Used by: part4_visualizations.R (ORA bar plots, dot plots, fGSEA)

enrichment_colors <- list(
  ## Gradient for significance: Green (low) -> Orange (high)
  ## Higher -log10(FDR) = more significant = orange
  gradient = c("#196B23", "#408F1C", "#8EB21E", "#D9B11C", "#E97132")
)

## ============================================================
## VOLCANO PLOT COLORS
## ============================================================
## Used by: part1_main_analysis.R, part4_visualizations.R

volcano_colors <- list(
  up   = "#D95F02",   ## Orange for up-regulated
  down = "#02C8D3",   ## Cyan for down-regulated
  ns   = "grey70"     ## Grey for non-significant
)

## ============================================================
## QC PLOT COLORS (PCA/MDS/Density)
## ============================================================
## Used by: part1_main_analysis.R, part4_visualizations.R

qc_plot_params <- list(
  ## Palette for DepotSex colors: any viridis option or "manual"
  depot_sex_palette = "manual",
  ## Manual colors (named vector by DepotSex). Used when palette = "manual".
  depot_sex_manual = c(
    "iWAT_F" = "#5B8FF9",
    "iWAT_M" = "#5AD8A6",
    "gWAT_F" = "#5D7092",
    "gWAT_M" = "#F6BD16"
  )
)

## ============================================================
## BARCODE PLOT PARAMETERS
## ============================================================
## Used by: part5_barcode_plots.R

barcode_params <- list(
  ## Should match gsea_params for consistency
  top_n_per_direction = gsea_params$top_n_per_direction,
  fdr_cutoff          = gsea_params$fdr_cutoff
)

## ============================================================
## HELPER FUNCTION: Get tissue-specific thresholds
## ============================================================

get_tissue_thresholds <- function(contrast_name) {
  if (grepl("iWAT", contrast_name, ignore.case = TRUE)) {
    return(list(logFC_cut = iWAT_logFC_cut, fdr_cut = iWAT_fdr_cut))
  } else if (grepl("gWAT", contrast_name, ignore.case = TRUE)) {
    return(list(logFC_cut = gWAT_logFC_cut, fdr_cut = gWAT_fdr_cut))
  } else {
    return(list(logFC_cut = default_logFC_cut, fdr_cut = default_fdr_cut))
  }
}

## ============================================================
## BIOLOGICAL INTERPRETATION PARAMETERS (Part 6)
## ============================================================
## Used by: part6_biological_interpretation.R

interpretation_params <- list(
  ## Cell type signature scoring
  min_markers_for_score = 3,    ## Minimum markers found to calculate score
  signature_method = "mean_logfc",  ## Options: "mean_logfc", "median_logfc", "ssgsea"

  ## Pathway convergence thresholds
  convergence_fdr = 0.05,       ## FDR cutoff for convergent pathways

  ## Leading edge analysis
  min_le_overlap = 2,           ## Minimum overlap for pathway-cell type reporting
  min_le_fraction = 0.05,       ## Minimum fraction of leading edge for visualization

  ## Publication figure parameters
  top_pathways_for_heatmap = 15,
  celltype_logfc_threshold = 0.3  ## Threshold for asterisk annotation
)

## ============================================================
## PRINT CONFIRMATION
## ============================================================

cat("[INFO] Loaded central parameters from parameters.R\n")
cat("       Reproducibility: MASTER_SEED=", MASTER_SEED, "\n", sep = "")
cat("       DEG thresholds: iWAT FDR<", iWAT_fdr_cut, " |logFC|>", iWAT_logFC_cut, "\n", sep = "")
cat("                       gWAT FDR<", gWAT_fdr_cut, " |logFC|>", gWAT_logFC_cut, "\n", sep = "")
cat("       ORA: top_n=", ora_params$top_n, ", simplify=", ora_params$simplify_cutoff, "\n", sep = "")
cat("       GSEA: top_n_per_direction=", gsea_params$top_n_per_direction, ", FDR<", gsea_params$fdr_cutoff, "\n", sep = "")
cat("       Heatmap: display_mode=", heatmap_params$display_mode, ", clip=", heatmap_params$lfc_clip, "\n", sep = "")
cat("       Interpretation: signature_method=", interpretation_params$signature_method, "\n", sep = "")
