## =========================================================
## KAT8 bulk RNA-seq - CENTRAL PARAMETERS FILE
## =========================================================
## All analysis cutoffs and thresholds are defined HERE ONLY
## Other scripts (part1-5) source this file to use these values
## =========================================================

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

  ## GO:BP simplification
  simplify_go     = TRUE,    ## Always simplify GO:BP to remove redundant terms
  simplify_cutoff = 0.7,

  ## Plotting parameters
  top_n_per_direction = 10,  ## Top pathways per direction for plots

  ## FDR cutoff for significance
  fdr_cutoff      = 0.05,

  ## Random seed for reproducibility (used by gseGO/gseKEGG)
  seed            = 12345
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

  ## GO:BP simplification cutoff (should match ora_params)
  gobp_simplify_cutoff = 0.7
)

## ============================================================
## HEATMAP PARAMETERS
## ============================================================
## Used by: part4_publication_plots.R

heatmap_params <- list(
  lfc_clip    = 2.0,         ## Clip log2FC at +/- this value
  lfc_label   = "log2FC (DESeq2)",
  reference   = "CTL=0 reference"
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
## PRINT CONFIRMATION
## ============================================================

cat("[INFO] Loaded central parameters from parameters.R\n")
cat("       DEG thresholds: iWAT FDR<", iWAT_fdr_cut, " |logFC|>", iWAT_logFC_cut, "\n", sep = "")
cat("                       gWAT FDR<", gWAT_fdr_cut, " |logFC|>", gWAT_logFC_cut, "\n", sep = "")
cat("       ORA: top_n=", ora_params$top_n, ", simplify=", ora_params$simplify_cutoff, "\n", sep = "")
cat("       GSEA: top_n_per_direction=", gsea_params$top_n_per_direction, ", FDR<", gsea_params$fdr_cutoff, "\n", sep = "")
