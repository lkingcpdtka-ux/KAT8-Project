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
## MODEL DESIGN (Unified DESeq2 Model)
## ============================================================
## Used by: part1_main_analysis.R
##
## DESIGN: ~ Sex + GroupDepot
##   - GroupDepot is a 4-level factor: iWAT_CTL, iWAT_KAT8KD, gWAT_CTL, gWAT_KAT8KD
##   - Sex as covariate absorbs sex-driven variance (not the focus)
##   - ALL 38 tissue samples in ONE model (shared dispersion estimates, more power)
##   - Per-depot KAT8 effects extracted via contrasts:
##       iWAT: contrast = c("GroupDepot", "iWAT_KAT8KD", "iWAT_CTL")
##       gWAT: contrast = c("GroupDepot", "gWAT_KAT8KD", "gWAT_CTL")
##
## JUSTIFICATION (Part 7.7 validation):
##   - Sex x Genotype interaction is minimal (<4% of genes, pi0 > 0.85)
##   - Sex-stratified logFC concordance >99% directional agreement
##   - Adding Sex as covariate gains +5-14% DEGs with >99% concordance
##
## WHY NOT ~ Sex * Depot * Genotype?
##   - Too many parameters (7+ coefficients from 38 samples)
##   - Underpowered for three-way interaction
##   - Cell-means (GroupDepot) approach is simpler and equally flexible

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
  plot_width  = 6.9,
  plot_height = 6.8,
  label_wrap_width = 32
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
  plot_width  = 6.9,
  plot_height = 6.8,
  label_wrap_width = 32
)

## ============================================================
## HEATMAP PARAMETERS
## ============================================================
## Used by: part4_visualizations.R

heatmap_params <- list(
  lfc_clip    = 2.0,         ## Clip z-scores at +/- this value

  ## Custom gradient: RdYlBu (ColorBrewer) diverging
  ## Down (blue) -> neutral (yellow) -> up (red)
  custom_gradient = c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027"),

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
  ## Gradient: RdYlBu diverging
  gradient = c("#4575B4", "#91BFDB", "#FFFFBF", "#FC8D59", "#D73027")
)

## ============================================================
## VOLCANO PLOT GENE SELECTION
## ============================================================
## Used by: part4_visualizations.R (ORA-informed labeling)
## Genes are selected by intersecting significant ORA pathways
## with DEGs, then ranked by a composite score:
##   score = z(pathway_count) + z(|logFC|) + z(-log10(FDR))

volcano_gene_selection <- list(
  ## >>> CHANGE THIS to control how many genes are labeled per direction <<<
  top_n_per_direction = 20,

  ## ORA pathway significance threshold for counting gene membership
  ora_padj_cutoff     = 0.05,

  ## Minimum pathway count for a gene to be considered ORA-supported
  min_pathway_count   = 1,

  ## Genes always labeled if they pass DE thresholds (biological anchors)
  mandatory_genes     = c("Kat8"),

  ## Fallback: when no ORA pathways exist for a direction, use DE-only ranking
  ## (top genes by padj, then by |logFC|)
  fallback_to_de      = TRUE
)

## ============================================================
## VOLCANO PLOT COLORS
## ============================================================
## Used by: part1_main_analysis.R, part4_visualizations.R

volcano_colors <- list(
  up   = "#D73027",   ## Red for up-regulated
  down = "#4575B4",   ## Blue for down-regulated
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
    "iWAT_F" = "#4575B4",  ## blue
    "iWAT_M" = "#91BFDB",  ## light blue
    "gWAT_F" = "#FC8D59",  ## orange
    "gWAT_M" = "#D73027"   ## red
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
## DEPOT CONCORDANCE PARAMETERS (Part 7)
## ============================================================
## Used by: part7_depot_concordance.R
## Compares iWAT vs gWAT KAT8-KD responses to identify
## shared (depot-general) vs depot-specific DEGs

depot_concordance_params <- list(
  ## Correlation method for logFC comparison
  cor_method = "pearson",

  ## Number of top genes to label on scatter plot
  top_n_label = 20,

  ## Plot dimensions
  plot_width  = 10,
  plot_height = 8
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
