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
## GENES OF INTEREST  <<< ADD YOUR GENES HERE, ONCE >>>
## ============================================================
## One list, honoured by every plot that labels individual genes:
##   part1  - genes-of-interest heatmap + EnhancedVolcano
##   part4  - volcano plots (always labelled if they pass DE thresholds)
##   part4b - cell-vs-tissue concordance quadrant plot
##
## These do NOT need to be transcription factors. The only TF-restricted
## analyses are part5b layers 1 and 3, which are limited to TFs by definition
## (GO:0003700 / TF target sets) -- everything else labels any gene you list.
##
## Seeded with KAT8 itself and its MSL / NSL complex partners: their mRNA is
## expected to stay FLAT under siRNA (which depletes KAT8 protein, not partner
## transcripts), so they act as a built-in negative control on the plots.
GENES_OF_INTEREST <- c(
  "Kat8",                                         ## catalytic subunit
  "Msl1", "Msl2", "Msl3",                         ## MSL complex
  "Kansl1", "Kansl2", "Kansl3", "Mcrs1", "Phf20", ## NSL/KANSL complex
  ## Medium-chain acyl-CoA synthetases: among the strongest tissue effects
  ## (Acsm3 iWAT log2FC -5.8, gWAT -3.1) and NOT expressed in 3T3-L1, so they
  ## are tissue-only and cannot be assessed for cell-autonomy.
  "Acsm3", "Acsm5"
)

## ============================================================
## CELL (3T3-L1) PARAMETERS
## ============================================================
## Used by: KAT8_bulk_cells_comprehensive.R, downstream_analysis/*
## Previously these lived INSIDE the cells script, which contradicted this
## file's "all cutoffs defined HERE ONLY" contract. Values are unchanged.
##
## IMPORTANT NOTE ON THE logFC GATE:
##   The cells are a MAINTENANCE experiment (KAT8 siRNA in already-
##   differentiated 3T3-L1). Effects are real but small: among FDR<0.05 genes
##   the median |log2FC| is ~0.37, so gating at 0.5 discards ~74% of genuine
##   DEGs (290 of 392) before any downstream analysis. FDR-only is therefore
##   the PRIMARY DEG definition for cells; the logFC-gated list is retained
##   for backward compatibility and reporting.
cells_params <- list(
  ## Prefilter: counts >= min_count in >= min_samples_per_group in AT LEAST
  ## ONE group (appropriate for a 2-group comparison).
  prefilter_min_count            = 10,
  prefilter_min_samples_per_group = 2,

  ## DE thresholds
  fdr_cutoff   = 0.05,
  logFC_cutoff = 0.5,     ## descriptive gate, NOT the primary DEG rule

  ## Which rule defines the AUTHORITATIVE cell DEG list used downstream:
  ##   "fdr_only"    - FDR < fdr_cutoff                (RECOMMENDED)
  ##   "fdr_and_lfc" - FDR < fdr_cutoff AND |logFC| >  (legacy behaviour)
  primary_deg_rule = "fdr_only"
)

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
  ## WARNING: simplify() is O(n�) complexity - can hang with many terms!
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

## Fraction of EXTRA labels reserved for the largest fold changes, on top of
## the composite-ranked ones. Guarantees the volcano labels both the most
## significant genes AND the biggest effects, instead of clustering them all
## at the top. 0 disables (composite ranking only).
volcano_fc_label_fraction <- 0.5

volcano_gene_selection <- list(
  ## >>> CHANGE THIS to control how many genes are labeled per direction <<<
  top_n_per_direction = 20,
  
  ## ORA pathway significance threshold for counting gene membership
  ora_padj_cutoff     = 0.05,
  
  ## Minimum pathway count for a gene to be considered ORA-supported
  min_pathway_count   = 1,
  
  ## Genes always labeled if they pass DE thresholds (biological anchors).
  ## Inherits the single GENES_OF_INTEREST list defined above.
  mandatory_genes     = GENES_OF_INTEREST,
  
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
  ## --- PCA construction (these change how the plot LOOKS) ---
  ## pca_ntop  : number of most-variable genes to use. 500 = DESeq2::plotPCA
  ##             default. Set to Inf to use every gene (the original setting).
  ## pca_scale : unit-scale each gene before PCA? FALSE = DESeq2 default.
  ##             TRUE gives every gene equal weight, so low-variance noise
  ##             contributes heavily and PC1/PC2 explain much less variance.
  ## The sign of a principal component is arbitrary -- a mirrored plot between
  ## runs is not an error.
  pca_ntop  = 500,
  pca_scale = FALSE,

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