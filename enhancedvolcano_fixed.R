## =============================================================================
## CORRECTED EnhancedVolcano Section (only color significant genes by direction)
## =============================================================================

## This replaces the EnhancedVolcano section in your per-contrast loop

cat("\n--- EnhancedVolcano: Genes of Interest (", cn, ") ---\n", sep = "")

## Prepare data for EnhancedVolcano
volcano_df <- tt %>%
  dplyr::filter(is.finite(logFC), is.finite(adj.P.Val), adj.P.Val > 0)

## Only label genes of interest that are SIGNIFICANT (pass both thresholds)
goi_significant <- volcano_df %>%
  dplyr::filter(gene_name %in% genes_of_interest,
                adj.P.Val < fdr_cut_tissue,
                abs(logFC) >= logFC_cut_tissue) %>%
  dplyr::pull(gene_name)

cat("[INFO] Genes of interest (total): ", length(genes_of_interest), "\n", sep = "")
cat("[INFO] Genes of interest (significant): ", length(goi_significant), "\n", sep = "")

if (nrow(volcano_df) > 0) {
  
  ## FIXED: Create custom color scheme using keyvals
  ## Only color significant genes (both thresholds met), split by direction
  ## All others (NS, FC-only, P-only) are grey
  
  keyvals_color <- rep("grey70", nrow(volcano_df))
  names(keyvals_color) <- rep("NS", nrow(volcano_df))
  
  ## Significant UP (orange)
  sig_up_idx <- which(volcano_df$adj.P.Val < fdr_cut_tissue & 
                      volcano_df$logFC >= logFC_cut_tissue)
  keyvals_color[sig_up_idx] <- "#E69F00"
  names(keyvals_color)[sig_up_idx] <- "Significant Up"
  
  ## Significant DOWN (teal/blue)
  sig_down_idx <- which(volcano_df$adj.P.Val < fdr_cut_tissue & 
                        volcano_df$logFC <= -logFC_cut_tissue)
  keyvals_color[sig_down_idx] <- "#0072B2"
  names(keyvals_color)[sig_down_idx] <- "Significant Down"
  
  ev_plot <- EnhancedVolcano(
    volcano_df,
    lab = volcano_df$gene_name,
    selectLab = goi_significant,
    x = "logFC",
    y = "adj.P.Val",
    title = paste0("Volcano: ", cn),
    subtitle = "ECM, metabolism, and remodeling genes labeled (significant only)",
    pCutoff = fdr_cut_tissue,
    FCcutoff = logFC_cut_tissue,
    pointSize = 2.0,
    labSize = 4.0,
    labCol = "black",
    labFace = "bold",
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = "gray40",
    colCustom = keyvals_color,  ## Use custom colors
    colAlpha = 0.8,
    legendPosition = "right",
    legendLabSize = 10,
    legendIconSize = 3.0,
    max.overlaps = Inf,
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  )
  
  ev_file <- paste0("EnhancedVolcano_GOI_", cn, "_", run_tag, ".png")
  ggsave(file.path(outdir, "plots", ev_file), plot = ev_plot, width = 10, height = 8, dpi = 300)
  cat("[OK] EnhancedVolcano saved: ", ev_file, "\n", sep = "")
  
} else {
  cat("[WARN] No valid data for EnhancedVolcano\n")
}
