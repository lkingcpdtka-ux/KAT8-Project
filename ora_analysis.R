#!/usr/bin/env Rscript

## =========================================================
## ORA-only pipeline (GO:BP + KEGG)
## Self-contained script: reads a DE table and runs ORA
## =========================================================

required_pkgs <- c(
  "dplyr",
  "ggplot2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (length(args) < hit + 1) return(default)
  args[hit + 1]
}

de_table <- get_arg("--de-table")
outdir <- get_arg("--outdir", "ora_results")
contrast_name <- get_arg("--contrast-name", "contrast")
run_tag <- get_arg("--run-tag", format(Sys.time(), "%Y%m%d_%H%M%S"))
logfc_cutoff <- as.numeric(get_arg("--logfc-cutoff", "1"))
fdr_cutoff <- as.numeric(get_arg("--fdr-cutoff", "0.05"))
min_genes <- as.integer(get_arg("--min-genes", "5"))

if (is.null(de_table)) {
  stop("Missing required argument: --de-table")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

tt <- read.csv(de_table, stringsAsFactors = FALSE)
if (!"gene_name" %in% colnames(tt)) {
  if ("X" %in% colnames(tt)) {
    tt$gene_name <- tt$X
  } else {
    stop("DE table must include gene_name column or rownames in column X.")
  }
}

tt <- tt %>%
  mutate(
    logFC = if ("log2FoldChange" %in% colnames(.)) log2FoldChange else logFC,
    adj.P.Val = if ("padj" %in% colnames(.)) padj else adj.P.Val
  )

run_pathway_analysis <- function(gene_list, direction, contrast_name,
                                 run_tag, outdir, fdr_cutoff = 0.05) {
  cat("\n--- ORA: ", direction, " genes (", contrast_name, ") ---\n", sep = "")
  cat("[INFO] Input gene count: ", length(gene_list), " genes\n", sep = "")

  if (length(gene_list) < min_genes) {
    cat("[WARN] Too few genes (", length(gene_list), ") for ORA\n", sep = "")
    return(NULL)
  }

  gene_entrez <- tryCatch({
    bitr(
      unique(gene_list),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Mm.eg.db,
      drop     = TRUE
    )
  }, error = function(e) {
    cat("[WARN] Gene ID conversion failed: ", conditionMessage(e), "\n", sep = "")
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })

  gene_entrez <- gene_entrez[!is.na(gene_entrez$ENTREZID), , drop = FALSE]
  gene_entrez <- gene_entrez[!duplicated(gene_entrez$ENTREZID), , drop = FALSE]

  if (nrow(gene_entrez) < min_genes) {
    cat("[WARN] Too few genes mapped to Entrez IDs\n")
    return(NULL)
  }

  entrez_ids <- gene_entrez$ENTREZID
  results <- list()

  ## GO:BP (Entrez)
  enrich_go <- enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    minGSSize     = 5,
    maxGSSize     = 500
  )
  if (!is.null(enrich_go) && nrow(enrich_go@result) > 0) {
    results$gobp <- simplify(enrich_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
  }

  ## KEGG (Entrez)
  enrich_kegg <- enrichKEGG(
    gene         = entrez_ids,
    organism     = "mmu",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.2,
    minGSSize    = 5,
    maxGSSize    = 500
  )
  if (!is.null(enrich_kegg) && nrow(enrich_kegg@result) > 0) {
    results$kegg <- enrich_kegg
  }

  if (length(results) == 0) {
    cat("[INFO] No significant ORA results\n")
    return(NULL)
  }

  for (db_name in names(results)) {
    enrich_obj <- results[[db_name]]

    sig_results <- enrich_obj@result %>%
      dplyr::filter(p.adjust < fdr_cutoff) %>%
      dplyr::arrange(p.adjust)

    if (nrow(sig_results) == 0) {
      cat("[INFO] No significant pathways in ", db_name, "\n", sep = "")
      next
    }

    sig_results <- sig_results[, !duplicated(colnames(sig_results)), drop = FALSE]

    table_file <- paste0("Pathway_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".csv")
    write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
    cat("[OK] Saved ORA table: ", table_file, "\n", sep = "")

    plot_data <- sig_results %>%
      dplyr::slice_head(n = 15)

    plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
      num <- as.numeric(x[1]); denom <- as.numeric(x[2])
      if (is.na(denom) || denom == 0) return(0)
      num / denom
    })

    plot_data <- plot_data %>%
      dplyr::arrange(GeneRatio_numeric) %>%
      dplyr::mutate(Description = factor(Description, levels = Description))

    fill_scale <- if (direction == "Up") {
      scale_fill_gradient(low = "#FDD49E", high = "#E69F00", name = "Adj.\nP-value")
    } else {
      scale_fill_gradient(low = "#9ECAE1", high = "#0072B2", name = "Adj.\nP-value")
    }

    p_pathway <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
      fill_scale +
      labs(
        title = paste0(toupper(db_name), " ORA (", direction, ")\n", contrast_name),
        x     = "Gene Ratio",
        y     = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.y     = element_text(size = 10, color = "black"),
        axis.text.x     = element_text(size = 10, color = "black", face = "bold"),
        axis.title.x    = element_text(size = 12, face = "bold"),
        legend.title    = element_text(size = 10, face = "bold"),
        legend.text     = element_text(size = 9),
        legend.position = "right",
        panel.border    = element_rect(color = "black", fill = NA, linewidth = 1)
      )

    plot_file <- paste0("Pathway_barplot_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".png")
    ggsave(
      file.path(outdir, "plots", plot_file),
      plot   = p_pathway,
      width  = 10,
      height = max(6, nrow(plot_data) * 0.3),
      dpi    = 300,
      bg     = "white"
    )
    cat("[OK] Saved ORA bar plot: ", plot_file, "\n", sep = "")
  }

  return(results)
}

up_genes <- tt %>%
  dplyr::filter(adj.P.Val < fdr_cutoff, logFC > logfc_cutoff) %>%
  dplyr::pull(gene_name)

down_genes <- tt %>%
  dplyr::filter(adj.P.Val < fdr_cutoff, logFC < -logfc_cutoff) %>%
  dplyr::pull(gene_name)

cat("[INFO] Up-regulated DEGs for ORA: ", length(up_genes), "\n", sep = "")
cat("[INFO] Down-regulated DEGs for ORA: ", length(down_genes), "\n", sep = "")

if (length(up_genes) >= min_genes) {
  run_pathway_analysis(
    gene_list     = up_genes,
    direction     = "Up",
    contrast_name = contrast_name,
    run_tag       = run_tag,
    outdir        = outdir,
    fdr_cutoff    = fdr_cutoff
  )
}

if (length(down_genes) >= min_genes) {
  run_pathway_analysis(
    gene_list     = down_genes,
    direction     = "Down",
    contrast_name = contrast_name,
    run_tag       = run_tag,
    outdir        = outdir,
    fdr_cutoff    = fdr_cutoff
  )
}

cat("\n=== ORA COMPLETE ===\n")
