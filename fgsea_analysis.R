#!/usr/bin/env Rscript

## =========================================================
## fgsea-only pipeline (GO:BP + KEGG)
## Self-contained script: reads a DE table and runs fgsea
## =========================================================

required_pkgs <- c(
  "dplyr",
  "ggplot2",
  "fgsea",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "GO.db",
  "msigdbr"
)

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(GO.db)
})

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (length(args) < hit + 1) return(default)
  args[hit + 1]
}

de_table <- get_arg("--de-table")
outdir <- get_arg("--outdir", "fgsea_results")
contrast_name <- get_arg("--contrast-name", "contrast")
run_tag <- get_arg("--run-tag", format(Sys.time(), "%Y%m%d_%H%M%S"))
fdr_cutoff <- as.numeric(get_arg("--fdr-cutoff", "0.05"))
min_size <- as.integer(get_arg("--min-size", "5"))
max_size <- as.integer(get_arg("--max-size", "500"))
n_perm <- as.integer(get_arg("--n-perm", "10000"))

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

ranked_genes <- NULL
if ("stat" %in% colnames(tt)) {
  ranked_genes <- tt %>%
    dplyr::filter(is.finite(stat)) %>%
    dplyr::mutate(rank_stat = stat) %>%
    dplyr::arrange(dplyr::desc(rank_stat))
} else if (all(c("log2FoldChange", "pvalue") %in% colnames(tt))) {
  ranked_genes <- tt %>%
    dplyr::filter(is.finite(log2FoldChange), is.finite(pvalue), pvalue > 0) %>%
    dplyr::mutate(rank_stat = sign(log2FoldChange) * -log10(pvalue)) %>%
    dplyr::arrange(dplyr::desc(rank_stat))
  cat("[WARN] Using sign(logFC) * -log10(pvalue) ranking (stat not found).\n")
} else {
  stop("DE table must include stat or (log2FoldChange + pvalue).")
}

ranked_vec <- setNames(ranked_genes$rank_stat, ranked_genes$gene_name)
cat("[INFO] Ranked gene list for fgsea: ", length(ranked_vec), " genes\n", sep = "")

if (length(ranked_vec) < 10) {
  stop("Too few genes for fgsea.")
}

run_fgsea <- function(ranked_vec, contrast_name, run_tag, outdir, fdr_cutoff) {
  results <- list()

  ## GO:BP gene sets
  cat("[INFO] Running fgsea for GO:BP...\n")
  go_bp_genes <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = keys(org.Mm.eg.db, keytype = "GOALL"),
    columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
    keytype = "GOALL"
  )

  go_bp_genes <- go_bp_genes %>%
    dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
    dplyr::select(GOALL, SYMBOL)

  go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)
  go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= min_size & sapply(go_bp_list, length) <= max_size]

  fgsea_go <- fgsea(
    pathways = go_bp_list,
    stats    = ranked_vec,
    minSize  = min_size,
    maxSize  = max_size,
    nPermSimple = n_perm
  )

  if (!is.null(fgsea_go) && nrow(fgsea_go) > 0) {
    go_terms <- AnnotationDbi::select(
      GO.db,
      keys = fgsea_go$pathway,
      columns = "TERM",
      keytype = "GOID"
    )

    fgsea_go <- fgsea_go %>%
      dplyr::left_join(
        go_terms %>% dplyr::rename(pathway = GOID, Description = TERM),
        by = "pathway"
      ) %>%
      dplyr::arrange(padj)

    results$gobp <- fgsea_go
  }

  ## KEGG gene sets
  cat("[INFO] Running fgsea for KEGG...\n")
  kegg_list <- NULL
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    kegg_msigdb <- suppressMessages(msigdbr::msigdbr(species = "Mus musculus", category = "C2"))
    kegg_msigdb <- kegg_msigdb %>% dplyr::filter(gs_subcat == "CP:KEGG")
    if (nrow(kegg_msigdb) > 0) {
      kegg_list <- split(kegg_msigdb$gene_symbol, kegg_msigdb$gs_name)
    }
  }
  if (is.null(kegg_list)) {
    kegg_genes <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys = keys(org.Mm.eg.db, keytype = "PATH"),
      columns = c("SYMBOL", "PATH"),
      keytype = "PATH"
    )
    kegg_genes <- kegg_genes %>% dplyr::filter(!is.na(SYMBOL), !is.na(PATH))
    kegg_list <- split(kegg_genes$SYMBOL, paste0("mmu", kegg_genes$PATH))
  }

  kegg_list <- kegg_list[sapply(kegg_list, length) >= min_size & sapply(kegg_list, length) <= max_size]

  fgsea_kegg <- fgsea(
    pathways = kegg_list,
    stats    = ranked_vec,
    minSize  = min_size,
    maxSize  = max_size,
    nPermSimple = n_perm
  )

  if (!is.null(fgsea_kegg) && nrow(fgsea_kegg) > 0) {
    fgsea_kegg <- fgsea_kegg %>%
      dplyr::mutate(
        Description = gsub("^KEGG_|^mmu", "", pathway),
        Description = gsub("_", " ", Description)
      ) %>%
      dplyr::arrange(padj)

    results$kegg <- fgsea_kegg
  }

  if (length(results) == 0) {
    cat("[INFO] No significant fgsea results\n")
    return(NULL)
  }

  for (db_name in names(results)) {
    fgsea_obj <- results[[db_name]]

    sig_results <- fgsea_obj %>%
      dplyr::filter(padj < fdr_cutoff) %>%
      dplyr::arrange(padj)

    if (nrow(sig_results) == 0) {
      cat("[INFO] No significant fgsea pathways in ", db_name, "\n", sep = "")
      next
    }

    table_file <- paste0("fgsea_", db_name, "_", contrast_name, "_", run_tag, ".csv")
    write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
    cat("[OK] Saved fgsea table: ", table_file, "\n", sep = "")

    plot_data_up <- sig_results %>%
      dplyr::filter(NES > 0) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::arrange(NES)

    plot_data_down <- sig_results %>%
      dplyr::filter(NES < 0) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::arrange(desc(NES))

    if (nrow(plot_data_up) > 0) {
      plot_data_up <- plot_data_up %>%
        dplyr::mutate(
          pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
          pathway_label = factor(pathway_label, levels = pathway_label)
        )

      p_fgsea_up <- ggplot(plot_data_up, aes(x = NES, y = pathway_label, fill = padj)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
        scale_fill_gradient(
          low  = "#FDD49E",
          high = "#E69F00",
          name = "Adj.\nP-value"
        ) +
        labs(
          title = paste0("fgsea ", toupper(db_name), " (Up): ", contrast_name),
          x = "Normalized Enrichment Score (NES)",
          y = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black", face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        )

      plot_file_up <- paste0("fgsea_plot_", db_name, "_", contrast_name, "_Up_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", plot_file_up), plot = p_fgsea_up, width = 12,
             height = max(6, nrow(plot_data_up) * 0.4), dpi = 300, bg = "white")
      cat("[OK] Saved fgsea UP plot: ", plot_file_up, "\n", sep = "")
    }

    if (nrow(plot_data_down) > 0) {
      plot_data_down <- plot_data_down %>%
        dplyr::mutate(
          pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
          pathway_label = factor(pathway_label, levels = pathway_label)
        )

      p_fgsea_down <- ggplot(plot_data_down, aes(x = NES, y = pathway_label, fill = padj)) +
        geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
        scale_fill_gradient(
          low  = "#9ECAE1",
          high = "#0072B2",
          name = "Adj.\nP-value"
        ) +
        labs(
          title = paste0("fgsea ", toupper(db_name), " (Down): ", contrast_name),
          x = "Normalized Enrichment Score (NES)",
          y = NULL
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black", face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        )

      plot_file_down <- paste0("fgsea_plot_", db_name, "_", contrast_name, "_Down_", run_tag, ".png")
      ggsave(file.path(outdir, "plots", plot_file_down), plot = p_fgsea_down, width = 12,
             height = max(6, nrow(plot_data_down) * 0.4), dpi = 300, bg = "white")
      cat("[OK] Saved fgsea DOWN plot: ", plot_file_down, "\n", sep = "")
    }
  }

  return(results)
}

run_fgsea(
  ranked_vec = ranked_vec,
  contrast_name = contrast_name,
  run_tag = run_tag,
  outdir = outdir,
  fdr_cutoff = fdr_cutoff
)

cat("\n=== FGSEA COMPLETE ===\n")
