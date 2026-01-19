#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 3: FGSEA PATHWAY ANALYSIS
## =========================================================
## This script performs Gene Set Enrichment Analysis (fgsea)
## for GO:BP and KEGG pathways using ranked gene lists from Part 1
##
## FIXES APPLIED:
## - Use DESeq2 Wald statistic for ranking
## - Fix msigdbr KEGG subcategory filtering
## - Better error handling
## - Direction-specific pathway bar plots
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "msigdbr")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "fgsea",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "GO.db"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(msigdbr)
})

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  source(file.path(utils_dir, "purge.R"))
  source(file.path(utils_dir, "dedupe.R"))
  use_save_core <- TRUE
} else {
  cat("[WARN] save_core not found; proceeding without it\n")
  use_save_core <- FALSE
}

## 2.5) Find most recent Part 1 run and use same directory ----
cat("\n=== FINDING PART 1 RESULTS ===\n")

## Look for run directories in savepoints
savepoint_dir <- file.path(getwd(), "savepoints")
if (!dir.exists(savepoint_dir)) {
  stop("savepoints directory not found. Please run part1_main_analysis.R first.")
}

## Find most recent run directory
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) == 0) {
  stop("No Part 1 run directories found in savepoints/")
}

## Sort by modification time and get most recent
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]

## Extract run tag from directory name
run_tag <- gsub("^RUN_", "", basename(outdir))

cat("[INFO] Using Part 1 results from: ", outdir, "\n", sep = "")
cat("[INFO] Run tag: ", run_tag, "\n", sep = "")

## Verify directory structure exists
if (!dir.exists(file.path(outdir, "tables"))) {
  stop("tables directory not found in Part 1 run: ", outdir)
}
if (!dir.exists(file.path(outdir, "plots"))) {
  stop("plots directory not found in Part 1 run: ", outdir)
}

## 3) Parameters --------------------------------------------
fdr_cut   <- 0.05

## Sanity check tracker
fgsea_sanity_tracker <- data.frame(
  Contrast = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_Pathways_Tested = integer(),
  N_Sig_Pathways = integer(),
  stringsAsFactors = FALSE
)

## 4) Main logic --------------------------------------------
tryCatch({

  ## 4.1) Load DE tables from Part 1 ------------------------
  cat("\n=== LOADING DE TABLES FROM PART 1 ===\n")

  ## Find all DE tables
  tables_dir <- file.path(outdir, "tables")
  if (!dir.exists(tables_dir)) {
    stop("tables directory not found in Part 1 run: ", outdir)
  }

  de_files <- list.files(tables_dir, pattern = "^DE_tissue_.*\\.csv$", full.names = TRUE)
  if (length(de_files) == 0) {
    stop("No DE tables found in: ", tables_dir)
  }

  cat("[INFO] Found ", length(de_files), " DE tables\n")

  ## 4.2) Build gene sets ------------------------------------
  cat("\n=== BUILDING GENE SETS FOR fgsea ===\n")

  ## GO:BP gene sets
  cat("[INFO] Building GO:BP gene sets from org.Mm.eg.db...\n")
  go_bp_genes <- tryCatch({
    AnnotationDbi::select(
      org.Mm.eg.db,
      keys = keys(org.Mm.eg.db, keytype = "GOALL"),
      columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
      keytype = "GOALL"
    )
  }, error = function(e) {
    cat("[ERROR] Failed to get GO:BP gene sets: ", conditionMessage(e), "\n")
    return(NULL)
  })

  go_bp_list <- NULL
  if (!is.null(go_bp_genes)) {
    ## Filter for BP ontology and remove NA symbols
    go_bp_genes <- go_bp_genes %>%
      dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
      dplyr::select(GOALL, SYMBOL)

    ## Convert to named list (pathway name -> gene vector)
    go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)

    ## Filter for gene set size
    go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= 5 & sapply(go_bp_list, length) <= 500]

    cat("[OK] GO:BP gene sets: ", length(go_bp_list), " pathways\n", sep = "")
  }

  ## KEGG gene sets (FIXED: use msigdbr correctly)
  cat("[INFO] Building KEGG gene sets from msigdbr...\n")
  kegg_list <- NULL
  tryCatch({
    ## Get all C2 (curated gene sets) for mouse
    msigdb_c2 <- msigdbr(species = "Mus musculus", category = "C2")

    ## Filter for KEGG pathways (gs_subcat is empty for KEGG, but gs_name starts with "KEGG_")
    kegg_msigdb <- msigdb_c2 %>%
      dplyr::filter(grepl("^KEGG_", gs_name))

    if (nrow(kegg_msigdb) == 0) {
      cat("[WARN] No KEGG pathways found in msigdbr C2\n")
    } else {
      ## Convert to named list
      kegg_list <- split(kegg_msigdb$gene_symbol, kegg_msigdb$gs_name)

      ## Filter for gene set size
      kegg_list <- kegg_list[sapply(kegg_list, length) >= 5 & sapply(kegg_list, length) <= 500]

      cat("[OK] KEGG gene sets: ", length(kegg_list), " pathways\n", sep = "")
    }
  }, error = function(e) {
    cat("[WARN] msigdbr KEGG failed: ", conditionMessage(e), "\n")
    cat("[INFO] Falling back to org.Mm.eg.db PATH mapping...\n")

    ## Fallback: build from org.Mm.eg.db
    tryCatch({
      kegg_genes <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = keys(org.Mm.eg.db, keytype = "PATH"),
        columns = c("SYMBOL", "PATH"),
        keytype = "PATH"
      )
      kegg_genes <- kegg_genes %>% dplyr::filter(!is.na(SYMBOL), !is.na(PATH))
      kegg_list <- split(kegg_genes$SYMBOL, paste0("mmu", kegg_genes$PATH))

      ## Filter for gene set size
      kegg_list <- kegg_list[sapply(kegg_list, length) >= 5 & sapply(kegg_list, length) <= 500]

      cat("[OK] KEGG gene sets (from org.Mm.eg.db): ", length(kegg_list), " pathways\n", sep = "")
    }, error = function(e2) {
      cat("[ERROR] KEGG fallback also failed: ", conditionMessage(e2), "\n")
    })
  })

  ## 4.3) fgsea analysis function ----------------------------
  run_fgsea_analysis <- function(ranked_genes, contrast_name,
                                 go_bp_list, kegg_list,
                                 run_tag, outdir, fdr_cutoff = 0.05) {

    cat("\n--- fgsea (GSEA): ", contrast_name, " ---\n", sep = "")
    cat("[INFO] Ranked gene list size: ", length(ranked_genes), " genes\n", sep = "")

    if (length(ranked_genes) < 10) {
      cat("[WARN] Too few genes for fgsea\n")
      return(NULL)
    }

    results <- list()

    ## GO:BP fgsea
    if (!is.null(go_bp_list) && length(go_bp_list) > 0) {
      tryCatch({
        cat("[INFO] Running fgsea for GO:BP (", length(go_bp_list), " gene sets)...\n", sep = "")

        fgsea_go <- fgsea(
          pathways = go_bp_list,
          stats    = ranked_genes,
          minSize  = 5,
          maxSize  = 500,
          nPermSimple = 10000
        )

        if (!is.null(fgsea_go) && nrow(fgsea_go) > 0) {
          ## Add GO term descriptions from GO.db
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
          n_sig <- sum(fgsea_go$padj < fdr_cutoff, na.rm = TRUE)
          cat("[OK] fgsea GO:BP: ", nrow(fgsea_go), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")

          ## Track in sanity checker
          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "GO:BP",
              N_Input_Genes = length(ranked_genes),
              N_Pathways_Tested = nrow(fgsea_go),
              N_Sig_Pathways = n_sig,
              stringsAsFactors = FALSE
            )
          )
        }
      }, error = function(e) {
        cat("[WARN] fgsea GO:BP failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "GO:BP",
            N_Input_Genes = length(ranked_genes),
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            stringsAsFactors = FALSE
          )
        )
      })
    } else {
      cat("[WARN] GO:BP gene sets not available\n")
    }

    ## KEGG fgsea
    if (!is.null(kegg_list) && length(kegg_list) > 0) {
      tryCatch({
        cat("[INFO] Running fgsea for KEGG (", length(kegg_list), " gene sets)...\n", sep = "")

        fgsea_kegg <- fgsea(
          pathways = kegg_list,
          stats    = ranked_genes,
          minSize  = 5,
          maxSize  = 500,
          nPermSimple = 10000
        )

        if (!is.null(fgsea_kegg) && nrow(fgsea_kegg) > 0) {
          ## Add pathway descriptions
          fgsea_kegg <- fgsea_kegg %>%
            dplyr::mutate(
              Description = gsub("^KEGG_|^mmu", "", pathway),  ## Clean up pathway names
              Description = gsub("_", " ", Description)
            ) %>%
            dplyr::arrange(padj)

          results$kegg <- fgsea_kegg
          n_sig <- sum(fgsea_kegg$padj < fdr_cutoff, na.rm = TRUE)
          cat("[OK] fgsea KEGG: ", nrow(fgsea_kegg), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")

          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "KEGG",
              N_Input_Genes = length(ranked_genes),
              N_Pathways_Tested = nrow(fgsea_kegg),
              N_Sig_Pathways = n_sig,
              stringsAsFactors = FALSE
            )
          )
        }
      }, error = function(e) {
        cat("[WARN] fgsea KEGG failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "KEGG",
            N_Input_Genes = length(ranked_genes),
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            stringsAsFactors = FALSE
          )
        )
      })
    } else {
      cat("[WARN] KEGG gene sets not available\n")
    }

    ## Save results with direction-specific colors
    if (length(results) > 0) {
      for (db_name in names(results)) {
        fgsea_obj <- results[[db_name]]

        sig_results <- fgsea_obj %>%
          dplyr::filter(padj < fdr_cutoff) %>%
          dplyr::arrange(padj)

        if (nrow(sig_results) == 0) {
          cat("[INFO] No significant fgsea pathways in ", db_name, "\n", sep = "")
          next
        }

        ## Save table
        table_file <- paste0("fgsea_", db_name, "_", contrast_name, "_", run_tag, ".csv")
        write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
        cat("[OK] Saved fgsea table: ", table_file, "\n", sep = "")

        ## Create plots - separate by direction (NES > 0 = Up, NES < 0 = Down)
        plot_data_up <- sig_results %>%
          dplyr::filter(NES > 0) %>%
          dplyr::slice_head(n = 10) %>%
          dplyr::arrange(NES)

        plot_data_down <- sig_results %>%
          dplyr::filter(NES < 0) %>%
          dplyr::slice_head(n = 10) %>%
          dplyr::arrange(desc(NES))

        ## Plot UP-regulated pathways
        if (nrow(plot_data_up) > 0) {
          plot_data_up <- plot_data_up %>%
            dplyr::mutate(
              pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
              pathway_label = factor(pathway_label, levels = pathway_label)
            )

          p_fgsea_up <- ggplot(plot_data_up, aes(x = NES, y = pathway_label, fill = padj)) +
            geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
            scale_fill_gradient(
              low  = "#FDD49E",  ## Light orange
              high = "#E69F00",  ## Dark orange
              name = "Adj.\nP-value"
            ) +
            labs(
              title = paste0("fgsea ", toupper(db_name), " (Up-regulated): ", contrast_name),
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

        ## Plot DOWN-regulated pathways
        if (nrow(plot_data_down) > 0) {
          plot_data_down <- plot_data_down %>%
            dplyr::mutate(
              pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
              pathway_label = factor(pathway_label, levels = pathway_label)
            )

          p_fgsea_down <- ggplot(plot_data_down, aes(x = NES, y = pathway_label, fill = padj)) +
            geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
            scale_fill_gradient(
              low  = "#9ECAE1",  ## Light blue
              high = "#0072B2",  ## Dark blue/teal
              name = "Adj.\nP-value"
            ) +
            labs(
              title = paste0("fgsea ", toupper(db_name), " (Down-regulated): ", contrast_name),
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
    }

    return(results)
  }

  ## 4.4) Process each DE table ------------------------------
  cat("\n=== RUNNING fgsea FOR ALL CONTRASTS ===\n")

  for (de_file in de_files) {

    ## Extract contrast name from filename
    de_filename <- basename(de_file)
    contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)

    cat("\n=== Processing contrast: ", contrast_name, " ===\n", sep = "")

    ## Load DE table
    de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

    ## Ensure we have the stat column (Wald statistic)
    if (!"stat" %in% colnames(de_table)) {
      cat("[WARN] Cannot find 'stat' column (Wald statistic) in ", de_file, "\n")
      next
    }

    ## Create ranked gene list using Wald statistic
    ranked_genes <- de_table %>%
      dplyr::filter(is.finite(stat)) %>%
      dplyr::arrange(dplyr::desc(stat))

    ## Create named vector (gene names -> Wald statistic)
    ranked_vec <- setNames(ranked_genes$stat, rownames(ranked_genes))

    cat("[INFO] Ranked gene list for fgsea: ", length(ranked_vec), " genes\n", sep = "")

    ## Run fgsea
    if (length(ranked_vec) >= 10) {
      run_fgsea_analysis(
        ranked_genes  = ranked_vec,
        contrast_name = contrast_name,
        go_bp_list    = go_bp_list,
        kegg_list     = kegg_list,
        run_tag       = run_tag,
        outdir        = outdir,
        fdr_cutoff    = fdr_cut
      )
    } else {
      cat("[WARN] Too few genes for fgsea\n")
    }
  }

  ## 4.5) Save sanity check table ----------------------------
  sanity_file <- paste0("fgsea_sanity_check_", run_tag, ".csv")
  write.csv(fgsea_sanity_tracker, file = file.path(outdir, "tables", sanity_file), row.names = FALSE)
  cat("\n[OK] Saved fgsea sanity check table: ", sanity_file, "\n", sep = "")

  cat("\n======================================\n")
  cat("=== PART 3 (fgsea) COMPLETE ===\n")
  cat("======================================\n\n")

}, error = function(e) {

  cat("\n======================================\n")
  cat("=== ERROR IN fgsea PIPELINE ===\n")
  cat("======================================\n")
  cat(conditionMessage(e), "\n\n")

  err_file <- paste0("ERROR_fgsea_", run_tag, ".txt")
  writeLines(
    c(
      "fgsea script error:",
      conditionMessage(e),
      "",
      "Traceback:",
      capture.output(traceback())
    ),
    con = file.path(outdir, "logs", err_file)
  )

  stop(e)
})

cat("\n=== Part 3 (fgsea) finished ===\n")
