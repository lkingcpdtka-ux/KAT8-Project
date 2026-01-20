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
  "GO.db",
  "clusterProfiler",
  "GOSemSim"
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
  library(clusterProfiler)
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
simplify_go_bp <- FALSE
simplify_go_cutoff <- 0.7
top_n_per_direction <- 10

## Library toggles (temporary)
run_go_bp <- FALSE
run_kegg <- FALSE
run_wikipathways <- TRUE
run_hallmark <- TRUE

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## Sanity check tracker
fgsea_sanity_tracker <- data.frame(
  Contrast = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_Pathways_Tested = integer(),
  N_Sig_Pathways = integer(),
  N_Sig_Up = integer(),
  N_Sig_Down = integer(),
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
  log_time("Starting gene set construction")
  
  ## GO:BP gene sets (disabled for now)
  go_bp_list <- NULL
  go_bp_term2gene <- NULL
  if (run_go_bp) {
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
    
    if (!is.null(go_bp_genes)) {
      ## Filter for BP ontology and remove NA symbols
      go_bp_genes <- go_bp_genes %>%
        dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
        dplyr::select(GOALL, SYMBOL)
      
      ## Convert to named list (pathway name -> gene vector)
      go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)
      go_bp_term2gene <- go_bp_genes %>% dplyr::select(GOALL, SYMBOL)
      
      ## Filter for gene set size
      go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= 5 & sapply(go_bp_list, length) <= 500]
      
      cat("[OK] GO:BP gene sets: ", length(go_bp_list), " pathways\n", sep = "")
    }
  } else {
    cat("[INFO] GO:BP gene sets disabled (run_go_bp = FALSE)\n")
  }
  
  ## KEGG gene sets (disabled for now)
  kegg_list <- NULL
  if (run_kegg) {
    cat("[INFO] Building KEGG gene sets from msigdbr...\n")
    tryCatch({
      ## Get all C2 (curated gene sets) for mouse
      msigdb_c2 <- msigdbr(species = "Mus musculus", collection = "C2")
      
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
  } else {
    cat("[INFO] KEGG gene sets disabled (run_kegg = FALSE)\n")
  }
  
  ## WikiPathways gene sets
  wikipathways_list <- NULL
  if (run_wikipathways) {
    cat("[INFO] Building WikiPathways gene sets from msigdbr...\n")
    tryCatch({
      msigdb_wp <- msigdbr(
        species = "Mus musculus",
        collection = "C2",
        subcollection = "CP:WIKIPATHWAYS"
      )
      if (nrow(msigdb_wp) == 0) {
        cat("[WARN] No WikiPathways gene sets found in msigdbr\n")
      } else {
        wikipathways_list <- split(msigdb_wp$gene_symbol, msigdb_wp$gs_name)
        wikipathways_list <- wikipathways_list[
          sapply(wikipathways_list, length) >= 5 & sapply(wikipathways_list, length) <= 500
        ]
        cat("[OK] WikiPathways gene sets: ", length(wikipathways_list), " pathways\n", sep = "")
      }
    }, error = function(e) {
      cat("[WARN] msigdbr WikiPathways failed: ", conditionMessage(e), "\n")
    })
  } else {
    cat("[INFO] WikiPathways gene sets disabled (run_wikipathways = FALSE)\n")
  }
  
  ## Hallmark gene sets
  hallmark_list <- NULL
  if (run_hallmark) {
    cat("[INFO] Building Hallmark gene sets from msigdbr...\n")
    tryCatch({
      msigdb_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
      if (nrow(msigdb_hallmark) == 0) {
        cat("[WARN] No Hallmark gene sets found in msigdbr\n")
      } else {
        hallmark_list <- split(msigdb_hallmark$gene_symbol, msigdb_hallmark$gs_name)
        hallmark_list <- hallmark_list[sapply(hallmark_list, length) >= 5 & sapply(hallmark_list, length) <= 500]
        cat("[OK] Hallmark gene sets: ", length(hallmark_list), " pathways\n", sep = "")
      }
    }, error = function(e) {
      cat("[WARN] msigdbr Hallmark failed: ", conditionMessage(e), "\n")
    })
  } else {
    cat("[INFO] Hallmark gene sets disabled (run_hallmark = FALSE)\n")
  }
  log_time("Completed gene set construction")
  
  ## 4.3) fgsea analysis function ----------------------------
  run_fgsea_analysis <- function(ranked_genes, contrast_name,
                                 go_bp_list, go_bp_term2gene, kegg_list, wikipathways_list, hallmark_list,
                                 run_tag, outdir, fdr_cutoff = 0.05,
                                 simplify_go = FALSE, simplify_cutoff = 0.7,
                                 top_n_per_direction = 10) {
    
    cat("\n--- fgsea (GSEA): ", contrast_name, " ---\n", sep = "")
    cat("[INFO] Ranked gene list size: ", length(ranked_genes), " genes\n", sep = "")
    
    if (length(ranked_genes) < 10) {
      cat("[WARN] Too few genes for fgsea\n")
      return(NULL)
    }
    
    results <- list()
    
    ## GO:BP fgsea
    if (run_go_bp && !is.null(go_bp_list) && length(go_bp_list) > 0) {
      log_time(paste0("Starting GO:BP fgsea for ", contrast_name))
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
          n_sig_up <- sum(fgsea_go$padj < fdr_cutoff & fgsea_go$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(fgsea_go$padj < fdr_cutoff & fgsea_go$NES < 0, na.rm = TRUE)
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
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
              stringsAsFactors = FALSE
            )
          )
          
          if (simplify_go) {
            cat("[INFO] GO:BP simplification skipped (fgsea-only logic enabled)\n")
          }
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
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed GO:BP fgsea for ", contrast_name))
    } else if (run_go_bp) {
      cat("[WARN] GO:BP gene sets not available\n")
    }
    
    ## KEGG fgsea
    if (run_kegg && !is.null(kegg_list) && length(kegg_list) > 0) {
      log_time(paste0("Starting KEGG fgsea for ", contrast_name))
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
          n_sig_up <- sum(fgsea_kegg$padj < fdr_cutoff & fgsea_kegg$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(fgsea_kegg$padj < fdr_cutoff & fgsea_kegg$NES < 0, na.rm = TRUE)
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
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
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
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed KEGG fgsea for ", contrast_name))
    } else if (run_kegg) {
      cat("[WARN] KEGG gene sets not available\n")
    }
    
    ## WikiPathways fgsea
    if (run_wikipathways && !is.null(wikipathways_list) && length(wikipathways_list) > 0) {
      log_time(paste0("Starting WikiPathways fgsea for ", contrast_name))
      tryCatch({
        cat("[INFO] Running fgsea for WikiPathways (", length(wikipathways_list), " gene sets)...\n", sep = "")
        
        fgsea_wp <- fgsea(
          pathways = wikipathways_list,
          stats    = ranked_genes,
          minSize  = 5,
          maxSize  = 500,
          nPermSimple = 10000
        )
        
        if (!is.null(fgsea_wp) && nrow(fgsea_wp) > 0) {
          fgsea_wp <- fgsea_wp %>%
            dplyr::mutate(
              Description = gsub("^WIKIPATHWAYS_|^WP", "", pathway),
              Description = gsub("_", " ", Description)
            ) %>%
            dplyr::arrange(padj)
          
          results$wikipathways <- fgsea_wp
          n_sig <- sum(fgsea_wp$padj < fdr_cutoff, na.rm = TRUE)
          n_sig_up <- sum(fgsea_wp$padj < fdr_cutoff & fgsea_wp$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(fgsea_wp$padj < fdr_cutoff & fgsea_wp$NES < 0, na.rm = TRUE)
          cat("[OK] fgsea WikiPathways: ", nrow(fgsea_wp), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")
          
          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "WikiPathways",
              N_Input_Genes = length(ranked_genes),
              N_Pathways_Tested = nrow(fgsea_wp),
              N_Sig_Pathways = n_sig,
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
              stringsAsFactors = FALSE
            )
          )
        }
      }, error = function(e) {
        cat("[WARN] fgsea WikiPathways failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "WikiPathways",
            N_Input_Genes = length(ranked_genes),
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed WikiPathways fgsea for ", contrast_name))
    } else if (run_wikipathways) {
      cat("[WARN] WikiPathways gene sets not available\n")
    }
    
    ## Hallmark fgsea
    if (run_hallmark && !is.null(hallmark_list) && length(hallmark_list) > 0) {
      log_time(paste0("Starting Hallmark fgsea for ", contrast_name))
      tryCatch({
        cat("[INFO] Running fgsea for Hallmark (", length(hallmark_list), " gene sets)...\n", sep = "")
        
        fgsea_hallmark <- fgsea(
          pathways = hallmark_list,
          stats    = ranked_genes,
          minSize  = 5,
          maxSize  = 500,
          nPermSimple = 10000
        )
        
        if (!is.null(fgsea_hallmark) && nrow(fgsea_hallmark) > 0) {
          fgsea_hallmark <- fgsea_hallmark %>%
            dplyr::mutate(
              Description = gsub("^HALLMARK_", "", pathway),
              Description = gsub("_", " ", Description)
            ) %>%
            dplyr::arrange(padj)
          
          results$hallmark <- fgsea_hallmark
          n_sig <- sum(fgsea_hallmark$padj < fdr_cutoff, na.rm = TRUE)
          n_sig_up <- sum(fgsea_hallmark$padj < fdr_cutoff & fgsea_hallmark$NES > 0, na.rm = TRUE)
          n_sig_down <- sum(fgsea_hallmark$padj < fdr_cutoff & fgsea_hallmark$NES < 0, na.rm = TRUE)
          cat("[OK] fgsea Hallmark: ", nrow(fgsea_hallmark), " pathways tested, ",
              n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")
          
          fgsea_sanity_tracker <<- rbind(
            fgsea_sanity_tracker,
            data.frame(
              Contrast = contrast_name,
              Database = "Hallmark",
              N_Input_Genes = length(ranked_genes),
              N_Pathways_Tested = nrow(fgsea_hallmark),
              N_Sig_Pathways = n_sig,
              N_Sig_Up = n_sig_up,
              N_Sig_Down = n_sig_down,
              stringsAsFactors = FALSE
            )
          )
        }
      }, error = function(e) {
        cat("[WARN] fgsea Hallmark failed: ", conditionMessage(e), "\n", sep = "")
        fgsea_sanity_tracker <<- rbind(
          fgsea_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Database = "Hallmark",
            N_Input_Genes = length(ranked_genes),
            N_Pathways_Tested = 0,
            N_Sig_Pathways = 0,
            N_Sig_Up = 0,
            N_Sig_Down = 0,
            stringsAsFactors = FALSE
          )
        )
      })
      log_time(paste0("Completed Hallmark fgsea for ", contrast_name))
    } else if (run_hallmark) {
      cat("[WARN] Hallmark gene sets not available\n")
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
        if (!"Description" %in% colnames(sig_results)) {
          sig_results$Description <- sig_results$pathway
        }
        
        ## Convert list columns to character (leadingEdge is a list)
        for (col_name in colnames(sig_results)) {
          if (is.list(sig_results[[col_name]])) {
            sig_results[[col_name]] <- sapply(sig_results[[col_name]], function(x) {
              if (is.null(x) || length(x) == 0) {
                return(NA_character_)
              } else {
                return(paste(x, collapse = ";"))
              }
            })
          }
        }
        
        ## Save table
        table_file <- paste0("fgsea_", db_name, "_", contrast_name, "_", run_tag, ".csv")
        write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
        cat("[OK] Saved fgsea table: ", table_file, "\n", sep = "")
        
        if (db_name == "gobp_simplified") {
          next
        }
        
        ## Create combined plot (balanced up/down)
        top_up <- sig_results %>%
          dplyr::filter(NES > 0) %>%
          dplyr::arrange(padj, dplyr::desc(NES)) %>%
          dplyr::slice_head(n = top_n_per_direction)
        
        top_down <- sig_results %>%
          dplyr::filter(NES < 0) %>%
          dplyr::arrange(padj, NES) %>%
          dplyr::slice_head(n = top_n_per_direction)
        
        plot_data <- dplyr::bind_rows(top_down, top_up) %>%
          dplyr::arrange(NES) %>%
          dplyr::mutate(
            pathway_label = ifelse(!is.na(Description) & Description != "", Description, pathway),
            pathway_label = factor(pathway_label, levels = pathway_label),
            Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")
          )
        
        if (nrow(plot_data) > 0) {
          ## Use consistent colors: orange for up, blue for down
          p_fgsea <- ggplot(plot_data, aes(x = NES, y = pathway_label, fill = Direction)) +
            geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
            scale_fill_manual(
              values = c("Up-regulated" = "#E69F00", "Down-regulated" = "#0072B2"),
              name = "Direction"
            ) +
            labs(
              title = paste0("fgsea ", toupper(db_name), ": ", contrast_name),
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
            ) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5)
          
          plot_file <- paste0("fgsea_plot_", db_name, "_", contrast_name, "_", run_tag, ".png")
          ggsave(file.path(outdir, "plots", plot_file), plot = p_fgsea, width = 12,
                 height = max(6, nrow(plot_data) * 0.35), dpi = 300, bg = "white")
          cat("[OK] Saved fgsea plot: ", plot_file, "\n", sep = "")
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
    log_time(paste0("Starting fgsea for contrast: ", contrast_name))
    
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
    n_pos <- sum(ranked_vec > 0, na.rm = TRUE)
    n_neg <- sum(ranked_vec < 0, na.rm = TRUE)
    n_zero <- sum(ranked_vec == 0, na.rm = TRUE)
    n_total <- length(ranked_vec)
    pct_pos <- if (n_total == 0) 0 else round(100 * n_pos / n_total, 1)
    pct_neg <- if (n_total == 0) 0 else round(100 * n_neg / n_total, 1)
    pct_zero <- if (n_total == 0) 0 else round(100 * n_zero / n_total, 1)
    cat("[INFO] Ranked stats distribution: +", n_pos, " (", pct_pos, "%), -",
        n_neg, " (", pct_neg, "%), 0=", n_zero, " (", pct_zero, "%)\n", sep = "")
    if (pct_neg < 5) {
      cat("[WARN] Low fraction of negative stats; down-regulated pathways may be scarce.\n")
    }
    
    cat("[INFO] Ranked gene list for fgsea: ", length(ranked_vec), " genes\n", sep = "")
    
    ## Run fgsea
    if (length(ranked_vec) >= 10) {
      run_fgsea_analysis(
        ranked_genes  = ranked_vec,
        contrast_name = contrast_name,
        go_bp_list    = go_bp_list,
        go_bp_term2gene = go_bp_term2gene,
        kegg_list     = kegg_list,
        wikipathways_list = wikipathways_list,
        hallmark_list = hallmark_list,
        run_tag       = run_tag,
        outdir        = outdir,
        fdr_cutoff    = fdr_cut,
        simplify_go = simplify_go_bp,
        simplify_cutoff = simplify_go_cutoff,
        top_n_per_direction = top_n_per_direction
      )
      log_time(paste0("Completed fgsea for contrast: ", contrast_name))
    } else {
      cat("[WARN] Too few genes for fgsea\n")
    }
  }
  
  ## 4.5) Save sanity check table ----------------------------
  sanity_file <- paste0("fgsea_sanity_check_", run_tag, ".csv")
  write.csv(fgsea_sanity_tracker, file = file.path(outdir, "tables", sanity_file), row.names = FALSE)
  cat("\n[OK] Saved fgsea sanity check table: ", sanity_file, "\n", sep = "")
  
  ## 4.6) Print comprehensive sanity check summary -----------
  cat("\n==========================================================\n")
  cat("=== SANITY CHECK SUMMARY: fgsea PATHWAY ANALYSIS ===\n")
  cat("==========================================================\n")
  cat("Note: fgsea uses ranked gene list as universe (no separate background needed)\n")
  cat("Ranking metric: DESeq2 Wald statistic\n\n")
  
  ## Summarize across all contrasts
  for (cn_name in unique(fgsea_sanity_tracker$Contrast)) {
    cat("--- ", cn_name, " ---\n", sep = "")
    subset_data <- fgsea_sanity_tracker[fgsea_sanity_tracker$Contrast == cn_name, ]
    
    for (i in seq_len(nrow(subset_data))) {
      row <- subset_data[i, ]
      cat("  ", row$Database, ":\n", sep = "")
      cat("    Ranked gene list size: ", row$N_Input_Genes, "\n", sep = "")
      cat("    Pathways tested: ", row$N_Pathways_Tested, "\n", sep = "")
      cat("    Significant pathways (padj<0.05): ", row$N_Sig_Pathways, "\n", sep = "")
      if (row$N_Sig_Pathways > 0) {
        cat("      Up-regulated (NES>0): ", row$N_Sig_Up, "\n", sep = "")
        cat("      Down-regulated (NES<0): ", row$N_Sig_Down, "\n", sep = "")
      }
    }
    cat("\n")
  }
  cat("==========================================================\n\n")
  
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
