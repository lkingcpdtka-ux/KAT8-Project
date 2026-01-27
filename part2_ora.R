#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 2: ORA PATHWAY ANALYSIS
## =========================================================
## This script performs Over-Representation Analysis (ORA)
## for GO:BP and KEGG pathways using DEG lists from Part 1
##
## FIXES APPLIED:
## - Proper ORA universe/background from all tested genes
## - Better error handling
## - Direction-specific pathway bar plots
## =========================================================

## 0) Working dir (optional) --------------------------------
## setwd("C:/Users/lking/OneDrive - Louisiana State University/PBRC/Bioinformatics/KAT8KD_RNAseq")

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot",
  "DOSE",
  "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(DOSE)
  library(AnnotationDbi)
})

## Note: Improved mapping not needed - current mapping rate is 90.8% (excellent!)

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
## ============================================================
## DEPOT-SPECIFIC CUTOFFS
## ============================================================
## Set different thresholds for iWAT vs gWAT if desired
## This allows relaxed cutoffs for gWAT while keeping iWAT stringent
## ============================================================

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
## IMPORTANT: ORA Universe/Background Setting
## ============================================================
## Set this to TRUE or FALSE to control background gene set
##
## use_universe = TRUE (RECOMMENDED - statistically correct)
##   - Uses all genes tested by DESeq2 as background
##   - More conservative (fewer significant pathways)
##   - Accounts for your specific gene filtering
##
## use_universe = FALSE (common default)
##   - No background specified (uses all genes in database)
##   - More liberal (more significant pathways)
##   - Standard approach in many publications
##
## CURRENT SETTING:
use_universe <- TRUE   ## Using background gene set (statistically rigorous)
## ============================================================

## Plot/output parameters (used in captions)
ora_enrich_params <- list(
  pvalue_cutoff  = 0.1,
  qvalue_cutoff  = 0.2,
  min_gs_size    = 5,
  max_gs_size    = 500,
  simplify_cutoff = 0.7
)

ora_plot_params <- list(
  top_n    = 15,
  order_by = "adj. p-value"
)

## ============================================================
## COMPLETE ORA PARAMETER DOCUMENTATION
## ============================================================
## This section documents all ORA parameters for reference
##
## METHOD: Over-Representation Analysis (Fisher's exact test)
##   - Tests if DEG list has more genes in a pathway than expected by chance
##   - Uses hypergeometric distribution
##   - Requires discrete gene lists (up/down analyzed separately)
##
## INPUT GENE SELECTION:
##   - logFC_cut: 1 (genes with |log2FC| > 1, i.e., 2-fold change)
##   - fdr_cut: 0.05 (Benjamini-Hochberg adjusted p-value)
##   - Minimum genes for analysis: 5 (pathways with <5 genes skipped)
##
## ENRICHMENT PARAMETERS (GO:BP):
##   - pvalueCutoff: 0.1 (initial p-value threshold for enrichment)
##   - qvalueCutoff: 0.2 (initial FDR threshold for enrichment)
##   - Final reporting cutoff: fdr_cut (0.05, more stringent)
##   - Gene set size: minGSSize=5, maxGSSize=500
##   - Simplification: Yes (cutoff=0.7, removes redundant GO terms)
##   - readable: TRUE (converts Entrez IDs back to gene symbols)
##
## ENRICHMENT PARAMETERS (KEGG):
##   - pvalueCutoff: 0.1
##   - qvalueCutoff: 0.2
##   - organism: "mmu" (Mus musculus)
##   - Gene set size: minGSSize=5, maxGSSize=500
##
## GENE ID CONVERSION:
##   - Input: Gene symbols (SYMBOL)
##   - Conversion: SYMBOL -> ENTREZID (via org.Mm.eg.db)
##   - Typical mapping rate: ~91% (excellent)
##   - Failed mappings are dropped
##
## BACKGROUND/UNIVERSE:
##   - If use_universe=TRUE: All genes tested by DESeq2
##   - If use_universe=FALSE: All genes in organism annotation
##   - Universe also converted to Entrez IDs
##
## DATABASES:
##   - GO:BP: Gene Ontology Biological Process (hierarchical, comprehensive)
##   - KEGG: Kyoto Encyclopedia of Genes and Genomes (curated pathways)
##
## OUTPUT:
##   - CSV tables: Significant pathways (FDR < 0.05)
##   - Bar plots: Top 15 pathways per direction/database
##   - Columns: pathway, description, GeneRatio, p.adjust, genes, etc.
##
## SANITY TRACKING:
##   - pathway_sanity_tracker logs: contrast, direction, database,
##     N_Input_Genes, N_Mapped_Entrez, N_Universe_Entrez, N_Sig_Pathways
## ============================================================

## Sanity check tracker
pathway_sanity_tracker <- data.frame(
  Contrast = character(),
  Direction = character(),
  Database = character(),
  N_Input_Genes = integer(),
  N_Mapped_Entrez = integer(),
  N_Universe_Entrez = integer(),
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
  
  ## 4.2) Optionally create ORA universe --------------------
  if (use_universe) {
    cat("\n=== CREATING ORA UNIVERSE ===\n")
    cat("[INFO] Using background gene set (statistically rigorous approach)\n")
    
    ## Load first DE table (any will work; they all have the same tested genes)
    de_first <- read.csv(de_files[1], row.names = 1, stringsAsFactors = FALSE)
    
    ## Universe = all genes tested by DESeq2
    universe_symbols <- rownames(de_first)
    cat("[INFO] Total genes tested by DESeq2: ", length(universe_symbols), "\n")

    ## Convert to Entrez IDs
    universe_entrez_df <- tryCatch({
      bitr(
        universe_symbols,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Mm.eg.db,
        drop     = TRUE
      )
    }, error = function(e) {
      cat("[ERROR] Universe gene ID conversion failed: ", conditionMessage(e), "\n")
      stop("Cannot create universe for ORA")
    })

    ## Get unique Entrez IDs
    universe_entrez <- unique(universe_entrez_df$ENTREZID)
    universe_entrez <- universe_entrez[!is.na(universe_entrez)]
    
    cat("[INFO] Universe Entrez IDs: ", length(universe_entrez), " (",
        round(100 * length(universe_entrez) / length(universe_symbols), 1), "% mapped)\n", sep = "")
  } else {
    cat("\n=== ORA UNIVERSE SETTING ===\n")
    cat("[INFO] NOT using background gene set (default/liberal approach)\n")
    cat("[INFO] To use background, set use_universe <- TRUE at top of script\n")
    universe_entrez <- NULL
  }
  
  ## 4.3) ORA analysis function ------------------------------
  run_ora_analysis <- function(gene_list, direction, contrast_name,
                               universe_entrez, run_tag, outdir, fdr_cutoff = 0.05,
                               logfc_cutoff = 1.0) {

    cat("\n--- ORA Pathway enrichment: ", direction, " genes (", contrast_name, ") ---\n", sep = "")
    cat("[INFO] Input gene count: ", length(gene_list), " genes\n", sep = "")

    if (length(gene_list) < 5) {
      cat("[WARN] Too few genes (", length(gene_list), ") for pathway analysis\n", sep = "")
      return(NULL)
    }

    ## Convert gene list to Entrez IDs
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

    mapping_rate <- if (length(gene_list) == 0) 0 else nrow(gene_entrez) / length(unique(gene_list))
    cat("[OK] Mapped ", nrow(gene_entrez), " of ", length(unique(gene_list)),
        " genes to Entrez (", round(mapping_rate * 100, 1), "%)\n", sep = "")
    cat("[INFO] Genes entering pathway enrichment: ", nrow(gene_entrez), "\n", sep = "")
    if (!is.null(universe_entrez)) {
      cat("[INFO] Universe size: ", length(universe_entrez), " Entrez IDs\n", sep = "")
    } else {
      cat("[INFO] Universe: using default (all genes in database)\n")
    }
    
    if (nrow(gene_entrez) < 3) {
      cat("[WARN] Too few genes mapped to Entrez IDs\n")
      ## Track in sanity checker
      pathway_sanity_tracker <<- rbind(
        pathway_sanity_tracker,
        data.frame(
          Contrast = contrast_name,
          Direction = direction,
          Database = "GO:BP",
          N_Input_Genes = length(gene_list),
          N_Mapped_Entrez = nrow(gene_entrez),
          N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
          N_Sig_Pathways = 0,
          stringsAsFactors = FALSE
        )
      )
      pathway_sanity_tracker <<- rbind(
        pathway_sanity_tracker,
        data.frame(
          Contrast = contrast_name,
          Direction = direction,
          Database = "KEGG",
          N_Input_Genes = length(gene_list),
          N_Mapped_Entrez = nrow(gene_entrez),
          N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
          N_Sig_Pathways = 0,
          stringsAsFactors = FALSE
        )
      )
      return(NULL)
    }
    
    entrez_ids <- gene_entrez$ENTREZID
    results <- list()
    
    ## GO:BP enrichment (conditionally with universe)
    tryCatch({
      enrich_go <- if (!is.null(universe_entrez)) {
        enrichGO(
          gene          = entrez_ids,
          OrgDb         = org.Mm.eg.db,
          ont           = "BP",
          pvalueCutoff  = ora_enrich_params$pvalue_cutoff,
          qvalueCutoff  = ora_enrich_params$qvalue_cutoff,
          readable      = TRUE,
          minGSSize     = ora_enrich_params$min_gs_size,
          maxGSSize     = ora_enrich_params$max_gs_size,
          universe      = universe_entrez
        )
      } else {
        enrichGO(
          gene          = entrez_ids,
          OrgDb         = org.Mm.eg.db,
          ont           = "BP",
          pvalueCutoff  = ora_enrich_params$pvalue_cutoff,
          qvalueCutoff  = ora_enrich_params$qvalue_cutoff,
          readable      = TRUE,
          minGSSize     = ora_enrich_params$min_gs_size,
          maxGSSize     = ora_enrich_params$max_gs_size
        )
      }
      if (!is.null(enrich_go) && nrow(enrich_go@result) > 0) {
        enrich_go_simp <- simplify(
          enrich_go,
          cutoff = ora_enrich_params$simplify_cutoff,
          by = "p.adjust",
          select_fun = min
        )
        results$gobp <- enrich_go_simp
        n_sig <- sum(enrich_go_simp@result$p.adjust < fdr_cutoff, na.rm = TRUE)
        cat("[OK] GO:BP enrichment: ", nrow(enrich_go_simp@result), " terms (simplified), ",
            n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")
        
        ## Track in sanity checker
        pathway_sanity_tracker <<- rbind(
          pathway_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Direction = direction,
            Database = "GO:BP",
            N_Input_Genes = length(gene_list),
            N_Mapped_Entrez = nrow(gene_entrez),
            N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
            N_Sig_Pathways = n_sig,
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("[INFO] No GO:BP terms found\n")
        pathway_sanity_tracker <<- rbind(
          pathway_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Direction = direction,
            Database = "GO:BP",
            N_Input_Genes = length(gene_list),
            N_Mapped_Entrez = nrow(gene_entrez),
            N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
            N_Sig_Pathways = 0,
            stringsAsFactors = FALSE
          )
        )
      }
    }, error = function(e) {
      cat("[WARN] GO:BP enrichment failed: ", conditionMessage(e), "\n", sep = "")
      pathway_sanity_tracker <<- rbind(
        pathway_sanity_tracker,
        data.frame(
          Contrast = contrast_name,
          Direction = direction,
          Database = "GO:BP",
          N_Input_Genes = length(gene_list),
          N_Mapped_Entrez = nrow(gene_entrez),
          N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
          N_Sig_Pathways = 0,
          stringsAsFactors = FALSE
        )
      )
    })
    
    ## KEGG enrichment (conditionally with universe)
    tryCatch({
      enrich_kegg <- if (!is.null(universe_entrez)) {
        enrichKEGG(
          gene         = entrez_ids,
          organism     = "mmu",
          pvalueCutoff = ora_enrich_params$pvalue_cutoff,
          qvalueCutoff = ora_enrich_params$qvalue_cutoff,
          minGSSize    = ora_enrich_params$min_gs_size,
          maxGSSize    = ora_enrich_params$max_gs_size,
          universe     = universe_entrez
        )
      } else {
        enrichKEGG(
          gene         = entrez_ids,
          organism     = "mmu",
          pvalueCutoff = ora_enrich_params$pvalue_cutoff,
          qvalueCutoff = ora_enrich_params$qvalue_cutoff,
          minGSSize    = ora_enrich_params$min_gs_size,
          maxGSSize    = ora_enrich_params$max_gs_size
        )
      }
      
      if (!is.null(enrich_kegg) && nrow(enrich_kegg@result) > 0) {
        ## Convert Entrez IDs back to gene symbols for KEGG results
        kegg_result <- enrich_kegg@result
        
        ## Parse geneID column (format: "entrez1/entrez2/entrez3")
        if ("geneID" %in% colnames(kegg_result)) {
          kegg_result$geneID_symbols <- sapply(kegg_result$geneID, function(gene_str) {
            entrez_vec <- unlist(strsplit(gene_str, "/"))
            symbol_vec <- gene_entrez$SYMBOL[match(entrez_vec, gene_entrez$ENTREZID)]
            symbol_vec <- symbol_vec[!is.na(symbol_vec)]
            paste(symbol_vec, collapse = "/")
          })
        }
        
        enrich_kegg@result <- kegg_result
        results$kegg <- enrich_kegg
        n_sig <- sum(enrich_kegg@result$p.adjust < fdr_cutoff, na.rm = TRUE)
        cat("[OK] KEGG enrichment: ", nrow(enrich_kegg@result), " pathways, ",
            n_sig, " significant (FDR<", fdr_cutoff, ")\n", sep = "")
        
        ## Track in sanity checker
        pathway_sanity_tracker <<- rbind(
          pathway_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Direction = direction,
            Database = "KEGG",
            N_Input_Genes = length(gene_list),
            N_Mapped_Entrez = nrow(gene_entrez),
            N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
            N_Sig_Pathways = n_sig,
            stringsAsFactors = FALSE
          )
        )
      } else {
        cat("[INFO] No KEGG pathways found\n")
        pathway_sanity_tracker <<- rbind(
          pathway_sanity_tracker,
          data.frame(
            Contrast = contrast_name,
            Direction = direction,
            Database = "KEGG",
            N_Input_Genes = length(gene_list),
            N_Mapped_Entrez = nrow(gene_entrez),
            N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
            N_Sig_Pathways = 0,
            stringsAsFactors = FALSE
          )
        )
      }
    }, error = function(e) {
      cat("[WARN] KEGG enrichment failed: ", conditionMessage(e), "\n", sep = "")
      pathway_sanity_tracker <<- rbind(
        pathway_sanity_tracker,
        data.frame(
          Contrast = contrast_name,
          Direction = direction,
          Database = "KEGG",
          N_Input_Genes = length(gene_list),
          N_Mapped_Entrez = nrow(gene_entrez),
          N_Universe_Entrez = ifelse(is.null(universe_entrez), NA, length(universe_entrez)),
          N_Sig_Pathways = 0,
          stringsAsFactors = FALSE
        )
      )
    })

    ## Save results tables + barplots
    if (length(results) > 0) {
      for (db_name in names(results)) {
        enrich_obj <- results[[db_name]]
        
        sig_results <- enrich_obj@result %>%
          dplyr::filter(p.adjust < fdr_cutoff) %>%
          dplyr::arrange(p.adjust)
        
        if (nrow(sig_results) == 0) {
          cat("[INFO] No significant pathways in ", db_name, "\n", sep = "")
          next
        }
        
        ## Remove duplicate columns before saving
        sig_results <- sig_results[, !duplicated(colnames(sig_results)), drop = FALSE]

        ## Record analysis thresholds for downstream captions
        sig_results$logFC_cutoff <- logfc_cutoff
        sig_results$fdr_cutoff <- fdr_cutoff
        
        ## Save table
        table_file <- paste0("ORA_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".csv")
        write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
        cat("[OK] Saved pathway table: ", table_file, " (", nrow(sig_results), " terms)\n", sep = "")
        
        ## Create bar plot (top 15 terms) - ordered by adjusted p-value, direction-specific colors
        plot_data <- sig_results %>%
          dplyr::slice_head(n = ora_plot_params$top_n)
        
        ## Parse GeneRatio to numeric
        plot_data$GeneRatio_numeric <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
          num <- as.numeric(x[1]); denom <- as.numeric(x[2])
          if (is.na(denom) || denom == 0) return(0)
          num / denom
        })
        
        ## Order by adjusted p-value (most significant at top)
        plot_data <- plot_data %>%
          dplyr::arrange(dplyr::desc(p.adjust)) %>%
          dplyr::mutate(Description = factor(Description, levels = Description))
        
        ## Direction-specific color gradient (orange for Up, blue for Down)
        if (direction == "Up") {
          fill_scale <- scale_fill_gradient(
            low   = "#FDD49E",  ## Light orange
            high  = "#E69F00",  ## Dark orange
            name  = "Adj.\nP-value"
          )
        } else {
          fill_scale <- scale_fill_gradient(
            low   = "#9ECAE1",  ## Light blue
            high  = "#0072B2",  ## Dark blue/teal
            name  = "Adj.\nP-value"
          )
        }
        
        ## Parameter caption
        universe_label <- if (is.null(universe_entrez)) {
          "Universe: all genes (default)"
        } else {
          paste0("Universe: ", length(universe_entrez), " DESeq2-tested genes")
        }
        param_caption <- paste(
          paste0(
            "FDR < ", fdr_cutoff,
            " | |log2FC| > ", logfc_cutoff,
            " | Top ", ora_plot_params$top_n, " pathways"
          ),
          paste0("Ordered by ", ora_plot_params$order_by, " | ", universe_label),
          sep = "\n"
        )

        p_pathway <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description, fill = p.adjust)) +
          geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
          fill_scale +
          labs(
            title = paste0(toupper(db_name), " Pathways (ORA - ", direction, ")\n", contrast_name),
            x     = "Gene Ratio",
            y     = NULL,
            caption = param_caption
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
            panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
            plot.caption    = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8))
          )
        
        plot_file <- paste0("ORA_barplot_", db_name, "_", contrast_name, "_", direction, "_", run_tag, ".png")
        ggsave(
          file.path(outdir, "plots", plot_file),
          plot   = p_pathway,
          width  = 10,
          height = max(6, nrow(plot_data) * 0.3),
          dpi    = 300,
          bg     = "white"
        )
        cat("[OK] Saved pathway bar plot: ", plot_file, "\n", sep = "")
      }
    }
    
    return(results)
  }
  
  ## 4.4) Process each DE table ------------------------------
  cat("\n=== RUNNING ORA FOR ALL CONTRASTS ===\n")
  
  for (de_file in de_files) {
    
    ## Extract contrast name from filename
    de_filename <- basename(de_file)
    contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)
    
    cat("\n=== Processing contrast: ", contrast_name, " ===\n", sep = "")

    ## Select depot-specific cutoffs
    if (grepl("iWAT", contrast_name, ignore.case = TRUE)) {
      logFC_cut <- iWAT_logFC_cut
      fdr_cut <- iWAT_fdr_cut
      cat("[INFO] Using iWAT cutoffs: FDR < ", fdr_cut, ", |logFC| > ", logFC_cut, "\n", sep = "")
    } else if (grepl("gWAT", contrast_name, ignore.case = TRUE)) {
      logFC_cut <- gWAT_logFC_cut
      fdr_cut <- gWAT_fdr_cut
      cat("[INFO] Using gWAT cutoffs: FDR < ", fdr_cut, ", |logFC| > ", logFC_cut, "\n", sep = "")
    } else {
      logFC_cut <- default_logFC_cut
      fdr_cut <- default_fdr_cut
      cat("[INFO] Using default cutoffs: FDR < ", fdr_cut, ", |logFC| > ", logFC_cut, "\n", sep = "")
    }

    ## Load DE table
    de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)

    ## Ensure column names are correct
    if (!"adj.P.Val" %in% colnames(de_table)) {
      if ("padj" %in% colnames(de_table)) {
        de_table$adj.P.Val <- de_table$padj
      } else {
        cat("[WARN] Cannot find FDR column in ", de_file, "\n")
        next
      }
    }

    if (!"logFC" %in% colnames(de_table)) {
      if ("log2FoldChange" %in% colnames(de_table)) {
        de_table$logFC <- de_table$log2FoldChange
      } else {
        cat("[WARN] Cannot find logFC column in ", de_file, "\n")
        next
      }
    }

    ## Get up-regulated genes (using depot-specific cutoffs)
    up_genes <- de_table %>%
      dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut, logFC > logFC_cut) %>%
      rownames()

    ## Get down-regulated genes (using depot-specific cutoffs)
    down_genes <- de_table %>%
      dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut, logFC < -logFC_cut) %>%
      rownames()

    cat("[INFO] Up-regulated DEGs: ", length(up_genes), "\n", sep = "")
    cat("[INFO] Down-regulated DEGs: ", length(down_genes), "\n", sep = "")
    
    ## Run ORA for up-regulated genes
    if (length(up_genes) >= 5) {
      run_ora_analysis(
        gene_list       = up_genes,
        direction       = "Up",
        contrast_name   = contrast_name,
        universe_entrez = universe_entrez,
        run_tag         = run_tag,
        outdir          = outdir,
        fdr_cutoff      = fdr_cut,
        logfc_cutoff    = logFC_cut
      )
    } else {
      cat("[WARN] Too few up-regulated genes for ORA\n")
    }

    ## Run ORA for down-regulated genes
    if (length(down_genes) >= 5) {
      run_ora_analysis(
        gene_list       = down_genes,
        direction       = "Down",
        contrast_name   = contrast_name,
        universe_entrez = universe_entrez,
        run_tag         = run_tag,
        outdir          = outdir,
        fdr_cutoff      = fdr_cut,
        logfc_cutoff    = logFC_cut
      )
    } else {
      cat("[WARN] Too few down-regulated genes for ORA\n")
    }
  }

  ## 4.4B) OPTIONAL: Re-run GWAT contrasts with RELAXED parameters ----
  ## Uncomment this section if GWAT has <50 down-regulated DEGs
  ##
  ## cat("\n\n==========================================================\n")
  ## cat("=== RELAXED ORA FOR GWAT (if needed) ===\n")
  ## cat("==========================================================\n\n")
  ##
  ## logFC_cut_relaxed <- 0.5   # More lenient (was 1.0)
  ## fdr_cut_relaxed   <- 0.1   # More lenient (was 0.05)
  ##
  ## cat("RELAXED PARAMETERS (GWAT only):\n")
  ## cat("  DEG thresholds: FDR < ", fdr_cut_relaxed, ", |logFC| > ", logFC_cut_relaxed, "\n", sep = "")
  ## cat("  Pathway p-cutoff: 0.15 (was 0.1)\n")
  ## cat("  Pathway q-cutoff: 0.25 (was 0.2)\n\n")
  ##
  ## gwat_de_files <- de_files[grepl("gWAT", de_files)]
  ##
  ## for (de_file in gwat_de_files) {
  ##   de_filename <- basename(de_file)
  ##   contrast_name <- gsub("^DE_tissue_(.+)_[0-9]{8}_[0-9]{6}\\.csv$", "\\1", de_filename)
  ##
  ##   cat("\n=== RELAXED: ", contrast_name, " ===\n", sep = "")
  ##   de_table <- read.csv(de_file, row.names = 1, stringsAsFactors = FALSE)
  ##
  ##   if (!"adj.P.Val" %in% colnames(de_table)) {
  ##     de_table$adj.P.Val <- de_table$padj
  ##   }
  ##   if (!"logFC" %in% colnames(de_table)) {
  ##     de_table$logFC <- de_table$log2FoldChange
  ##   }
  ##
  ##   ## Get down genes with RELAXED cutoffs
  ##   down_genes_relaxed <- de_table %>%
  ##     dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut_relaxed,
  ##                   logFC < -logFC_cut_relaxed) %>%
  ##     rownames()
  ##
  ##   cat("[INFO] Down DEGs (relaxed): ", length(down_genes_relaxed), "\n", sep = "")
  ##
  ##   if (length(down_genes_relaxed) >= 5) {
  ##     ## Convert to Entrez
  ##     gene_entrez <- tryCatch({
  ##       bitr(unique(down_genes_relaxed), fromType = "SYMBOL", toType = "ENTREZID",
  ##            OrgDb = org.Mm.eg.db, drop = TRUE)
  ##     }, error = function(e) {
  ##       data.frame(SYMBOL = character(), ENTREZID = character())
  ##     })
  ##
  ##     if (nrow(gene_entrez) >= 3) {
  ##       ## GO:BP with RELAXED parameters
  ##       enrich_go_relaxed <- if (!is.null(universe_entrez)) {
  ##         enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
  ##                  pvalueCutoff = 0.15, qvalueCutoff = 0.25, readable = TRUE,
  ##                  minGSSize = 3, maxGSSize = 500, universe = universe_entrez)
  ##       } else {
  ##         enrichGO(gene = gene_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
  ##                  pvalueCutoff = 0.15, qvalueCutoff = 0.25, readable = TRUE,
  ##                  minGSSize = 3, maxGSSize = 500)
  ##       }
  ##
  ##       if (!is.null(enrich_go_relaxed) && nrow(enrich_go_relaxed@result) > 0) {
  ##         enrich_go_simp <- simplify(enrich_go_relaxed, cutoff = 0.7)
  ##         sig_results <- enrich_go_simp@result %>%
  ##           dplyr::filter(p.adjust < 0.25) %>%
  ##           dplyr::arrange(p.adjust)
  ##
  ##         if (nrow(sig_results) > 0) {
  ##           table_file <- paste0("ORA_RELAXED_gobp_", contrast_name, "_Down_", run_tag, ".csv")
  ##           write.csv(sig_results, file = file.path(outdir, "tables", table_file), row.names = FALSE)
  ##           cat("[OK] RELAXED ORA: ", nrow(sig_results), " pathways (saved as ", table_file, ")\n", sep = "")
  ##         }
  ##       }
  ##     }
  ##   }
  ## }

  ## 4.5) Save sanity check table ----------------------------
  sanity_file <- paste0("ORA_sanity_check_", run_tag, ".csv")
  write.csv(pathway_sanity_tracker, file = file.path(outdir, "tables", sanity_file), row.names = FALSE)
  cat("\n[OK] Saved ORA sanity check table: ", sanity_file, "\n", sep = "")
  
  ## 4.6) Print comprehensive sanity check summary -----------
  cat("\n==========================================================\n")
  cat("=== SANITY CHECK SUMMARY: ORA PATHWAY ANALYSIS ===\n")
  cat("==========================================================\n")
  cat("Universe setting: use_universe = ", use_universe, "\n", sep = "")
  if (use_universe) {
    cat("  (Conservative: uses all DESeq2-tested genes as background)\n")
  } else {
    cat("  (Liberal: uses all genes in database as background)\n")
  }
  cat("\n")
  
  ## Summarize across all contrasts
  for (cn_name in unique(pathway_sanity_tracker$Contrast)) {
    cat("--- ", cn_name, " ---\n", sep = "")
    subset_data <- pathway_sanity_tracker[pathway_sanity_tracker$Contrast == cn_name, ]
    
    for (i in seq_len(nrow(subset_data))) {
      row <- subset_data[i, ]
      cat("  ", row$Direction, " | ", row$Database, ":\n", sep = "")
      cat("    Input genes: ", row$N_Input_Genes, "\n", sep = "")
      cat("    Mapped to Entrez: ", row$N_Mapped_Entrez,
          " (", round(100 * row$N_Mapped_Entrez / row$N_Input_Genes, 1), "%)\n", sep = "")
      if (!is.na(row$N_Universe_Entrez)) {
        cat("    Universe size: ", row$N_Universe_Entrez, "\n", sep = "")
      }
      cat("    Significant pathways (padj<0.05): ", row$N_Sig_Pathways, "\n", sep = "")
    }
    cat("\n")
  }
  cat("==========================================================\n\n")
  
  cat("\n======================================\n")
  cat("=== PART 2 (ORA) COMPLETE ===\n")
  cat("======================================\n\n")

  ## 4.7) Save session info for reproducibility -------------
  session_file <- file.path(outdir, "logs", paste0("sessionInfo_ORA_", run_tag, ".txt"))
  sink(session_file)
  cat("=== R SESSION INFORMATION ===\n\n")
  print(sessionInfo())
  cat("\n=== KEY PACKAGE VERSIONS ===\n\n")
  key_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "dplyr", "ggplot2")
  print(installed.packages()[intersect(key_pkgs, rownames(installed.packages())), c("Version", "Built")])
  sink()
  cat("[INFO] Session info saved to: ", session_file, "\n", sep = "")

}, error = function(e) {
  
  cat("\n======================================\n")
  cat("=== ERROR IN ORA PIPELINE ===\n")
  cat("======================================\n")
  cat(conditionMessage(e), "\n\n")
  
  err_file <- paste0("ERROR_ORA_", run_tag, ".txt")
  writeLines(
    c(
      "ORA script error:",
      conditionMessage(e),
      "",
      "Traceback:",
      capture.output(traceback())
    ),
    con = file.path(outdir, "logs", err_file)
  )
  
  stop(e)
})

cat("\n=== Part 2 (ORA) finished ===\n")
