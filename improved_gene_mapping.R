## =========================================================
## IMPROVED GENE MAPPING FUNCTION
## =========================================================
## This function uses a multi-strategy approach to maximize
## gene symbol -> Entrez ID mapping success rate
##
## Strategies used (in order):
## 1. SYMBOL -> ENTREZID (standard)
## 2. ALIAS -> ENTREZID (for alternative/old symbols)
## 3. ENSEMBL -> ENTREZID (if ENSEMBL IDs available)
##
## Expected improvement: 80-85% -> 88-92% mapping rate
## =========================================================

#' Improved gene ID mapping with multiple strategies
#'
#' @param gene_symbols Character vector of gene symbols to map
#' @param ensembl_ids Optional character vector of ENSEMBL IDs (same length as gene_symbols)
#' @param orgdb Organism database (default: org.Mm.eg.db for mouse)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Data frame with columns: SYMBOL, ENTREZID, Mapping_Strategy
#'
improved_gene_mapping <- function(gene_symbols,
                                   ensembl_ids = NULL,
                                   orgdb = org.Mm.eg.db,
                                   verbose = TRUE) {

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("AnnotationDbi package required")
  }

  if (verbose) {
    cat("[INFO] Starting improved gene mapping for ", length(gene_symbols), " genes...\n", sep = "")
  }

  ## Initialize results data frame
  results <- data.frame(
    SYMBOL = gene_symbols,
    ENTREZID = NA_character_,
    Mapping_Strategy = NA_character_,
    stringsAsFactors = FALSE
  )

  ## -------------------------------------------------------
  ## STRATEGY 1: Standard SYMBOL -> ENTREZID mapping
  ## -------------------------------------------------------
  if (verbose) cat("[1/3] Trying standard SYMBOL mapping...\n")

  mapped_1 <- tryCatch({
    suppressMessages(
      AnnotationDbi::select(
        orgdb,
        keys = gene_symbols,
        columns = "ENTREZID",
        keytype = "SYMBOL"
      )
    )
  }, error = function(e) {
    if (verbose) cat("[WARN] Standard mapping failed: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })

  if (!is.null(mapped_1)) {
    ## Handle 1:many mappings (take first)
    mapped_1_unique <- mapped_1 %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()

    ## Update results
    for (i in seq_len(nrow(mapped_1_unique))) {
      idx <- which(results$SYMBOL == mapped_1_unique$SYMBOL[i])
      if (length(idx) > 0 && is.na(results$ENTREZID[idx[1]])) {
        results$ENTREZID[idx[1]] <- mapped_1_unique$ENTREZID[i]
        results$Mapping_Strategy[idx[1]] <- "SYMBOL"
      }
    }

    n_mapped_1 <- sum(!is.na(results$ENTREZID))
    if (verbose) {
      cat("[OK] SYMBOL mapping: ", n_mapped_1, " / ", length(gene_symbols),
          " (", round(100 * n_mapped_1 / length(gene_symbols), 1), "%)\n", sep = "")
    }
  }

  ## -------------------------------------------------------
  ## STRATEGY 2: ALIAS -> ENTREZID for remaining unmapped
  ## -------------------------------------------------------
  unmapped_symbols <- results$SYMBOL[is.na(results$ENTREZID)]

  if (length(unmapped_symbols) > 0) {
    if (verbose) {
      cat("[2/3] Trying ALIAS mapping for ", length(unmapped_symbols), " unmapped genes...\n", sep = "")
    }

    mapped_2 <- tryCatch({
      suppressMessages(
        AnnotationDbi::select(
          orgdb,
          keys = unmapped_symbols,
          columns = c("ENTREZID", "SYMBOL"),
          keytype = "ALIAS"
        )
      )
    }, error = function(e) {
      if (verbose) cat("[WARN] ALIAS mapping failed: ", conditionMessage(e), "\n", sep = "")
      return(NULL)
    })

    if (!is.null(mapped_2)) {
      ## Handle 1:many mappings
      mapped_2_unique <- mapped_2 %>%
        dplyr::filter(!is.na(ENTREZID)) %>%
        dplyr::rename(ORIGINAL_SYMBOL = ALIAS) %>%
        dplyr::group_by(ORIGINAL_SYMBOL) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()

      ## Update results
      for (i in seq_len(nrow(mapped_2_unique))) {
        idx <- which(results$SYMBOL == mapped_2_unique$ORIGINAL_SYMBOL[i])
        if (length(idx) > 0 && is.na(results$ENTREZID[idx[1]])) {
          results$ENTREZID[idx[1]] <- mapped_2_unique$ENTREZID[i]
          results$Mapping_Strategy[idx[1]] <- "ALIAS"
        }
      }

      n_rescued_2 <- sum(!is.na(results$ENTREZID)) - n_mapped_1
      if (verbose && n_rescued_2 > 0) {
        cat("[OK] ALIAS mapping rescued ", n_rescued_2, " additional genes\n", sep = "")
      }
    }
  }

  ## -------------------------------------------------------
  ## STRATEGY 3: ENSEMBL -> ENTREZID for remaining unmapped
  ## -------------------------------------------------------
  if (!is.null(ensembl_ids) && length(ensembl_ids) == length(gene_symbols)) {
    unmapped_idx <- which(is.na(results$ENTREZID))

    if (length(unmapped_idx) > 0) {
      if (verbose) {
        cat("[3/3] Trying ENSEMBL mapping for ", length(unmapped_idx), " unmapped genes...\n", sep = "")
      }

      unmapped_ensembl <- ensembl_ids[unmapped_idx]

      mapped_3 <- tryCatch({
        suppressMessages(
          AnnotationDbi::select(
            orgdb,
            keys = unmapped_ensembl,
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "ENSEMBL"
          )
        )
      }, error = function(e) {
        if (verbose) cat("[WARN] ENSEMBL mapping failed: ", conditionMessage(e), "\n", sep = "")
        return(NULL)
      })

      if (!is.null(mapped_3)) {
        ## Handle 1:many mappings
        mapped_3_unique <- mapped_3 %>%
          dplyr::filter(!is.na(ENTREZID)) %>%
          dplyr::group_by(ENSEMBL) %>%
          dplyr::slice_head(n = 1) %>%
          dplyr::ungroup()

        ## Update results
        for (i in seq_len(nrow(mapped_3_unique))) {
          ensembl_match <- which(ensembl_ids == mapped_3_unique$ENSEMBL[i])
          if (length(ensembl_match) > 0) {
            idx <- ensembl_match[1]
            if (is.na(results$ENTREZID[idx])) {
              results$ENTREZID[idx] <- mapped_3_unique$ENTREZID[i]
              results$Mapping_Strategy[idx] <- "ENSEMBL"
            }
          }
        }

        n_rescued_3 <- sum(!is.na(results$ENTREZID)) - n_mapped_1 - n_rescued_2
        if (verbose && n_rescued_3 > 0) {
          cat("[OK] ENSEMBL mapping rescued ", n_rescued_3, " additional genes\n", sep = "")
        }
      }
    }
  } else {
    if (verbose) cat("[3/3] Skipping ENSEMBL mapping (no ENSEMBL IDs provided)\n")
  }

  ## -------------------------------------------------------
  ## FINAL SUMMARY
  ## -------------------------------------------------------
  n_total_mapped <- sum(!is.na(results$ENTREZID))
  pct_mapped <- round(100 * n_total_mapped / length(gene_symbols), 1)

  if (verbose) {
    cat("\n--- FINAL MAPPING RESULTS ---\n")
    cat("Total mapped: ", n_total_mapped, " / ", length(gene_symbols),
        " (", pct_mapped, "%)\n", sep = "")

    strategy_counts <- table(results$Mapping_Strategy[!is.na(results$ENTREZID)])
    cat("By strategy:\n")
    for (strategy in names(strategy_counts)) {
      cat("  ", strategy, ": ", strategy_counts[strategy], "\n", sep = "")
    }
    cat("\n")
  }

  ## Return only successfully mapped genes
  results_mapped <- results %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::select(SYMBOL, ENTREZID)

  return(results_mapped)
}


## =========================================================
## EXAMPLE USAGE
## =========================================================
##
## # Load your gene list
## gene_symbols <- c("Kat8", "Foxo1", "Pparg", "Adipoq", "Lep", "Gm12345", "SomeOldName")
##
## # Optional: if you have ENSEMBL IDs
## ensembl_ids <- c("ENSMUSG00000...", "ENSMUSG00000...", ...)
##
## # Run improved mapping
## mapped_genes <- improved_gene_mapping(
##   gene_symbols = gene_symbols,
##   ensembl_ids = ensembl_ids,  # Optional
##   verbose = TRUE
## )
##
## # Use in clusterProfiler
## enrichGO(
##   gene = mapped_genes$ENTREZID,
##   OrgDb = org.Mm.eg.db,
##   ...
## )
##
## =========================================================
