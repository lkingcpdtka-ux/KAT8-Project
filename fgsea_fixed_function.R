## =============================================================================
## CORRECTED fgsea Function (fixes GO:BP and KEGG errors)
## =============================================================================

run_fgsea_analysis <- function(ranked_genes, contrast_name, run_tag, outdir, fdr_cutoff = 0.05) {
  
  cat("\n--- fgsea (GSEA): ", contrast_name, " ---\n", sep = "")
  cat("[INFO] Ranked gene list size: ", length(ranked_genes), " genes\n", sep = "")
  
  if (length(ranked_genes) < 10) {
    cat("[WARN] Too few genes for fgsea\n")
    return(NULL)
  }
  
  results <- list()
  
  ## ===================
  ## GO:BP gene sets (FIXED)
  ## ===================
  tryCatch({
    cat("[INFO] Running fgsea for GO:BP...\n")
    
    ## Get GO:BP gene sets from org.Mm.eg.db
    go_bp_genes <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys = keys(org.Mm.eg.db, keytype = "GOALL"),
      columns = c("SYMBOL", "GOALL", "ONTOLOGYALL"),
      keytype = "GOALL"
    )
    
    ## Filter for BP ontology and remove NA symbols
    go_bp_genes <- go_bp_genes %>%
      dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
      dplyr::select(GOALL, SYMBOL)
    
    ## Convert to named list (pathway name -> gene vector)
    go_bp_list <- split(go_bp_genes$SYMBOL, go_bp_genes$GOALL)
    
    ## Filter for minimum gene set size
    go_bp_list <- go_bp_list[sapply(go_bp_list, length) >= 5 & sapply(go_bp_list, length) <= 500]
    
    cat("[INFO] Testing ", length(go_bp_list), " GO:BP gene sets\n", sep = "")
    
    ## Run fgsea
    fgsea_go <- fgsea(
      pathways = go_bp_list,
      stats    = ranked_genes,
      minSize  = 5,
      maxSize  = 500,
      nPermSimple = 10000
    )
    
    if (!is.null(fgsea_go) && nrow(fgsea_go) > 0) {
      ## FIXED: Get GO term descriptions from GO.db instead of org.Mm.eg.db
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
  
  ## ===================
  ## KEGG gene sets (FIXED - alternative approach)
  ## ===================
  tryCatch({
    cat("[INFO] Running fgsea for KEGG...\n")
    
    ## FIXED: Use enrichKEGG to get pathway gene sets, then convert to fgsea format
    ## Get all mouse genes with Entrez IDs
    all_entrez <- keys(org.Mm.eg.db, keytype = "ENTREZID")
    
    ## Get KEGG pathways using enrichKEGG (this works reliably)
    kegg_pathways <- enrichKEGG(
      gene = all_entrez[1:min(100, length(all_entrez))],  ## dummy gene list just to get pathway structure
      organism = "mmu",
      pvalueCutoff = 1,  ## get all pathways
      qvalueCutoff = 1
    )
    
    ## Extract pathway to gene mappings from KEGG database
    ## Alternative: use msigdbr for KEGG if available
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      cat("[INFO] Using msigdbr for KEGG pathways\n")
      kegg_msigdb <- msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
      kegg_list <- split(kegg_msigdb$gene_symbol, kegg_msigdb$gs_name)
    } else {
      ## Fallback: build from org.Mm.eg.db
      cat("[INFO] Building KEGG pathways from org.Mm.eg.db\n")
      kegg_genes <- AnnotationDbi::select(
        org.Mm.eg.db,
        keys = keys(org.Mm.eg.db, keytype = "PATH"),
        columns = c("SYMBOL", "PATH"),
        keytype = "PATH"
      )
      kegg_genes <- kegg_genes %>% dplyr::filter(!is.na(SYMBOL), !is.na(PATH))
      kegg_list <- split(kegg_genes$SYMBOL, paste0("mmu", kegg_genes$PATH))
    }
    
    ## Filter for gene set size
    kegg_list <- kegg_list[sapply(kegg_list, length) >= 5 & sapply(kegg_list, length) <= 500]
    
    cat("[INFO] Testing ", length(kegg_list), " KEGG pathways\n", sep = "")
    
    if (length(kegg_list) > 0) {
      ## Run fgsea
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
  
  ## ===================
  ## Save results with DIRECTION-SPECIFIC COLORS (matching pathway bar plots)
  ## ===================
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
      
      ## Create plot - separate by direction (NES > 0 = Up, NES < 0 = Down)
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
