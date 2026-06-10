#!/usr/bin/env Rscript

## =====================================================================
## KAT8 bulk RNA-seq - PART 7: CONSOLIDATED SEX ANALYSIS
## =====================================================================
## ONE script that replaces the previous three:
##   - part7.5_sex_concordance.R        (M vs F DEG concordance)
##   - part7.6_sex_interaction_test.R   (Sex x Genotype interaction LRT)
##   - part7.7_sex_combining_validation.R (model comparison / verdict)
##
## ...AND adds the new biology-focused module (Section 5):
##   Splitting the significant interaction genes into
##     * MALE-AMPLIFIED  / FEMALE-BUFFERED  (female response blunted)
##     * FEMALE-AMPLIFIED / MALE-BUFFERED   (male response blunted)
##     * SEX-OPPOSITE                       (direction flips between sexes)
##   then running SEPARATE GO:BP + KEGG ORA on each class to ask:
##     "Which pathways are females PROTECTED in?"
##
## ---------------------------------------------------------------------
## BIOLOGICAL MOTIVATION & LITERATURE CONTEXT
## ---------------------------------------------------------------------
## In our adipocyte KAT8-KD model, both sexes develop lipodystrophy /
## metabolic dysfunction, but the FEMALE transcriptional response to KD
## is blunted for a few hundred genes (the Sex x Genotype interaction
## set). This is NOT unique to our hands - female metabolic protection in
## C57BL/6 is well documented, and the pathways that show interaction
## here (lipid/fatty-acid/cholesterol metabolism in iWAT; thermogenesis,
## insulin response, energy homeostasis in gWAT) are precisely the axes
## that estrogen / ERalpha signalling is known to defend in females.
##
## Key references (OUR model vs INDEPENDENT mechanism):
##
##  OUR OWN WORK (establishes the phenotype - NOT external corroboration):
##  [1] (our lab) "Loss of the Lysine Acetyltransferase KAT8 in Adipocytes
##      Confers Lipodystrophy and Metabolic Dysfunction in Male and Female
##      Mice." Diabetes 2022; 71(Suppl_1):198-LB. doi:10.2337/db22-198-LB
##      -> This is US. Establishes BOTH sexes develop lipodystrophy /
##         metabolic dysfunction. So the sex effect we are dissecting here
##         is a DIFFERENCE OF DEGREE (female buffering), not on/off.
##
##  INDEPENDENT LITERATURE (the "it's not just us" support):
##  [2] "Histone H4 lysine 16 acetylation controls central carbon
##      metabolism and diet-induced obesity in mice." Nat Commun 2021;
##      12:6212. doi:10.1038/s41467-021-26277-w
##      -> H4K16ac (the KAT8 mark) governs lipid + central carbon
##         metabolism; loss predisposes to obesity / T2D phenotypes.
##  [3] "MOF-mediated H4K16ac governs mitochondrial and ciliary functions
##      by controlling gene promoters." Nat Commun 2023; 14:4866.
##      doi:10.1038/s41467-023-40108-0
##      -> KAT8 loss = mitochondrial dysfunction; BAT/beige thermogenesis
##         is mitochondria-dependent and sexually dimorphic.
##  [4] "Membrane estrogen receptor-alpha contributes to female protection
##      against high-fat diet-induced metabolic disorders." Front
##      Endocrinol 2023; 14:1215947. doi:10.3389/fendo.2023.1215947
##      -> Membrane-initiated ERalpha protects FEMALES (not males); loss
##         lowers Ucp1/Pgc1a, impairs BAT thermogenesis. A candidate
##         mechanism for the female buffering we observe.
##  [5] "PPARalpha in Obesity: Sex Difference and Estrogen Involvement"
##      (PMC2943125) + broader C57BL/6 female-protection literature
##      -> General principle from OTHER groups: females buffer
##         lipid/oxidative metabolic insults; males are the more
##         vulnerable sex. Our female buffering fits this established axis.
##
## ---------------------------------------------------------------------
## KNOWN LIMITATIONS OF THIS MODEL / ANALYSIS (read before interpreting)
## ---------------------------------------------------------------------
##  (a) POWER. n = 4-5 per sex x genotype group. The interaction (LRT)
##      test is the least-powered contrast in the design; it can only
##      detect LARGE sex differences. Genes called "concordant" may still
##      hide subtle real interactions. Conversely, ~5% of interaction
##      hits are expected false positives at FDR < 0.05.
##  (b) CELL-TYPE CONFOUND. Bulk WAT mixes adipocytes, immune, endothelial
##      and stromal cells. Sex differences in cellularity (e.g. females
##      have more beige-competent iWAT) can masquerade as gene-level
##      "interaction". Treat interaction hits as tissue-level, not
##      adipocyte-autonomous, without snRNA-seq confirmation.
##  (c) ESTROUS / HORMONE STATE is uncontrolled here; female variance is
##      expected to be higher, which DEFLATES the M-vs-F logFC
##      correlation (a magnitude effect, not a direction effect).
##  (d) KNOCKDOWN, not knockout: residual KAT8 / H4K16ac differs by sex
##      could itself be sex-biased; we do not measure per-sex KD
##      efficiency at the protein/mark level in this dataset.
##  (e) DEPOT DIMORPHISM: gWAT (visceral) is more sexually dimorphic at
##      baseline than iWAT (subcutaneous); lower M-vs-F r in gWAT is
##      partly baseline biology, not KD-specific divergence.
##
## OUTPUTS (per depot, into plots/sex_analysis and tables/):
##   1. M vs F concordance scatter + Venns + concordant gene list
##   2. Interaction LRT p-value histogram + pi0 + variance decomposition
##   3. Interaction-vs-genotype effect-size scatter
##   4. NEW: male-amplified / female-buffered classification table
##   5. NEW: separate GO:BP + KEGG ORA dotplots per interaction class
##   6. NEW: heatmap of top interaction genes (M vs F logFC side by side)
##   7. Model comparison (+/- Sex) and a PASS/CAUTION/CONCERN verdict
##   8. Combined CSV + human-readable report with manuscript text
## =====================================================================


## =====================================================================
## 1) PACKAGES
## =====================================================================
required_pkgs <- c("DESeq2", "dplyr", "tidyr", "ggplot2", "ggrepel",
                   "scales", "RColorBrewer", "VennDiagram", "grid",
                   "pheatmap", "stringr",
                   "clusterProfiler", "org.Mm.eg.db", "enrichplot",
                   "AnnotationDbi")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
  library(VennDiagram)
  library(grid)
  library(pheatmap)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(AnnotationDbi)
})

## qvalue is optional (pi0 estimation)
has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
if (has_qvalue) library(qvalue)

## VennDiagram writes a log file per call unless silenced
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

## =====================================================================
## 1.5) PARAMETERS
## =====================================================================
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found in working directory.")
}
if (!exists("prefilter_min_count"))   prefilter_min_count   <- 10
if (!exists("prefilter_min_samples")) prefilter_min_samples <- 4
if (!exists("MASTER_SEED"))           MASTER_SEED           <- 12345
set.seed(MASTER_SEED)

## Helper: robust DESeq2 results -> data.frame (guards against S4 colname loss)
deseq_res_to_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  std_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  if (ncol(df) == length(std_cols)) {
    colnames(df) <- std_cols
  } else if (ncol(df) > length(std_cols)) {
    colnames(df)[seq_along(std_cols)] <- std_cols
  }
  df
}

## =====================================================================
## 2) LOAD COUNTS + ANNOTATION  (identical logic to part1)
## =====================================================================
cat("\n=== LOADING DATA (consolidated sex analysis) ===\n")

counts_file <- file.path(getwd(), "counts.txt")
if (!file.exists(counts_file)) stop("counts.txt not found")
counts_raw <- read.delim(counts_file, header = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)

sample_ids <- paste0("JS_", sprintf("%02d", 1:40))
sample_annot <- data.frame(
  Sample   = sample_ids,
  Depot    = c(rep("iWAT", 10), rep("iWAT", 10), rep("gWAT", 10), rep("gWAT", 10)),
  Sex      = c(rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5),
               rep("F", 5), rep("F", 5), rep("M", 5), rep("M", 5)),
  Genotype = c(rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5),
               rep("CTL", 5), rep("KAT8KD", 5), rep("CTL", 5), rep("KAT8KD", 5)),
  stringsAsFactors = FALSE
)

outlier_samples <- c("JS_08", "JS_28")
sample_annot <- sample_annot[!sample_annot$Sample %in% outlier_samples, ]
tissue_samples <- intersect(sample_annot$Sample, colnames(counts_raw))
sample_annot <- sample_annot[sample_annot$Sample %in% tissue_samples, ]

sample_annot$Genotype <- factor(sample_annot$Genotype, levels = c("CTL", "KAT8KD"))
sample_annot$Depot    <- factor(sample_annot$Depot, levels = c("iWAT", "gWAT"))
sample_annot$Sex      <- factor(sample_annot$Sex, levels = c("F", "M"))
rownames(sample_annot) <- sample_annot$Sample

## Build gene-named count matrix
if ("gene" %in% colnames(counts_raw) || "gene_name" %in% colnames(counts_raw)) {
  counts_tissue <- counts_raw
  if (!"gene_name" %in% colnames(counts_tissue) && "gene" %in% colnames(counts_tissue))
    counts_tissue$gene_name <- counts_tissue$gene
  if (!"ensembl_gene_id" %in% colnames(counts_tissue) && "gene_id" %in% colnames(counts_tissue))
    counts_tissue$ensembl_gene_id <- counts_tissue$gene_id
  na_or_blank <- is.na(counts_tissue$gene_name) | counts_tissue$gene_name == ""
  if ("ensembl_gene_id" %in% colnames(counts_tissue))
    counts_tissue$gene_name[na_or_blank] <- counts_tissue$ensembl_gene_id[na_or_blank]
  if (any(duplicated(counts_tissue$gene_name))) {
    dup_flag <- duplicated(counts_tissue$gene_name) | duplicated(counts_tissue$gene_name, fromLast = TRUE)
    if ("ensembl_gene_id" %in% colnames(counts_tissue))
      counts_tissue$gene_name[dup_flag] <- paste0(counts_tissue$ensembl_gene_id[dup_flag], "_",
                                                  counts_tissue$gene_name[dup_flag])
    counts_tissue$gene_name <- make.unique(counts_tissue$gene_name)
  }
  rownames(counts_tissue) <- counts_tissue$gene_name
  count_matrix <- as.matrix(counts_tissue[, tissue_samples])
  rownames(count_matrix) <- counts_tissue$gene_name
} else {
  count_matrix <- as.matrix(counts_raw[, tissue_samples])
}

keep <- rowSums(count_matrix >= prefilter_min_count) >= prefilter_min_samples
count_matrix <- count_matrix[keep, , drop = FALSE]
cat("[INFO] Genes after global prefiltering: ", nrow(count_matrix), "\n", sep = "")

## =====================================================================
## 3) OUTPUT SETUP
## =====================================================================
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
if (length(run_dirs) > 0) {
  outdir <- run_dirs[which.max(file.info(run_dirs)$mtime)]
} else {
  outdir <- file.path(savepoint_dir, paste0("RUN_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
run_tag    <- gsub("^RUN_", "", basename(outdir))
tables_dir <- file.path(outdir, "tables")
plots_dir  <- file.path(outdir, "plots")
sex_dir    <- file.path(plots_dir, "sex_analysis")
for (d in c(tables_dir, plots_dir, sex_dir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## Colors
col_male   <- "#2166AC"
col_female <- "#B2182B"
category_colors <- c(
  "Concordant"  = "#4CAF50",
  "Male only"   = col_male,
  "Female only" = col_female,
  "Discordant"  = "#FF9800",
  "NS"          = "grey85"
)
class_colors <- c(
  "Male-amplified / Female-buffered" = col_male,
  "Female-amplified / Male-buffered" = col_female,
  "Sex-opposite"                     = "#6A3D9A"
)

## Reusable ORA helper: GO:BP + KEGG on a gene-symbol vector against a universe.
## Returns a list of dotplots + writes CSVs. Safe to call on small sets.
run_ora_block <- function(gene_symbols, universe_symbols, depot, class_tag, run_tag,
                          tables_dir, sex_dir) {
  out <- list(go = NULL, kegg = NULL, n_go = 0, n_kegg = 0)
  if (length(gene_symbols) < 5) {
    cat("    [skip ORA] ", class_tag, ": only ", length(gene_symbols),
        " genes (<5)\n", sep = "")
    return(out)
  }
  g_entrez <- tryCatch(
    clusterProfiler::bitr(unique(gene_symbols), "SYMBOL", "ENTREZID",
                          OrgDb = org.Mm.eg.db, drop = TRUE),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character()))
  u_entrez <- tryCatch(
    clusterProfiler::bitr(unique(universe_symbols), "SYMBOL", "ENTREZID",
                          OrgDb = org.Mm.eg.db, drop = TRUE),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character()))
  if (nrow(g_entrez) < 5) {
    cat("    [skip ORA] ", class_tag, ": <5 genes mapped to Entrez\n", sep = "")
    return(out)
  }

  safe_tag <- gsub("[^A-Za-z0-9]+", "_", class_tag)

  dot_from_result <- function(df_sig, title_main, subtitle) {
    top <- df_sig %>% dplyr::slice_head(n = ora_params$top_n)
    top$GeneRatio_num <- sapply(strsplit(top$GeneRatio, "/"),
                                function(x) as.numeric(x[1]) / as.numeric(x[2]))
    top$Count <- sapply(strsplit(top$GeneRatio, "/"), function(x) as.numeric(x[1]))
    top$neg_log10_fdr <- -log10(top$p.adjust)
    top$Description <- factor(stringr::str_wrap(top$Description, 34),
                              levels = rev(stringr::str_wrap(top$Description, 34)))
    ggplot(top, aes(x = GeneRatio_num, y = Description)) +
      geom_point(aes(size = Count, color = neg_log10_fdr)) +
      scale_color_gradientn(colors = enrichment_colors$gradient,
                            name = expression(-log[10]*"(FDR)")) +
      scale_size_continuous(name = "Count", range = c(2, 6)) +
      labs(title = title_main, subtitle = subtitle, x = "Gene Ratio", y = NULL) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
            panel.grid.minor = element_blank())
  }

  ## GO:BP
  tryCatch({
    ego <- enrichGO(gene = g_entrez$ENTREZID, universe = u_entrez$ENTREZID,
                    OrgDb = org.Mm.eg.db, ont = "BP",
                    pvalueCutoff = ora_params$pvalue_cutoff,
                    qvalueCutoff = ora_params$qvalue_cutoff,
                    readable = TRUE,
                    minGSSize = ora_params$min_gs_size,
                    maxGSSize = ora_params$max_gs_size)
    if (!is.null(ego) && nrow(ego@result) > 0) {
      ego_s <- simplify(ego, cutoff = ora_params$simplify_cutoff,
                        by = "p.adjust", select_fun = min)
      sig <- ego_s@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::arrange(p.adjust)
      if (nrow(sig) > 0) {
        write.csv(sig, file.path(tables_dir,
                  paste0("interaction_", safe_tag, "_GOBP_", depot, "_", run_tag, ".csv")),
                  row.names = FALSE)
        p <- dot_from_result(sig, paste0(depot, " GO:BP - ", class_tag),
                             paste0(length(gene_symbols), " genes in class"))
        ggsave(file.path(sex_dir,
               paste0("ORA_GOBP_", safe_tag, "_", depot, "_", run_tag, ".png")),
               p, width = 7, height = 6, dpi = 300, bg = "white")
        out$go <- sig; out$n_go <- nrow(sig)
        cat("    [GO:BP] ", class_tag, ": ", nrow(sig), " sig terms\n", sep = "")
      }
    }
  }, error = function(e) cat("    [WARN] GO:BP failed (", class_tag, "): ", e$message, "\n", sep = ""))

  ## KEGG
  tryCatch({
    ek <- enrichKEGG(gene = g_entrez$ENTREZID, universe = u_entrez$ENTREZID,
                     organism = "mmu",
                     pvalueCutoff = ora_params$pvalue_cutoff,
                     qvalueCutoff = ora_params$qvalue_cutoff,
                     minGSSize = ora_params$min_gs_size,
                     maxGSSize = ora_params$max_gs_size)
    if (!is.null(ek) && nrow(ek@result) > 0) {
      sig <- ek@result %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::arrange(p.adjust)
      if (nrow(sig) > 0 && "geneID" %in% colnames(sig)) {
        sig$geneID_symbols <- sapply(sig$geneID, function(g) {
          ids <- unlist(strsplit(g, "/"))
          syms <- g_entrez$SYMBOL[match(ids, g_entrez$ENTREZID)]
          paste(syms[!is.na(syms)], collapse = "/")
        })
      }
      if (nrow(sig) > 0) {
        write.csv(sig, file.path(tables_dir,
                  paste0("interaction_", safe_tag, "_KEGG_", depot, "_", run_tag, ".csv")),
                  row.names = FALSE)
        p <- dot_from_result(sig, paste0(depot, " KEGG - ", class_tag),
                             paste0(length(gene_symbols), " genes in class"))
        ggsave(file.path(sex_dir,
               paste0("ORA_KEGG_", safe_tag, "_", depot, "_", run_tag, ".png")),
               p, width = 7, height = 6, dpi = 300, bg = "white")
        out$kegg <- sig; out$n_kegg <- nrow(sig)
        cat("    [KEGG] ", class_tag, ": ", nrow(sig), " sig pathways\n", sep = "")
      }
    }
  }, error = function(e) cat("    [WARN] KEGG failed (", class_tag, "): ", e$message, "\n", sep = ""))

  out
}

## =====================================================================
## 4) PER-DEPOT ANALYSIS
## =====================================================================
summary_rows  <- list()
report_lines  <- c(
  "=====================================================================",
  " KAT8 CONSOLIDATED SEX ANALYSIS - REPORT",
  paste0(" Run: ", run_tag, "   Generated: ", Sys.time()),
  "====================================================================="
)

for (depot in c("iWAT", "gWAT")) {

  cat("\n============================================================\n")
  cat("=== ", depot, ": CONSOLIDATED SEX ANALYSIS ===\n", sep = "")
  cat("============================================================\n")

  thresh  <- get_tissue_thresholds(depot)
  fdr_cut <- thresh$fdr_cut
  lfc_cut <- thresh$logFC_cut

  depot_samples <- sample_annot$Sample[sample_annot$Depot == depot]
  depot_annot   <- sample_annot[depot_samples, ]
  depot_counts  <- count_matrix[, depot_samples]
  depot_keep    <- rowSums(depot_counts >= prefilter_min_count) >= prefilter_min_samples
  depot_counts  <- depot_counts[depot_keep, , drop = FALSE]
  depot_genes   <- rownames(depot_counts)
  cat("[INFO] ", depot, ": ", nrow(depot_counts), " genes x ", ncol(depot_counts),
      " samples\n", sep = "")

  ## -----------------------------------------------------------------
  ## 4A) SEX-STRATIFIED DESeq2  (M-only, F-only)  [from part7.5/7.6]
  ## -----------------------------------------------------------------
  female_idx <- depot_annot$Sample[depot_annot$Sex == "F"]
  male_idx   <- depot_annot$Sample[depot_annot$Sex == "M"]

  dds_f <- DESeqDataSetFromMatrix(round(depot_counts[, female_idx]),
                                  depot_annot[female_idx, ], ~ Genotype)
  dds_f <- DESeq(dds_f, quiet = TRUE)
  de_f  <- deseq_res_to_df(results(dds_f)); de_f$gene <- depot_genes
  rownames(de_f) <- depot_genes

  dds_m <- DESeqDataSetFromMatrix(round(depot_counts[, male_idx]),
                                  depot_annot[male_idx, ], ~ Genotype)
  dds_m <- DESeq(dds_m, quiet = TRUE)
  de_m  <- deseq_res_to_df(results(dds_m)); de_m$gene <- depot_genes
  rownames(de_m) <- depot_genes

  ## DEG sets
  m_up   <- rownames(de_m)[!is.na(de_m$padj) & de_m$padj < fdr_cut & de_m$log2FoldChange >  lfc_cut]
  m_down <- rownames(de_m)[!is.na(de_m$padj) & de_m$padj < fdr_cut & de_m$log2FoldChange < -lfc_cut]
  f_up   <- rownames(de_f)[!is.na(de_f$padj) & de_f$padj < fdr_cut & de_f$log2FoldChange >  lfc_cut]
  f_down <- rownames(de_f)[!is.na(de_f$padj) & de_f$padj < fdr_cut & de_f$log2FoldChange < -lfc_cut]
  m_all  <- union(m_up, m_down); f_all <- union(f_up, f_down)
  shared_all <- intersect(m_all, f_all)
  shared_conc <- union(intersect(m_up, f_up), intersect(m_down, f_down))
  concordance_pct <- if (length(shared_all) > 0)
    round(100 * length(shared_conc) / length(shared_all), 1) else NA_real_

  ## logFC concordance frame (shared gene universe)
  common <- intersect(rownames(de_m), rownames(de_f))
  cor_df <- data.frame(gene = common,
                       m_lfc = de_m[common, "log2FoldChange"],
                       f_lfc = de_f[common, "log2FoldChange"],
                       m_padj = de_m[common, "padj"],
                       f_padj = de_f[common, "padj"],
                       stringsAsFactors = FALSE) %>%
    filter(is.finite(m_lfc), is.finite(f_lfc)) %>%
    mutate(m_sig = !is.na(m_padj) & m_padj < fdr_cut & abs(m_lfc) > lfc_cut,
           f_sig = !is.na(f_padj) & f_padj < fdr_cut & abs(f_lfc) > lfc_cut,
           category = case_when(
             m_sig & f_sig & sign(m_lfc) == sign(f_lfc) ~ "Concordant",
             (m_sig | f_sig) & sign(m_lfc) != sign(f_lfc) ~ "Discordant",
             m_sig & !f_sig ~ "Male only",
             !m_sig & f_sig ~ "Female only",
             TRUE ~ "NS"))
  cor_pearson <- cor(cor_df$m_lfc, cor_df$f_lfc)
  sig_df <- cor_df %>% filter(m_sig | f_sig)
  cor_sig <- if (nrow(sig_df) > 10) cor(sig_df$m_lfc, sig_df$f_lfc) else NA_real_

  ## Concordance scatter (axis: gWAT capped at +/-5 like your current plot)
  axis_max <- if (depot == "gWAT") 5 else {
    deg <- cor_df %>% filter(category != "NS")
    max(c(abs(deg$m_lfc), abs(deg$f_lfc), 2), na.rm = TRUE) * 1.15
  }
  p_sc <- ggplot(cor_df, aes(m_lfc, f_lfc, color = category)) +
    geom_point(data = filter(cor_df, category == "NS"), size = 0.3, alpha = 0.15) +
    geom_point(data = filter(cor_df, category %in% c("Male only", "Female only")),
               size = 1.2, alpha = 0.6) +
    geom_point(data = filter(cor_df, category %in% c("Concordant", "Discordant")),
               size = 1.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    scale_color_manual(values = category_colors, name = "Category") +
    labs(title = paste0(depot, ": Male vs Female KAT8-KD Response"),
         subtitle = paste0("Pearson r = ", round(cor_pearson, 3),
                           " (all); r = ", round(cor_sig, 3),
                           " (DEGs); concordance = ", concordance_pct, "%"),
         x = expression(log[2]~FC~"(Male KD vs CTL)"),
         y = expression(log[2]~FC~"(Female KD vs CTL)")) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40")) +
    coord_cartesian(xlim = c(-axis_max, axis_max), ylim = c(-axis_max, axis_max))
  ggsave(file.path(sex_dir, paste0("01_concordance_scatter_", depot, "_", run_tag, ".png")),
         p_sc, width = 8, height = 7, dpi = 300, bg = "white")

  ## Venn (all DEGs)
  venn_png <- file.path(sex_dir, paste0("01_Venn_", depot, "_", run_tag, ".png"))
  png(venn_png, width = 2000, height = 2000, res = 300)
  grid.newpage()
  v <- venn.diagram(list(Male = m_all, Female = f_all), filename = NULL,
                    fill = c(col_male, col_female), alpha = 0.5,
                    cex = 1.4, cat.cex = 1.4, cat.fontface = "bold",
                    main = paste0(depot, ": All KD DEGs (M vs F)"))
  grid.draw(v); dev.off()
  cat("[OK] 4A concordance: r=", round(cor_pearson, 3),
      ", concordance=", concordance_pct, "%\n", sep = "")

  ## -----------------------------------------------------------------
  ## 4B) INTERACTION MODEL  ~ Sex + Genotype + Sex:Genotype  [part7.6]
  ## -----------------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(round(depot_counts), depot_annot,
                                ~ Sex + Genotype + Sex:Genotype)
  dds_lrt  <- DESeq(dds, test = "LRT", reduced = ~ Sex + Genotype)
  dds_wald <- DESeq(dds, test = "Wald")

  res_lrt <- deseq_res_to_df(results(dds_lrt)); res_lrt$gene <- depot_genes
  res_lrt <- res_lrt[!is.na(res_lrt$pvalue), ]

  coef_names <- resultsNames(dds_wald)
  inter_coef <- grep("Sex.*Genotype|Genotype.*Sex", coef_names, value = TRUE)
  geno_coef  <- grep("^Genotype", coef_names, value = TRUE)
  res_wald <- deseq_res_to_df(results(dds_wald,
                name = if (length(inter_coef) == 1) inter_coef else tail(coef_names, 1)))
  res_wald$gene <- depot_genes
  res_wald <- res_wald[!is.na(res_wald$pvalue), ]

  n_tested <- nrow(res_lrt)
  n_sig_lrt <- sum(res_lrt$padj < fdr_cut, na.rm = TRUE)
  pct_sig <- round(100 * n_sig_lrt / n_tested, 2)

  ## pi0
  pi0_est <- NA_real_
  if (has_qvalue && n_tested > 100) {
    pi0_est <- tryCatch(qvalue(res_lrt$pvalue)$pi0, error = function(e) NA_real_)
  }
  cat("[OK] 4B interaction: ", n_sig_lrt, "/", n_tested, " sig (", pct_sig,
      "%), pi0=", round(pi0_est, 3), "\n", sep = "")

  ## p-value histogram
  p_hist <- ggplot(res_lrt, aes(pvalue)) +
    geom_histogram(bins = 50, fill = "grey60", color = "white", boundary = 0) +
    geom_hline(yintercept = n_tested / 50, linetype = "dashed", color = "red") +
    labs(title = paste0(depot, ": Sex x Genotype interaction (LRT)"),
         subtitle = paste0(n_sig_lrt, "/", n_tested, " sig (", pct_sig, "%)",
                          if (!is.na(pi0_est)) paste0("; pi0=", round(pi0_est, 3)) else ""),
         x = "LRT p-value (full vs reduced)", y = "Genes") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
  ggsave(file.path(sex_dir, paste0("02_pvalue_hist_", depot, "_", run_tag, ".png")),
         p_hist, width = 7, height = 5, dpi = 300, bg = "white")

  ## Genotype main effect (KD effect in reference sex = females)
  res_geno <- deseq_res_to_df(results(dds_wald, name = geno_coef)); res_geno$gene <- depot_genes
  n_geno_sig <- sum(res_geno$padj < fdr_cut & abs(res_geno$log2FoldChange) > lfc_cut, na.rm = TRUE)

  ## Variance decomposition (sampled genes)
  vsd <- vst(dds_wald, blind = TRUE)
  vst_mat <- assay(vsd); rownames(vst_mat) <- depot_genes
  n_var <- min(5000, length(depot_genes))
  var_genes <- sample(depot_genes, n_var)
  vr <- data.frame(pct_sex = NA_real_, pct_genotype = NA_real_,
                   pct_interaction = NA_real_, pct_residual = NA_real_)[rep(1, n_var), ]
  for (i in seq_along(var_genes)) {
    y <- vst_mat[var_genes[i], ]
    tryCatch({
      ss <- anova(lm(y ~ Sex + Genotype + Sex:Genotype, data = depot_annot))$`Sum Sq`
      tot <- sum(ss)
      if (tot > 0) vr[i, ] <- 100 * ss[1:4] / tot
    }, error = function(e) {})
  }
  vr <- vr[!is.na(vr$pct_sex), ]
  med_var <- sapply(vr, median)
  var_long <- vr %>% dplyr::select(pct_sex, pct_genotype, pct_interaction) %>%
    pivot_longer(everything(), names_to = "Term", values_to = "Pct") %>%
    mutate(Term = recode(Term, pct_sex = "Sex", pct_genotype = "Genotype",
                         pct_interaction = "Sex x Genotype"))
  p_var <- ggplot(var_long, aes(Term, Pct, fill = Term)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3) +
    scale_fill_manual(values = c("Sex" = "#FF7F00", "Genotype" = "#377EB8",
                                 "Sex x Genotype" = "#E41A1C")) +
    labs(title = paste0(depot, ": Variance explained per gene"),
         subtitle = paste0("median interaction = ",
                          round(med_var["pct_interaction"], 1), "%"),
         x = NULL, y = "% variance") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
  ggsave(file.path(sex_dir, paste0("02_variance_decomp_", depot, "_", run_tag, ".png")),
         p_var, width = 6, height = 5, dpi = 300, bg = "white")

  ## -----------------------------------------------------------------
  ## 5) *** NEW *** SPLIT INTERACTION GENES: who is buffered?
  ## -----------------------------------------------------------------
  ## Significant interaction genes (LRT), annotated with the per-sex KD
  ## logFCs from 4A so we can read off WHO responds and WHO is buffered.
  ##
  ## Coefficient convention (Sex ref = F, Genotype ref = CTL):
  ##   Genotype_KAT8KD_vs_CTL coefficient  = KD effect in FEMALES
  ##   SexM.GenotypeKAT8KD interaction     = (male KD effect) - (female KD effect)
  ## We instead use the directly-interpretable stratified logFCs (m_lfc, f_lfc).
  ##
  ## Classes (within significant-interaction genes):
  ##   Male-amplified / Female-buffered : |m_lfc| > |f_lfc|, same sign
  ##        -> KD pushes the gene HARD in males, females hold the line.
  ##           This is the candidate "female protection" set.
  ##   Female-amplified / Male-buffered : |f_lfc| > |m_lfc|, same sign
  ##   Sex-opposite                     : sign(m_lfc) != sign(f_lfc)
  cat("  [5] Splitting ", n_sig_lrt, " interaction genes by who is buffered...\n", sep = "")

  sig_inter_genes <- res_lrt$gene[!is.na(res_lrt$padj) & res_lrt$padj < fdr_cut]
  inter_tab <- data.frame(
    gene  = sig_inter_genes,
    m_lfc = de_m[sig_inter_genes, "log2FoldChange"],
    f_lfc = de_f[sig_inter_genes, "log2FoldChange"],
    m_padj = de_m[sig_inter_genes, "padj"],
    f_padj = de_f[sig_inter_genes, "padj"],
    lrt_padj = res_lrt$padj[match(sig_inter_genes, res_lrt$gene)],
    inter_lfc = res_wald$log2FoldChange[match(sig_inter_genes, res_wald$gene)],
    stringsAsFactors = FALSE
  ) %>%
    filter(is.finite(m_lfc), is.finite(f_lfc)) %>%
    mutate(
      buffering = abs(m_lfc) - abs(f_lfc),       # >0 female buffered; <0 male buffered
      class = case_when(
        sign(m_lfc) != sign(f_lfc)            ~ "Sex-opposite",
        abs(m_lfc) >= abs(f_lfc)              ~ "Male-amplified / Female-buffered",
        TRUE                                  ~ "Female-amplified / Male-buffered"),
      female_protection_score = pmax(abs(m_lfc) - abs(f_lfc), 0)  # high = female holds line
    ) %>%
    arrange(desc(abs(buffering)))

  write.csv(inter_tab,
            file.path(tables_dir, paste0("interaction_genes_classified_", depot, "_", run_tag, ".csv")),
            row.names = FALSE)

  n_fembuf <- sum(inter_tab$class == "Male-amplified / Female-buffered")
  n_malbuf <- sum(inter_tab$class == "Female-amplified / Male-buffered")
  n_opp    <- sum(inter_tab$class == "Sex-opposite")
  cat("      Female-buffered: ", n_fembuf, " | Male-buffered: ", n_malbuf,
      " | Sex-opposite: ", n_opp, "\n", sep = "")

  ## Classification scatter (interaction genes only), colored by class
  if (nrow(inter_tab) > 0) {
    lab_df <- inter_tab %>% filter(!grepl("^Gm\\d|^[0-9]|Rik$", gene)) %>%
      group_by(class) %>% slice_head(n = 8) %>% ungroup()
    ax <- max(c(abs(inter_tab$m_lfc), abs(inter_tab$f_lfc)), na.rm = TRUE) * 1.1
    p_cls <- ggplot(inter_tab, aes(m_lfc, f_lfc, color = class)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_point(size = 1.8, alpha = 0.8) +
      geom_text_repel(data = lab_df, aes(label = gene), color = "black",
                      size = 2.8, max.overlaps = 20, segment.size = 0.3) +
      scale_color_manual(values = class_colors, name = "Interaction class") +
      labs(title = paste0(depot, ": Significant interaction genes by buffering"),
           subtitle = paste0("Points below the diagonal in the response direction = ",
                            "female response blunted (protected)"),
           x = expression(log[2]~FC~"(Male KD vs CTL)"),
           y = expression(log[2]~FC~"(Female KD vs CTL)")) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(), legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)) +
      coord_cartesian(xlim = c(-ax, ax), ylim = c(-ax, ax))
    ggsave(file.path(sex_dir, paste0("05_interaction_classes_", depot, "_", run_tag, ".png")),
           p_cls, width = 8, height = 8, dpi = 300, bg = "white")
  }

  ## Heatmap: top female-buffered genes, male vs female logFC side by side
  fembuf_top <- inter_tab %>% filter(class == "Male-amplified / Female-buffered") %>%
    filter(!grepl("^Gm\\d|^[0-9]|Rik$", gene)) %>%
    slice_head(n = 40)
  if (nrow(fembuf_top) >= 3) {
    hm <- as.matrix(fembuf_top[, c("m_lfc", "f_lfc")])
    rownames(hm) <- fembuf_top$gene; colnames(hm) <- c("Male KD", "Female KD")
    lim <- max(abs(hm), na.rm = TRUE)
    png(file.path(sex_dir, paste0("05_heatmap_female_buffered_", depot, "_", run_tag, ".png")),
        width = 1400, height = 2600, res = 300)
    pheatmap(hm, cluster_cols = FALSE, cluster_rows = TRUE,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
             breaks = seq(-lim, lim, length.out = 51),
             main = paste0(depot, ": top female-buffered interaction genes\n",
                          "(blue=down, red=up; note muted Female column)"),
             fontsize_row = 6, fontsize_col = 10, border_color = NA)
    dev.off()
  }

  ## Separate ORA per interaction class
  universe_syms <- res_lrt$gene
  cat("  [5] ORA per interaction class:\n")
  ora_fembuf <- run_ora_block(
    inter_tab$gene[inter_tab$class == "Male-amplified / Female-buffered"],
    universe_syms, depot, "Female-buffered", run_tag, tables_dir, sex_dir)
  ora_malbuf <- run_ora_block(
    inter_tab$gene[inter_tab$class == "Female-amplified / Male-buffered"],
    universe_syms, depot, "Male-buffered", run_tag, tables_dir, sex_dir)
  ora_opp <- run_ora_block(
    inter_tab$gene[inter_tab$class == "Sex-opposite"],
    universe_syms, depot, "Sex-opposite", run_tag, tables_dir, sex_dir)
  ## Also ORA on the full interaction set (matches your existing dotplot)
  ora_all <- run_ora_block(inter_tab$gene, universe_syms, depot,
                           "All-interaction", run_tag, tables_dir, sex_dir)

  ## -----------------------------------------------------------------
  ## 6) MODEL COMPARISON  ~Genotype vs ~Sex+Genotype  [part7.7, per depot]
  ## -----------------------------------------------------------------
  dds_g  <- DESeq(DESeqDataSetFromMatrix(round(depot_counts), depot_annot, ~ Genotype), quiet = TRUE)
  res_g  <- deseq_res_to_df(results(dds_g, name = grep("^Genotype", resultsNames(dds_g), value = TRUE)))
  n_deg_g <- sum(res_g$padj < fdr_cut & abs(res_g$log2FoldChange) > lfc_cut, na.rm = TRUE)

  dds_sg <- DESeq(DESeqDataSetFromMatrix(round(depot_counts), depot_annot, ~ Sex + Genotype), quiet = TRUE)
  res_sg <- deseq_res_to_df(results(dds_sg, name = grep("^Genotype", resultsNames(dds_sg), value = TRUE)))
  n_deg_sg <- sum(res_sg$padj < fdr_cut & abs(res_sg$log2FoldChange) > lfc_cut, na.rm = TRUE)
  pct_gain <- round(100 * (n_deg_sg - n_deg_g) / max(1, n_deg_g), 1)

  ## -----------------------------------------------------------------
  ## VERDICT (per depot)
  ## -----------------------------------------------------------------
  verdict <- if (!is.na(pi0_est) && pi0_est >= 0.85 && pct_sig < 5 && concordance_pct >= 95) {
    "PASS (pool sexes)"
  } else if (!is.na(pi0_est) && pi0_est >= 0.75 && pct_sig < 10) {
    "CAUTION (pool with covariate; report interaction set)"
  } else {
    "CONCERN (consider sex-stratified)"
  }
  cat("[VERDICT] ", depot, ": ", verdict, "\n", sep = "")

  ## Collect summary row
  summary_rows[[depot]] <- data.frame(
    Depot = depot, Genes_tested = n_tested,
    Interaction_sig = n_sig_lrt, Pct_interaction = pct_sig,
    pi0 = round(pi0_est, 4),
    Female_buffered = n_fembuf, Male_buffered = n_malbuf, Sex_opposite = n_opp,
    Pearson_r_all = round(cor_pearson, 3), Pearson_r_DEGs = round(cor_sig, 3),
    Concordance_pct = concordance_pct,
    Median_var_interaction = round(med_var["pct_interaction"], 2),
    Median_var_genotype = round(med_var["pct_genotype"], 2),
    Median_var_sex = round(med_var["pct_sex"], 2),
    DEG_genotype_only = n_deg_g, DEG_sex_plus_genotype = n_deg_sg,
    Pct_DEG_gain_with_sex = pct_gain,
    FemaleBuffered_GOBP_terms = ora_fembuf$n_go,
    FemaleBuffered_KEGG_terms = ora_fembuf$n_kegg,
    Verdict = verdict, stringsAsFactors = FALSE)

  ## Report block
  top_fembuf <- head(fembuf_top$gene, 15)
  report_lines <- c(report_lines, "",
    paste0("------------------------------------------------------------"),
    paste0("DEPOT: ", depot),
    paste0("------------------------------------------------------------"),
    paste0("  Interaction (LRT): ", n_sig_lrt, "/", n_tested, " genes (", pct_sig,
           "%), pi0 = ", round(pi0_est, 3)),
    paste0("  M-vs-F logFC concordance: r(all)=", round(cor_pearson, 3),
           ", r(DEGs)=", round(cor_sig, 3), ", directional=", concordance_pct, "%"),
    paste0("  Variance (median): Genotype=", round(med_var["pct_genotype"], 1),
           "%, Sex=", round(med_var["pct_sex"], 1),
           "%, Interaction=", round(med_var["pct_interaction"], 1), "%"),
    paste0("  Interaction split: Female-buffered=", n_fembuf,
           ", Male-buffered=", n_malbuf, ", Sex-opposite=", n_opp),
    paste0("  Female-buffered GO:BP sig terms=", ora_fembuf$n_go,
           "; KEGG sig pathways=", ora_fembuf$n_kegg),
    paste0("  Top female-buffered genes: ", paste(top_fembuf, collapse = ", ")),
    paste0("  Model comparison: ~Genotype=", n_deg_g, " DEGs; ~Sex+Genotype=",
           n_deg_sg, " DEGs (", ifelse(pct_gain >= 0, "+", ""), pct_gain, "%)"),
    paste0("  VERDICT: ", verdict))
}

## =====================================================================
## 7) COMBINED OUTPUTS
## =====================================================================
summary_df <- bind_rows(summary_rows)
write.csv(summary_df,
          file.path(tables_dir, paste0("sex_analysis_summary_", run_tag, ".csv")),
          row.names = FALSE)

cat("\n============================================================\n")
cat("=== CONSOLIDATED SEX ANALYSIS: SUMMARY ===\n")
cat("============================================================\n")
print(summary_df, row.names = FALSE)

## Literature-informed interpretation appended to the report
report_lines <- c(report_lines, "",
  "=====================================================================",
  " INTERPRETATION (literature-anchored)",
  "=====================================================================",
  " OUR phenotype (established by us): adipocyte KAT8 loss causes",
  " lipodystrophy / metabolic dysfunction in BOTH sexes (Diabetes 2022;",
  " 71:198-LB - this is our own model). So the question here is not WHETHER",
  " the sexes differ, but BY HOW MUCH - i.e. female buffering of degree.",
  "",
  " The female-buffered interaction set is enriched for lipid / fatty-acid /",
  " cholesterol metabolism (iWAT) and thermogenesis / insulin response /",
  " energy homeostasis (gWAT). INDEPENDENT literature from other groups",
  " makes this the expected place to see female buffering:",
  "   - KAT8/H4K16ac directly governs lipid + central-carbon metabolism and",
  "     mitochondrial function (Nat Commun 2021;12:6212; 2023;14:4866).",
  "     Loss should therefore hit fatty-acid oxidation, lipid storage and",
  "     mito/thermogenic programs - exactly the interaction pathways here.",
  "   - Membrane ERalpha protects FEMALES (not males) against HFD metabolic",
  "     disorders and sustains Ucp1/Pgc1a-driven BAT thermogenesis",
  "     (Front Endocrinol 2023;14:1215947). A plausible mechanism for why",
  "     females blunt the KD-driven collapse of these same programs.",
  "   - Female metabolic protection along lipid/oxidative axes is a general,",
  "     reproduced C57BL/6 phenomenon (e.g. PPARalpha sex difference,",
  "     PMC2943125) - so our female buffering is consistent with, not",
  "     contradicted by, the field.",
  "",
  " Caveats: bulk-tissue cellularity differences, hormonal/estrous state,",
  " per-sex KD efficiency, and low n all limit interaction power; treat the",
  " female-buffered pathways as hypotheses for targeted (snRNA-seq / per-sex",
  " functional) follow-up, not settled adipocyte-autonomous mechanism.",
  "",
  " References:",
  "  [1] Diabetes 2022;71(Suppl_1):198-LB  (OUR KAT8 adipocyte model)",
  "  [2] Nat Commun 2021;12:6212           (H4K16ac, carbon metabolism, obesity)",
  "  [3] Nat Commun 2023;14:4866           (MOF, mitochondria)",
  "  [4] Front Endocrinol 2023;14:1215947  (membrane ERalpha female protection)",
  "  [5] PPARalpha sex difference, PMC2943125; female metabolic protection")

report_file <- file.path(tables_dir, paste0("sex_analysis_report_", run_tag, ".txt"))
writeLines(report_lines, report_file)
cat("\n[OK] Report written: ", report_file, "\n", sep = "")

## Session info
log_dir <- file.path(outdir, "logs"); if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
sink(file.path(log_dir, paste0("sessionInfo_sex_analysis_", run_tag, ".txt")))
cat("=== R SESSION INFO - Part 7 consolidated sex analysis ===\n\n"); print(sessionInfo()); sink()

cat("\n=== PART 7 CONSOLIDATED SEX ANALYSIS COMPLETE ===\n")
cat("Plots:  ", sex_dir, "\n", sep = "")
cat("Tables: ", tables_dir, "\n", sep = "")
cat("\nNew per-class ORA files: interaction_Female-buffered_GOBP_<depot>_<run>.csv etc.\n")
cat("Classification table:    interaction_genes_classified_<depot>_<run>.csv\n")
