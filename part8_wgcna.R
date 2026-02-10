#!/usr/bin/env Rscript

## =========================================================
## KAT8 bulk RNA-seq - PART 8: WGCNA CO-EXPRESSION ANALYSIS
## =========================================================
## Weighted Gene Co-expression Network Analysis to identify
## co-regulated gene modules and their trait associations
##
## KEY ANALYSES:
## 1. Network Construction: Soft-thresholded adjacency matrix
## 2. Module Detection: Dynamic tree cutting
## 3. Module-Trait Correlation: Link modules to Tissue/Sex/Genotype
## 4. Hub Gene Identification: Top connectivity genes per module
## 5. Module Enrichment: GO/KEGG for key modules
##
## REQUIRES:
##   - Part 1 tissue results (savepoints/RUN_*/tables/DE_tissue_*.csv)
##   - Raw count data or VST matrix
##   - WGCNA package
## =========================================================

## 1) Packages ----------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "scales", "tidyr", "tibble", "RColorBrewer")
to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

required_bioc_pkgs <- c("DESeq2", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
to_install_bioc <- setdiff(required_bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(WGCNA)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

## Allow multi-threading for WGCNA
allowWGCNAThreads()

## 1.5) Load central parameters -----------------------------
params_file <- file.path(getwd(), "parameters.R")
if (file.exists(params_file)) {
  source(params_file)
} else {
  stop("parameters.R not found. Please ensure it exists in the project root.")
}

set.seed(MASTER_SEED)

## 2) save_core utilities -----------------------------------
utils_dir <- file.path(getwd(), "save_core")
if (file.exists(file.path(utils_dir, "save.R"))) {
  source(file.path(utils_dir, "save.R"))
  source(file.path(utils_dir, "findsave.R"))
  use_save_core <- TRUE
} else {
  cat("[WARN] save_core not found; proceeding without it\n")
  use_save_core <- FALSE
}

log_time <- function(message) {
  cat("[TIME] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", message, "\n", sep = "")
}

## ============================================================
## 3) LOAD DATA
## ============================================================

cat("\n=== LOADING COUNT DATA AND METADATA ===\n")

## Load raw counts
counts_raw <- read.delim("counts.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

## Build metadata from column names (tissue samples)
## Assumes tissue sample naming convention from part1_main_analysis.R
## Adjust this section to match your actual sample naming
all_samples <- colnames(counts_raw)
all_samples <- all_samples[!all_samples %in% c("gene_id", "gene")]

## Load metadata if available, otherwise construct from sample names
metadata_file <- file.path(getwd(), "data", "metadata.csv")
if (file.exists(metadata_file)) {
  sample_meta <- read.csv(metadata_file, stringsAsFactors = FALSE)
  cat("[OK] Loaded metadata from ", metadata_file, "\n", sep = "")
} else {
  cat("[WARN] No metadata.csv found. Attempting to construct from Part 1 outputs.\n")
  cat("[INFO] You may need to provide a metadata file with columns: Sample, Tissue, Sex, Genotype\n")

  ## Try to find the DESeq2 colData from an RDS file
  savepoint_dir <- file.path(getwd(), "savepoints")
  rds_files <- list.files(savepoint_dir, pattern = "dds.*\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(rds_files) > 0) {
    dds_obj <- readRDS(rds_files[length(rds_files)])
    sample_meta <- as.data.frame(colData(dds_obj))
    cat("[OK] Extracted metadata from saved DESeq2 object\n")
  } else {
    stop("Cannot find metadata. Please provide data/metadata.csv or ensure DESeq2 RDS exists.")
  }
}

## Ensure required columns
required_cols <- c("Tissue", "Genotype")
if (!all(required_cols %in% colnames(sample_meta))) {
  stop("Metadata must contain columns: ", paste(required_cols, collapse = ", "))
}

## Identify sample column
if ("Sample" %in% colnames(sample_meta)) {
  rownames(sample_meta) <- sample_meta$Sample
} else if (!is.null(rownames(sample_meta))) {
  sample_meta$Sample <- rownames(sample_meta)
}

## Match samples between counts and metadata
tissue_samples <- intersect(sample_meta$Sample, colnames(counts_raw))
if (length(tissue_samples) < 6) {
  stop("Too few matched samples (", length(tissue_samples), "). Check sample names.")
}

cat("[OK] Matched ", length(tissue_samples), " samples between counts and metadata\n", sep = "")

## Build count matrix
gene_col <- if ("gene" %in% colnames(counts_raw)) "gene" else "gene_id"
count_matrix <- counts_raw[, tissue_samples]
rownames(count_matrix) <- make.unique(counts_raw[[gene_col]])
count_matrix <- as.matrix(count_matrix)

sample_meta <- sample_meta[tissue_samples, , drop = FALSE]

## Set up output
savepoint_dir <- file.path(getwd(), "savepoints")
run_dirs <- list.dirs(savepoint_dir, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^RUN_", basename(run_dirs))]
run_dirs_sorted <- run_dirs[order(file.info(run_dirs)$mtime, decreasing = TRUE)]
outdir <- run_dirs_sorted[1]
run_tag <- gsub("^RUN_", "", basename(outdir))

tables_dir <- file.path(outdir, "tables")
plots_dir <- file.path(outdir, "plots")
wgcna_dir <- file.path(plots_dir, "wgcna")
if (!dir.exists(wgcna_dir)) dir.create(wgcna_dir, recursive = TRUE)

## ============================================================
## 4) PREPARE EXPRESSION MATRIX
## ============================================================

cat("\n=== PREPARING EXPRESSION DATA FOR WGCNA ===\n")

## DESeq2 normalization (VST)
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = sample_meta,
  design    = ~ 1  ## No design for WGCNA - we want unbiased co-expression
)

## Prefilter: require minimum counts
keep <- rowSums(counts(dds) >= 10) >= ncol(dds) * 0.5
dds <- dds[keep, ]
cat("[INFO] Genes after prefilter: ", nrow(dds), "\n", sep = "")

## VST transform
vst_obj <- vst(dds, blind = TRUE)
expr_mat <- assay(vst_obj)

## Filter to most variable genes (WGCNA works best with variable genes)
gene_vars <- apply(expr_mat, 1, var)
var_threshold <- quantile(gene_vars, probs = wgcna_params$min_variance_quantile)
keep_var <- gene_vars >= var_threshold

## Cap at max_genes
if (sum(keep_var) > wgcna_params$max_genes) {
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:wgcna_params$max_genes]
  expr_mat <- expr_mat[top_var_genes, ]
} else {
  expr_mat <- expr_mat[keep_var, ]
}

cat("[INFO] Genes for WGCNA: ", nrow(expr_mat), "\n", sep = "")

## WGCNA expects samples in rows, genes in columns
datExpr <- t(expr_mat)

## Check for outlier samples
cat("\n=== CHECKING FOR OUTLIER SAMPLES ===\n")
sampleTree <- hclust(dist(datExpr), method = "average")

png(file.path(wgcna_dir, paste0("SampleDendrogram_", run_tag, ".png")),
    width = 12, height = 6, units = "in", res = 300)
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1, cex.main = 1.5)
dev.off()
cat("[OK] Saved sample dendrogram\n")

## Check for genes/samples with too many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("[WARN] Removed problematic samples/genes\n")
}

## ============================================================
## 5) PICK SOFT-THRESHOLDING POWER
## ============================================================

cat("\n=== PICKING SOFT-THRESHOLDING POWER ===\n")
log_time("Starting power estimation")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = wgcna_params$network_type, verbose = 2)

## Auto-select power (first to achieve R^2 > 0.8)
if (is.null(wgcna_params$soft_power)) {
  r2_threshold <- 0.80
  candidates <- sft$fitIndices[sft$fitIndices$SFT.R.sq > r2_threshold, ]
  if (nrow(candidates) > 0) {
    soft_power <- candidates$Power[1]
  } else {
    ## Fallback: use power with highest R^2
    soft_power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
  }
} else {
  soft_power <- wgcna_params$soft_power
}

cat("[INFO] Selected soft-thresholding power: ", soft_power, "\n", sep = "")
cat("[INFO] Scale-free R^2 at this power: ",
    round(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == soft_power], 3), "\n", sep = "")

## Plot scale-free topology fit
png(file.path(wgcna_dir, paste0("SoftThreshold_", run_tag, ".png")),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))

## Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R²)",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "red", lty = 2)

## Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = 0.9, col = "red")

dev.off()
cat("[OK] Saved soft threshold plot\n")

## ============================================================
## 6) NETWORK CONSTRUCTION & MODULE DETECTION
## ============================================================

cat("\n=== BUILDING CO-EXPRESSION NETWORK ===\n")
log_time("Starting network construction")

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = wgcna_params$network_type,
  TOMType = "signed",
  minModuleSize = wgcna_params$min_module_size,
  reassignThreshold = 0,
  mergeCutHeight = wgcna_params$merge_cut_height,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3,
  maxBlockSize = ncol(datExpr)
)

log_time("Network construction complete")

## Convert numeric labels to colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

cat("[INFO] Number of modules: ", length(unique(moduleColors)) - 1, " (excluding grey/unassigned)\n", sep = "")
cat("[INFO] Module sizes:\n")
module_sizes <- sort(table(moduleColors), decreasing = TRUE)
for (mod in names(module_sizes)) {
  cat("  ", mod, ": ", module_sizes[mod], " genes\n", sep = "")
}

## Plot dendrogram with module colors
png(file.path(wgcna_dir, paste0("ModuleDendrogram_", run_tag, ".png")),
    width = 12, height = 8, units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module Colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()
cat("[OK] Saved module dendrogram\n")

## ============================================================
## 7) MODULE-TRAIT CORRELATIONS
## ============================================================

cat("\n=== MODULE-TRAIT CORRELATIONS ===\n")

## Calculate module eigengenes
MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

## Build trait matrix (numeric encoding)
trait_matrix <- data.frame(row.names = rownames(datExpr))

## Encode traits as numeric
if ("Genotype" %in% colnames(sample_meta)) {
  trait_matrix$Genotype <- as.numeric(factor(sample_meta[rownames(datExpr), "Genotype"])) - 1
}
if ("Tissue" %in% colnames(sample_meta)) {
  trait_matrix$Tissue <- as.numeric(factor(sample_meta[rownames(datExpr), "Tissue"])) - 1
}
if ("Sex" %in% colnames(sample_meta)) {
  trait_matrix$Sex <- as.numeric(factor(sample_meta[rownames(datExpr), "Sex"])) - 1
}

## Correlate eigengenes with traits
moduleTraitCor <- cor(MEs, trait_matrix, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("[INFO] Module-trait correlations computed\n")
cat("\nSignificant module-trait associations (p < ", wgcna_params$trait_pvalue_cut, "):\n", sep = "")

for (i in seq_len(nrow(moduleTraitCor))) {
  for (j in seq_len(ncol(moduleTraitCor))) {
    if (moduleTraitPvalue[i, j] < wgcna_params$trait_pvalue_cut) {
      cat("  ", rownames(moduleTraitCor)[i], " ~ ", colnames(moduleTraitCor)[j],
          ": r = ", round(moduleTraitCor[i, j], 3),
          " (p = ", format(moduleTraitPvalue[i, j], digits = 3), ")\n", sep = "")
    }
  }
}

## Save module-trait correlation table
mt_df <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("Module") %>%
  left_join(
    as.data.frame(moduleTraitPvalue) %>%
      rownames_to_column("Module") %>%
      rename_with(~ paste0(.x, "_pvalue"), -Module),
    by = "Module"
  )

mt_file <- paste0("wgcna_module_trait_correlations_", run_tag, ".csv")
write.csv(mt_df, file = file.path(tables_dir, mt_file), row.names = FALSE)
cat("[OK] Saved: ", mt_file, "\n", sep = "")

## Heatmap of module-trait correlations
png(file.path(wgcna_dir, paste0("ModuleTraitHeatmap_", run_tag, ".png")),
    width = 8, height = max(6, nrow(moduleTraitCor) * 0.4), units = "in", res = 300)

## Build text matrix for display
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(trait_matrix),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()
cat("[OK] Saved module-trait heatmap\n")

## ============================================================
## 8) HUB GENE IDENTIFICATION
## ============================================================

cat("\n=== IDENTIFYING HUB GENES ===\n")

## Calculate module membership (kME) for each gene
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
colnames(geneModuleMembership) <- gsub("^ME", "kME_", colnames(geneModuleMembership))

## For each module, find hub genes (highest kME)
hub_genes_all <- list()

for (mod_color in unique(moduleColors)) {
  if (mod_color == "grey") next  ## Skip unassigned

  mod_genes <- colnames(datExpr)[moduleColors == mod_color]
  kme_col <- paste0("kME_", mod_color)

  if (!(kme_col %in% colnames(geneModuleMembership))) next

  mod_kme <- geneModuleMembership[mod_genes, kme_col, drop = FALSE]
  mod_kme$gene <- rownames(mod_kme)
  mod_kme <- mod_kme %>%
    arrange(desc(abs(get(kme_col)))) %>%
    slice_head(n = wgcna_params$hub_n)

  hub_genes_all[[mod_color]] <- data.frame(
    Module = mod_color,
    Gene = mod_kme$gene,
    kME = mod_kme[[kme_col]],
    stringsAsFactors = FALSE
  )

  cat("  ", mod_color, " module: top hub = ", mod_kme$gene[1],
      " (kME = ", round(mod_kme[[kme_col]][1], 3), ")\n", sep = "")
}

hub_df <- bind_rows(hub_genes_all)
hub_file <- paste0("wgcna_hub_genes_", run_tag, ".csv")
write.csv(hub_df, file = file.path(tables_dir, hub_file), row.names = FALSE)
cat("[OK] Saved hub genes: ", hub_file, "\n", sep = "")

## ============================================================
## 9) MODULE ENRICHMENT (GO:BP for genotype-associated modules)
## ============================================================

cat("\n=== MODULE ENRICHMENT ANALYSIS ===\n")

## Find modules significantly associated with Genotype
genotype_modules <- rownames(moduleTraitCor)[
  moduleTraitPvalue[, "Genotype"] < wgcna_params$trait_pvalue_cut
]
genotype_modules <- gsub("^ME", "", genotype_modules)

if (length(genotype_modules) == 0) {
  cat("[INFO] No modules significantly associated with Genotype. Enriching top 3 largest modules.\n")
  non_grey <- names(module_sizes)[names(module_sizes) != "grey"]
  genotype_modules <- head(non_grey, 3)
}

cat("[INFO] Enriching modules: ", paste(genotype_modules, collapse = ", "), "\n", sep = "")

module_enrichment <- list()

for (mod_color in genotype_modules) {
  cat("\n--- Module: ", mod_color, " ---\n", sep = "")

  mod_genes <- colnames(datExpr)[moduleColors == mod_color]
  cat("  Genes: ", length(mod_genes), "\n", sep = "")

  ## Convert to Entrez IDs
  entrez_df <- tryCatch({
    bitr(mod_genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db, drop = TRUE)
  }, error = function(e) {
    data.frame(SYMBOL = character(), ENTREZID = character())
  })

  if (nrow(entrez_df) < 5) {
    cat("  [WARN] Too few genes mapped to Entrez\n")
    next
  }

  ## GO:BP enrichment
  tryCatch({
    ego <- enrichGO(
      gene = entrez_df$ENTREZID,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      readable = TRUE,
      minGSSize = 10,
      maxGSSize = 500
    )

    if (!is.null(ego) && nrow(ego@result) > 0) {
      sig_paths <- ego@result %>% filter(p.adjust < 0.05)
      cat("  GO:BP significant terms: ", nrow(sig_paths), "\n", sep = "")

      if (nrow(sig_paths) > 0) {
        module_enrichment[[mod_color]] <- sig_paths %>%
          mutate(Module = mod_color) %>%
          slice_head(n = 20)

        enrich_file <- paste0("wgcna_enrichment_", mod_color, "_", run_tag, ".csv")
        write.csv(sig_paths, file = file.path(tables_dir, enrich_file), row.names = FALSE)
        cat("  [OK] Saved: ", enrich_file, "\n", sep = "")

        ## Top terms
        cat("  Top 5 GO terms:\n")
        for (k in seq_len(min(5, nrow(sig_paths)))) {
          cat("    ", sig_paths$Description[k], " (p.adj = ",
              format(sig_paths$p.adjust[k], digits = 3), ")\n", sep = "")
        }
      }
    }
  }, error = function(e) {
    cat("  [WARN] GO enrichment failed: ", conditionMessage(e), "\n", sep = "")
  })
}

## ============================================================
## 10) SAVE MODULE ASSIGNMENTS
## ============================================================

cat("\n=== SAVING MODULE ASSIGNMENTS ===\n")

module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

## Add kME for assigned module
module_df$kME <- sapply(seq_len(nrow(module_df)), function(i) {
  gene <- module_df$Gene[i]
  mod <- module_df$Module[i]
  kme_col <- paste0("kME_", mod)
  if (kme_col %in% colnames(geneModuleMembership)) {
    geneModuleMembership[gene, kme_col]
  } else {
    NA
  }
})

module_assign_file <- paste0("wgcna_module_assignments_", run_tag, ".csv")
write.csv(module_df, file = file.path(tables_dir, module_assign_file), row.names = FALSE)
cat("[OK] Saved: ", module_assign_file, "\n", sep = "")

## Save eigengenes
me_df <- as.data.frame(MEs) %>% rownames_to_column("Sample")
me_file <- paste0("wgcna_eigengenes_", run_tag, ".csv")
write.csv(me_df, file = file.path(tables_dir, me_file), row.names = FALSE)
cat("[OK] Saved: ", me_file, "\n", sep = "")

## ============================================================
## 11) SESSION INFO
## ============================================================

session_file <- file.path(outdir, "logs", paste0("sessionInfo_wgcna_", run_tag, ".txt"))
if (dir.exists(dirname(session_file))) {
  sink(session_file)
  cat("=== R SESSION INFORMATION - Part 8 (WGCNA) ===\n\n")
  print(sessionInfo())
  sink()
}

cat("\n=== Part 8 (WGCNA) COMPLETE ===\n")
cat("Modules detected: ", length(unique(moduleColors)) - 1, "\n", sep = "")
cat("Output directory: ", wgcna_dir, "\n", sep = "")
