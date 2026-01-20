# KAT8 Knockdown RNA-seq Analysis Pipeline

Comprehensive RNA-seq differential expression and pathway enrichment analysis pipeline for KAT8 knockdown in adipose tissue (iWAT and gWAT, male and female mice).

## Overview

This pipeline performs three main analyses:

1. **Part 1: DESeq2 Differential Expression** (`part1_main_analysis.R`)
   - Quality control filtering (count ≥10 in ≥3 samples per group)
   - Differential expression analysis using DESeq2
   - Visualization: PCA, heatmaps, volcano plots, MA plots
   - Gene order locked across all heatmaps for cross-comparison

2. **Part 2: Over-Representation Analysis** (`part2_ora.R`)
   - ORA for GO:BP and KEGG pathways
   - Separate analysis for up- and down-regulated genes
   - GO term simplification to reduce redundancy
   - Statistical background: all DESeq2-tested genes

3. **Part 3: Gene Set Enrichment Analysis** (`part3_fgsea.R`)
   - FGSEA for GO:BP, KEGG, WikiPathways, and MSigDB Hallmark
   - Ranked by DESeq2 Wald statistic
   - GO:BP simplification using semantic similarity
   - Leading edge gene export

## Requirements

### Software
- **R version**: 4.0 or higher (recommended: 4.3+)
- **RAM**: Minimum 8 GB (16 GB recommended for large datasets)

### R Packages

#### Part 1 (DESeq2 Analysis)
```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer"))
BiocManager::install(c("DESeq2", "ComplexHeatmap"))
```

#### Part 2 (ORA Analysis)
```r
install.packages("dplyr")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot"))
```

#### Part 3 (FGSEA Analysis)
```r
install.packages(c("dplyr", "ggplot2"))
BiocManager::install(c("fgsea", "clusterProfiler", "org.Mm.eg.db", "GOSemSim", "msigdbr"))
```

Install BiocManager first if not already installed:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

## Input Files

### Required Directory Structure
```
KAT8-Project/
├── data/
│   ├── counts.csv           # Raw count matrix (genes × samples)
│   └── metadata.csv         # Sample metadata
└── output/                  # Created automatically
    ├── plots/
    ├── tables/
    └── logs/
```

### Input File Formats

**counts.csv**: Gene count matrix
- First column: Gene symbols (as row names)
- Remaining columns: Sample counts
- Example:
  ```
  Gene,Sample1,Sample2,Sample3,...
  Gapdh,5432,6123,5890,...
  Actb,8901,9234,8765,...
  ```

**metadata.csv**: Sample metadata
- Required columns:
  - `Sample`: Sample identifier (must match count matrix column names)
  - `Tissue`: Depot (iWAT or gWAT)
  - `Sex`: Sex (M or F)
  - `Genotype`: Genotype (CTL or KD)
- Example:
  ```
  Sample,Tissue,Sex,Genotype
  Sample1,iWAT,F,CTL
  Sample2,iWAT,F,KD
  Sample3,gWAT,M,CTL
  ```

## Usage

### Quick Start (Run Full Pipeline)

Run all three parts in sequence:

```bash
Rscript run_all.R
```

or from R console:

```r
source("run_all.R")
```

The master script will:
- Run all parts in order
- Handle errors gracefully
- Track timing for each part
- Save a timing log to `output/logs/`

### Run Individual Parts

If you need to run parts separately:

```r
# Part 1: DESeq2 analysis (run first)
source("part1_main_analysis.R")

# Part 2: ORA analysis (requires Part 1 output)
source("part2_ora.R")

# Part 3: FGSEA analysis (requires Part 1 output)
source("part3_fgsea.R")
```

**Note**: Part 2 and Part 3 require Part 1 to complete successfully first, as they use the differential expression results.

## Output Structure

### Part 1: DESeq2 Analysis

**Tables** (`output/tables/`):
- `DE_tissue_<contrast>_<timestamp>.csv` - Differential expression results for each contrast
- `DE_all_contrasts_<timestamp>.csv` - Combined results for all contrasts

**Plots** (`output/plots/`):
- `PCA_<timestamp>.png` - Principal component analysis
- `GOI_heatmap_<contrast>_<timestamp>.png` - Genes of interest heatmap (per contrast)
- `GOI_heatmap_OVERALL_<timestamp>.png` - Combined heatmap across all contrasts
- `Volcano_<contrast>_<timestamp>.png` - Volcano plots
- `MA_plot_<contrast>_<timestamp>.png` - MA plots

**Logs** (`output/logs/`):
- `sessionInfo_<timestamp>.txt` - R session info for reproducibility

### Part 2: ORA Analysis

**Tables** (`output/tables/`):
- `ORA_gobp_<contrast>_<direction>_<timestamp>.csv` - GO:BP enrichment results
- `ORA_kegg_<contrast>_<direction>_<timestamp>.csv` - KEGG enrichment results
- `ORA_sanity_check_<timestamp>.csv` - Quality control summary

**Plots** (`output/plots/`):
- `ORA_plot_gobp_<contrast>_<direction>_<timestamp>.png` - GO:BP bar plots
- `ORA_plot_kegg_<contrast>_<direction>_<timestamp>.png` - KEGG bar plots

**Logs** (`output/logs/`):
- `sessionInfo_ORA_<timestamp>.txt` - R session info

### Part 3: FGSEA Analysis

**Tables** (`output/tables/`):
- `fgsea_gobp_<contrast>_<timestamp>.csv` - GO:BP FGSEA results
- `fgsea_kegg_<contrast>_<timestamp>.csv` - KEGG FGSEA results
- `fgsea_wikipathways_<contrast>_<timestamp>.csv` - WikiPathways results
- `fgsea_hallmark_<contrast>_<timestamp>.csv` - MSigDB Hallmark results
- `fgsea_sanity_check_<timestamp>.csv` - Quality control summary

**Plots** (`output/plots/`):
- `fgsea_plot_<database>_<contrast>_<timestamp>.png` - NES bar plots for each database

**Logs** (`output/logs/`):
- `sessionInfo_fgsea_<timestamp>.txt` - R session info

## Key Parameters

### Part 1: QC Filtering
- **Minimum count**: 10 reads per gene
- **Minimum samples**: 3 samples per group
- **DEG cutoffs**: FDR < 0.05, |log2FC| > 1

### Part 2: ORA
- **Gene universe**: All DESeq2-tested genes (conservative)
- **p-value cutoff**: 0.05
- **q-value cutoff**: 0.1
- **GO simplification**: Enabled (semantic similarity cutoff: 0.7)

### Part 3: FGSEA
- **Ranking metric**: DESeq2 Wald statistic
- **FDR cutoff**: 0.05
- **GO simplification**: Enabled (Wang method, cutoff: 0.7)
- **Min gene set size**: 10 genes
- **Max gene set size**: 500 genes

## Contrasts Analyzed

The pipeline analyzes four contrasts (KD vs CTL):

1. **iWAT_F_KD_vs_CTL** - iWAT female KD vs control
2. **iWAT_M_KD_vs_CTL** - iWAT male KD vs control
3. **gWAT_F_KD_vs_CTL** - gWAT female KD vs control
4. **gWAT_M_KD_vs_CTL** - gWAT male KD vs control

**Note**: Gene order in heatmaps is locked to iWAT_F reference for cross-comparison.

## Troubleshooting

### Low Number of Pathways Found

If you see few enriched pathways (especially for down-regulated genes in gWAT):

1. **Check DEG counts**: Verify you have sufficient DEGs
   - Run diagnostic: `source("diagnostic_gene_mapping.R")`
   - Recommended: ≥50 DEGs per direction

2. **Check gene mapping rate**:
   - Expected mapping rate: 80-90%
   - Unmapped genes are often non-coding (Gm*, *Rik) and lack Entrez IDs

3. **Biological interpretation**:
   - Some tissues naturally show scattered regulation (not organized into pathways)
   - This is not a technical failure but biological reality

### Memory Issues

If R crashes or runs out of memory:

1. Close other applications
2. Increase R memory limit:
   ```r
   memory.limit(size = 16000)  # Windows only
   ```
3. Run parts separately instead of using `run_all.R`
4. Consider using an HPC cluster for very large datasets

### Package Installation Errors

If Bioconductor packages fail to install:

1. Update R to the latest version
2. Update Bioconductor:
   ```r
   BiocManager::install(version = "3.18")
   ```
3. Install dependencies manually if needed

### File Not Found Errors

Ensure your working directory is set to the project root:

```r
setwd("/path/to/KAT8-Project")
getwd()  # Verify
```

## Reproducibility

Each script saves session information including:
- R version
- Package versions
- Operating system details

Session info files are saved to `output/logs/` for reproducibility.

## Citation

If you use this pipeline in your research, please cite:

- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014) Genome Biology
- **clusterProfiler**: Yu, G., et al. (2012) OMICS
- **fgsea**: Korotkevich, G., et al. (2021) bioRxiv

## Contact

For questions or issues with this pipeline, please open an issue on the project repository or contact the KAT8 project team.

## License

This project is licensed under the MIT License.

---

**Last Updated**: 2026-01-20
**Pipeline Version**: 1.0
