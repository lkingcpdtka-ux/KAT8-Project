# KAT8 RNA-seq Analysis Workflow

This repository contains a three-part analysis pipeline for KAT8 bulk RNA-seq data (iWAT and gWAT tissue).

## Overview

The analysis has been split into three standalone scripts:

1. **Part 1: Main Gene-Level Analysis** (`part1_main_analysis.R`)
2. **Part 2: ORA Pathway Analysis** (`part2_ora.R`)
3. **Part 3: fgsea Pathway Analysis** (`part3_fgsea.R`)

## Why Split Into Three Parts?

The original monolithic script had several issues:
- **Error in enrichGO**: Missing universe/background causing crashes
- **KEGG subcategory bug**: msigdbr filtering was incorrect
- **Difficult debugging**: Pathway errors would halt the entire pipeline
- **Long runtime**: All analyses had to complete in one session

The split approach provides:
- ✅ **Modularity**: Run only the analyses you need
- ✅ **Error isolation**: Pathway analysis errors don't affect gene-level results
- ✅ **Faster iteration**: Re-run pathway analysis without repeating DESeq2
- ✅ **Better fixes**: Each part has targeted fixes for known issues

## Part 1: Main Gene-Level Analysis

**File**: `part1_main_analysis.R`

**What it does**:
- Quality control (density plots, PCA, MDS)
- DESeq2 differential expression analysis
- Volcano plots with top DEGs labeled
- Heatmaps of top DEGs (ComplexHeatmap)
- Genes of interest analysis (focused heatmaps and volcanos)
- Overall contrast (all depots/sexes combined)

**What it does NOT do**:
- Pathway enrichment (ORA or fgsea)

**Outputs**:
- `savepoints/RUN_YYYYMMDD_HHMMSS/tables/DE_tissue_*.csv` - Differential expression results
- `savepoints/RUN_YYYYMMDD_HHMMSS/plots/Volcano_*.png` - Volcano plots
- `savepoints/RUN_YYYYMMDD_HHMMSS/plots/Heatmap_*.png` - Heatmaps
- `savepoints/RUN_YYYYMMDD_HHMMSS/plots/PCA_*.png` - PCA plot
- `savepoints/RUN_YYYYMMDD_HHMMSS/plots/MDS_*.png` - MDS plot

**Runtime**: ~5-10 minutes

**Usage**:
```r
source("part1_main_analysis.R")
```

## Part 2: ORA Pathway Analysis

**File**: `part2_ora.R`

**What it does**:
- Reads DE results from Part 1
- Performs Over-Representation Analysis (ORA) for:
  - GO Biological Processes (GO:BP)
  - KEGG pathways
- Creates pathway bar plots (direction-specific colors)
- Generates sanity check reports

**Key Fixes Applied**:
- ✅ **Proper universe/background**: Uses all tested genes as background (fixes enrichGO error)
- ✅ **Better error handling**: Won't crash if one contrast fails
- ✅ **Improved logging**: Tracks gene mapping and pathway counts

**Outputs**:
- `results_ora/RUN_YYYYMMDD_HHMMSS/tables/ORA_*.csv` - Pathway enrichment tables
- `results_ora/RUN_YYYYMMDD_HHMMSS/plots/ORA_barplot_*.png` - Pathway bar plots
- `results_ora/RUN_YYYYMMDD_HHMMSS/tables/ORA_sanity_check_*.csv` - QC report

**Runtime**: ~2-5 minutes

**Usage**:
```r
# Run AFTER part1_main_analysis.R
source("part2_ora.R")
```

## Part 3: fgsea Pathway Analysis

**File**: `part3_fgsea.R`

**What it does**:
- Reads DE results from Part 1
- Performs Gene Set Enrichment Analysis (fgsea) for:
  - GO Biological Processes (GO:BP)
  - KEGG pathways
- Uses DESeq2 Wald statistic for gene ranking
- Creates direction-specific pathway plots (NES-based)

**Key Fixes Applied**:
- ✅ **Correct msigdbr filtering**: Fixed KEGG subcategory issue
- ✅ **Wald statistic ranking**: Uses proper DESeq2 test statistic (not p-value-based)
- ✅ **Fallback for KEGG**: Uses org.Mm.eg.db if msigdbr fails

**Outputs**:
- `results_fgsea/RUN_YYYYMMDD_HHMMSS/tables/fgsea_*.csv` - fgsea results
- `results_fgsea/RUN_YYYYMMDD_HHMMSS/plots/fgsea_plot_*_Up_*.png` - Up-regulated pathways
- `results_fgsea/RUN_YYYYMMDD_HHMMSS/plots/fgsea_plot_*_Down_*.png` - Down-regulated pathways
- `results_fgsea/RUN_YYYYMMDD_HHMMSS/tables/fgsea_sanity_check_*.csv` - QC report

**Runtime**: ~5-10 minutes

**Usage**:
```r
# Run AFTER part1_main_analysis.R
source("part3_fgsea.R")
```

## Complete Workflow

### Step 1: Run Part 1 (Required)
```r
source("part1_main_analysis.R")
```

Wait for completion. This generates all gene-level analyses and DE tables.

### Step 2: Run Part 2 (Optional - ORA)
```r
source("part2_ora.R")
```

This reads the most recent Part 1 results and performs ORA enrichment.

### Step 3: Run Part 3 (Optional - fgsea)
```r
source("part3_fgsea.R")
```

This reads the most recent Part 1 results and performs fgsea enrichment.

**Note**: Parts 2 and 3 can be run in any order or independently. They both read from Part 1 outputs.

## Contrasts Analyzed

- `iWAT_F_KD_vs_CTL` - iWAT Female: KAT8KD vs Control
- `iWAT_M_KD_vs_CTL` - iWAT Male: KAT8KD vs Control
- `gWAT_F_KD_vs_CTL` - gWAT Female: KAT8KD vs Control
- `gWAT_M_KD_vs_CTL` - gWAT Male: KAT8KD vs Control
- `OVERALL_KD_vs_CTL` - All depots/sexes: KAT8KD vs Control

## Thresholds

- **Log2 Fold Change**: |logFC| > 1
- **FDR**: < 0.05
- **Prefiltering**: Genes with counts ≥10 in at least min(group size) samples

## Required Files

- `counts.txt` - Raw count matrix (genes × samples)
- `save_core/` - Utility functions for run management

## Genes of Interest

The analysis focuses on:
- **Collagens/ECM**: Col4a1, Col4a2, Col15a1, Col5a3, Col6a6, etc.
- **Adipocyte markers**: Lep, Ppargc1a, Sorbs1, Srebf1, Ppara, Pparg, etc.
- **Matrix remodeling**: Fn1, Mmp3, Timp4, Mmp12, Mmp14, Mmp16

## Troubleshooting

### "No Part 1 run directories found"
Run `part1_main_analysis.R` first. Parts 2 and 3 depend on its outputs.

### "Cannot find stat column"
Your Part 1 DE tables may be from an older version. Re-run Part 1.

### ORA enrichGO error persists
Check that you're using the latest `part2_ora.R` with universe parameter.

### KEGG pathways not found
Part 3 will automatically fall back to org.Mm.eg.db if msigdbr fails.

## Color Scheme

All visualizations use a consistent color scheme:
- **Up-regulated**: Orange (#E69F00)
- **Down-regulated**: Teal/Blue (#0072B2)
- **Not significant**: Grey (#grey70)

## Output Organization

```
KAT8-Project/
├── part1_main_analysis.R
├── part2_ora.R
├── part3_fgsea.R
├── savepoints/
│   └── RUN_YYYYMMDD_HHMMSS/     # Part 1 outputs
│       ├── tables/
│       ├── plots/
│       └── logs/
├── results_ora/
│   └── RUN_YYYYMMDD_HHMMSS/     # Part 2 outputs
│       ├── tables/
│       ├── plots/
│       └── logs/
└── results_fgsea/
    └── RUN_YYYYMMDD_HHMMSS/     # Part 3 outputs
        ├── tables/
        ├── plots/
        └── logs/
```

## Version History

### v2.0 (Current - Split Version)
- Split into three independent scripts
- Fixed enrichGO error with proper universe
- Fixed msigdbr KEGG subcategory issue
- Improved error handling and logging
- Added sanity check reports

### v1.0 (Original - Monolithic)
- Single script with all analyses
- Known issues: enrichGO crashes, KEGG errors

## Questions?

Contact: [Your contact info]
