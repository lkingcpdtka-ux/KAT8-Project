# KAT8 RNA-seq Analysis Workflow

This repository contains a six-part analysis pipeline for KAT8 bulk RNA-seq data (iWAT and gWAT tissue).

## Overview

The analysis has been split into modular standalone scripts:

1. **Part 1: Main Gene-Level Analysis** (`part1_main_analysis.R`)
2. **Part 2: ORA Pathway Analysis** (`part2_ora.R`)
3. **Part 3: fGSEA Pathway Analysis** (`part3_fgsea.R`)
4. **Part 4: Publication Plots** (`part4_publication_plots.R`)
5. **Part 5: Barcode Plots** (`part5_barcode_plots.R`)
6. **Part 6: Biological Interpretation** (`part6_biological_interpretation.R`) *NEW*

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

**Outputs** (all in `savepoints/RUN_YYYYMMDD_HHMMSS/`):
- `tables/DE_tissue_*.csv` - Differential expression results
- `plots/Volcano_*.png` - Volcano plots (including overall contrast)
- `plots/Heatmap_*.png` - Heatmaps (including overall contrast)
- `plots/Heatmap_GOI_*.png` - Genes of interest heatmaps
- `plots/EnhancedVolcano_GOI_*.png` - Genes of interest volcano plots
- `plots/PCA_*.png` - PCA plot
- `plots/MDS_*.png` - MDS plot

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

**Outputs** (added to SAME directory as Part 1):
- `tables/ORA_*.csv` - Pathway enrichment tables
- `plots/ORA_barplot_*.png` - Pathway bar plots
- `tables/ORA_sanity_check_*.csv` - QC report

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

**Outputs** (added to SAME directory as Part 1):
- `tables/fgsea_*.csv` - fgsea results
- `plots/fgsea_plot_*_Up_*.png` - Up-regulated pathways
- `plots/fgsea_plot_*_Down_*.png` - Down-regulated pathways
- `tables/fgsea_sanity_check_*.csv` - QC report

**Runtime**: ~5-10 minutes

**Usage**:
```r
# Run AFTER part1_main_analysis.R
source("part3_fgsea.R")
```

## Part 6: Biological Interpretation

**File**: `part6_biological_interpretation.R`

**What it does**:
- **Pathway Convergence Analysis**: Identifies pathways significant in BOTH ORA and fGSEA (highest confidence)
- **Cell Type Marker Overlay**: Scores enrichment of adipose-specific cell type signatures (adipocytes, macrophages, T cells, endothelial, fibroblasts, etc.)
- **Leading Edge × Cell Type Integration**: Identifies which cell types' markers appear in pathway leading edges
- **Publication Figures**: Generates manuscript-ready visualizations

**Key Features**:
- ✅ **14 cell type signatures** curated for adipose tissue biology
- ✅ **Pathway theme categorization**: Auto-categorizes pathways into Inflammation, Lipid Metabolism, ECM Remodeling, etc.
- ✅ **Fisher's exact tests**: Statistical testing for marker enrichment in DEGs
- ✅ **Interpretation guidelines**: Built-in guidance for manuscript language

**Important Caveat**:
Bulk RNA-seq CANNOT distinguish cell composition from cell state. The script explicitly guides appropriate language:
- **GOOD**: "Macrophage marker genes were enriched among upregulated genes"
- **BAD**: "Macrophage abundance was increased"

**Outputs** (added to SAME directory as Part 1):
- `tables/convergent_pathways_*.csv` - High-confidence pathways (ORA ∩ fGSEA)
- `tables/cell_type_signatures_*.csv` - Cell type marker enrichment scores
- `tables/leading_edge_celltype_*.csv` - Pathway × cell type integration
- `plots/interpretation/Fig_celltype_signatures_*.png` - Heatmap of cell type scores
- `plots/interpretation/Fig_pathway_convergence_*.png` - ORA vs fGSEA agreement
- `plots/interpretation/Fig_pathway_celltype_*.png` - Pathway × cell type tile plot

**Runtime**: ~2-5 minutes

**Usage**:
```r
# Run AFTER part1, part2, and part3
source("part6_biological_interpretation.R")
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

### Step 4: Run Part 6 (Optional - Biological Interpretation)
```r
source("part6_biological_interpretation.R")
```

This integrates ORA and fGSEA results, scores cell type signatures, and generates interpretation-ready figures. **Requires Parts 1, 2, and 3 to be completed first.**

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

## Frequently Asked Questions

### Why do we need to re-run DESeq2 for the overall contrast?

**Short answer**: Different design matrices require separate DESeq2 objects.

**Detailed explanation**:
- **Per-stratum contrasts** use design `~ 0 + GroupFull` where GroupFull has 8 levels (iWAT_F_CTL, iWAT_F_KAT8KD, iWAT_M_CTL, etc.). This allows us to make specific comparisons within each depot/sex combination.
- **Overall contrast** uses design `~ Depot + Sex + Genotype` to test for the main effect of genotype while controlling for depot and sex as covariates.

These are fundamentally different statistical models, so we need to create a new DESeq2 object with the appropriate design for each analysis. This is standard practice in DESeq2.

### Why does the overall contrast have a heatmap but not genes-of-interest plots?

The overall contrast shows ALL samples across depots and sexes, so:
- ✅ **Heatmap works**: Shows top DEGs with depot/sex/genotype annotations
- ❌ **GOI plots would be misleading**: Genes of interest are specific to certain depot/sex combinations, so an "overall" view wouldn't be biologically meaningful

## Troubleshooting

### "No Part 1 run directories found"
Run `part1_main_analysis.R` first. Parts 2 and 3 depend on its outputs.

### "Cannot find stat column"
Your Part 1 DE tables may be from an older version. Re-run Part 1.

### "could not find function save_run_file"
This is now fixed in the latest version. The function is wrapped in error handling.

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

**All three parts now use the SAME run directory!** This keeps all results from one analysis session together.

```
KAT8-Project/
├── part1_main_analysis.R
├── part2_ora.R
├── part3_fgsea.R
├── part4_publication_plots.R
├── part5_barcode_plots.R
├── part6_biological_interpretation.R
└── savepoints/
    └── RUN_YYYYMMDD_HHMMSS/     # Single directory for all outputs
        ├── tables/
        │   ├── DE_tissue_*.csv                  # Part 1: DE results
        │   ├── DEG_summary_*.csv                # Part 1: Summary
        │   ├── ORA_*.csv                        # Part 2: ORA results
        │   ├── ORA_sanity_check_*.csv           # Part 2: QC
        │   ├── fgsea_*.csv                      # Part 3: fgsea results
        │   ├── fgsea_sanity_check_*.csv         # Part 3: QC
        │   ├── convergent_pathways_*.csv        # Part 6: ORA ∩ fGSEA
        │   ├── cell_type_signatures_*.csv       # Part 6: Cell type scores
        │   └── leading_edge_celltype_*.csv      # Part 6: LE × cell type
        ├── plots/
        │   ├── Volcano_*.png                    # Part 1: Volcano plots
        │   ├── Heatmap_*.png                    # Part 1: Heatmaps
        │   ├── ORA_barplot_*.png                # Part 2: ORA bar plots
        │   ├── fgsea_plot_*.png                 # Part 3: fgsea plots
        │   └── interpretation/                  # Part 6: Interpretation figs
        │       ├── Fig_celltype_signatures_*.png
        │       ├── Fig_pathway_convergence_*.png
        │       └── Fig_pathway_celltype_*.png
        └── logs/
            ├── outliers_*.txt                   # Part 1: Outlier rationale
            └── params_*.txt                     # Part 1: Parameters
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
