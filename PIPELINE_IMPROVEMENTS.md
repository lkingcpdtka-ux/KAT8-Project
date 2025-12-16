# KAT8 Bulk RNA-seq Pipeline - Improvements & Documentation

## Overview

Your bulk RNA-seq pipeline has been significantly enhanced to be **publication-ready** with comprehensive pathway analysis, improved visualizations, and better reproducibility features.

---

## 🎯 Key Improvements Made

### 1. **Pathway Analysis (NEW)**
Added comprehensive pathway enrichment analysis for up- and down-regulated genes in each contrast:

- **Hallmark Gene Sets** (MSigDB): Canonical pathways and biological states
- **KEGG Pathways**: Metabolic and signaling pathways
- **GO Biological Process**: Functional annotations (simplified to reduce redundancy)

**Features:**
- Separate analysis for upregulated and downregulated genes
- Publication-ready bar plots showing:
  - Top 15 pathways by adjusted p-value
  - Gene ratio (proportion of query genes in pathway)
  - Color-coded by significance
  - Clean, professional styling
- CSV tables with full enrichment statistics
- Automatic handling of small gene sets (minimum 5 genes)

### 2. **Publication-Ready Figures**
All plots now meet publication standards:
- **300 DPI resolution** (upgraded from 150 DPI)
- Consistent font sizes and styling
- Professional color schemes (viridis palettes)
- Bold, clear axis labels
- High contrast for print/digital media
- White backgrounds for compatibility

### 3. **Enhanced Robustness**
- **File existence checks**: Verifies `counts.txt` before loading
- **save_core fallback**: Works with or without save_core utilities
- **Better error handling**: Graceful degradation when pathways unavailable
- **Input validation**: Checks for sample columns, missing data
- **Outlier documentation**: Clear notes on why samples removed

### 4. **Improved Reproducibility**
- **sessionInfo()**: Saved to logs folder with:
  - R version and package versions
  - Analysis parameters (thresholds, cutoffs)
  - Run timestamp and identifiers
  - Outlier samples removed
- **Comprehensive logging**: Better tracking of analysis steps
- **Structured output**: Organized plots/tables/logs directories

### 5. **Better DEG Tracking**
Enhanced summary table now includes:
- Total DEGs (FDR < 0.05)
- Number of upregulated genes
- Number of downregulated genes
- Per-contrast breakdown

---

## 📊 Output Structure

```
output/
└── [run_tag]/
    ├── plots/
    │   ├── logCPM_density_*.png
    │   ├── PCA_*.png
    │   ├── MDS_*.png
    │   ├── Volcano_tissue_*.png
    │   ├── Heatmap_tissue_*.png
    │   ├── Pathway_barplot_hallmark_*_Up.png
    │   ├── Pathway_barplot_hallmark_*_Down.png
    │   ├── Pathway_barplot_kegg_*_Up.png
    │   ├── Pathway_barplot_kegg_*_Down.png
    │   ├── Pathway_barplot_gobp_*_Up.png
    │   └── Pathway_barplot_gobp_*_Down.png
    ├── tables/
    │   ├── DE_tissue_*.csv
    │   ├── DEG_summary_tissue_*.csv
    │   ├── Pathway_hallmark_*_Up.csv
    │   ├── Pathway_hallmark_*_Down.csv
    │   ├── Pathway_kegg_*_Up.csv
    │   ├── Pathway_kegg_*_Down.csv
    │   ├── Pathway_gobp_*_Up.csv
    │   └── Pathway_gobp_*_Down.csv
    └── logs/
        └── sessionInfo_*.txt
```

---

## 🔧 Required Packages

### CRAN Packages
- dplyr
- limma
- edgeR
- ggplot2
- RColorBrewer
- ggrepel
- viridis
- pheatmap
- grid

### Bioconductor Packages (NEW)
- clusterProfiler
- org.Mm.eg.db (mouse annotation)
- enrichplot
- DOSE
- pathview
- msigdbr

**Note:** The script will automatically install missing packages.

---

## 📋 How to Use

### Prerequisites
1. Place your `counts.txt` file in the working directory
2. Ensure save_core utilities exist (optional - has fallback)
3. Have internet connection for first run (downloads gene sets)

### Run the Pipeline
```R
# Option 1: From command line
Rscript TEST

# Option 2: From RStudio
source("TEST")
```

### Customize Analysis
Edit these parameters in the script:

```R
## Line ~547-548: Differential expression thresholds
logFC_cut_tissue <- 1      # Log2 fold-change cutoff
fdr_cut_tissue   <- 0.05   # FDR cutoff

## Line ~205-206: Outlier samples
outlier_samples <- c("JS_08", "JS_28")

## Line ~568: Pathway analysis FDR
fdr_cutoff = 0.05  # In run_pathway_analysis function
```

---

## 🐛 Issues Fixed

1. ✅ **Missing pathway analysis** - Now includes Hallmark, KEGG, GO:BP
2. ✅ **save_core dependencies** - Added fallback for missing utilities
3. ✅ **No file checks** - Added validation for counts.txt
4. ✅ **Undocumented outliers** - Added explanatory comments
5. ✅ **Missing sessionInfo** - Now saved to logs folder
6. ✅ **Low resolution plots** - Upgraded to 300 DPI
7. ✅ **No reproducibility info** - Full session info saved
8. ✅ **Inconsistent styling** - Standardized theme across plots
9. ✅ **Poor error messages** - Added informative warnings/errors
10. ✅ **Limited QC reporting** - Enhanced console output

---

## 📈 Typical Workflow Results

For each of the 4 contrasts (iWAT_F, iWAT_M, gWAT_F, gWAT_M):

1. **Differential Expression Table** (CSV)
   - All genes with statistics
   - logFC, p-value, FDR
   - Gene names and IDs

2. **Volcano Plot** (PNG)
   - Color-coded by significance
   - Top DEGs labeled
   - Clear threshold lines

3. **Heatmap** (PNG)
   - Top DEGs clustered
   - Sample annotations
   - Row-scaled expression

4. **Pathway Analysis** (6 bar plots + 6 CSV tables)
   - Upregulated pathways (Hallmark, KEGG, GO)
   - Downregulated pathways (Hallmark, KEGG, GO)
   - Top 15 pathways per database
   - Full statistics tables

---

## 🎨 Figure Customization

### Pathway Bar Plots
Modify the `run_pathway_analysis` function (~line 694):
```R
## Change number of pathways shown
plot_data <- sig_results %>%
  dplyr::slice_head(n = 15)  # Change to 10, 20, etc.

## Change color scheme
scale_fill_viridis(
  option = "plasma"  # Try: viridis, magma, inferno, cividis
)

## Adjust plot dimensions
ggsave(..., width = 10, height = 8)  # Customize size
```

### Volcano Plots
Already optimized, but you can modify:
- Point sizes (line ~864: `size = 2`)
- Label count (line ~822-823: `top_n_fdr`, `top_n_fc`)
- Color palette (line ~866-869)

---

## 🔬 Best Practices for Publication

1. **Methods Section Text:**
   ```
   Differential expression analysis was performed using limma-voom with
   TMM normalization. Genes with FDR < 0.05 and |log2FC| > 1 were
   considered significantly differentially expressed. Pathway enrichment
   analysis was conducted using clusterProfiler against Hallmark, KEGG,
   and GO Biological Process gene sets (FDR < 0.05).
   ```

2. **Always include sessionInfo** in supplementary materials

3. **Report filtering criteria:**
   - CPM > 1 in at least 50% of samples
   - Outlier samples removed (document reasoning)

4. **Cite packages:**
   - limma: Ritchie et al., 2015
   - edgeR: Robinson et al., 2010
   - clusterProfiler: Yu et al., 2012
   - MSigDB: Liberzon et al., 2015

---

## ⚠️ Important Notes

### Outlier Samples
Currently removes `JS_08` and `JS_28`. Update line ~205 if:
- You have different outliers
- You want to include all samples
- You need to document specific QC failures

### Gene ID Conversion
The pipeline converts gene symbols → Entrez IDs for pathway analysis:
- Some genes may not map (reported in console)
- Check warnings if conversion rate is low (<80%)
- Ensure gene names are standard HUGO symbols

### Memory Requirements
For large datasets (>50K genes, >100 samples):
- Increase R memory: `R --max-mem-size=8G`
- Consider filtering more aggressively
- Reduce heatmap gene count (line ~931-933)

### Internet Connection
First run downloads MSigDB gene sets (~50MB):
- Requires active internet connection
- Cached for subsequent runs
- Can take 2-5 minutes

---

## 🚀 Future Enhancements (Optional)

Consider adding:
1. **GSEA** (Gene Set Enrichment Analysis) for ranked gene lists
2. **Sample correlation heatmap** for QC
3. **RLE plots** (Relative Log Expression)
4. **Batch effect visualization** if applicable
5. **Multi-panel figures** for manuscripts
6. **Interactive volcano plots** (e.g., plotly)
7. **Venn diagrams** comparing contrasts
8. **Leading edge analysis** for core enrichment genes

---

## 📞 Troubleshooting

### Error: "counts.txt not found"
- Ensure file is in working directory: `getwd()`
- Check file name (case-sensitive on Linux/Mac)
- Use absolute path if needed: `counts_file <- "/full/path/to/counts.txt"`

### Warning: "Too few genes for pathway analysis"
- Check your thresholds (logFC_cut, fdr_cut)
- Normal for some contrasts with subtle effects
- Consider lowering FDR cutoff for exploratory analysis

### Error: "Cannot convert genes to Entrez ID"
- Check internet connection
- Verify gene names are standard symbols
- Try: `bitr(head(gene_list), "SYMBOL", "ENTREZID", org.Mm.eg.db)`

### Plots look cramped
- Increase figure dimensions in ggsave()
- Reduce number of pathways shown (top 10 instead of 15)
- Adjust font sizes in theme()

---

## ✅ Quality Checklist

Before submitting for publication:
- [ ] All plots at 300+ DPI
- [ ] Consistent color schemes
- [ ] Clear axis labels and legends
- [ ] sessionInfo saved and reviewed
- [ ] Outlier removal documented
- [ ] Methods section written
- [ ] Package citations included
- [ ] Supplementary tables formatted
- [ ] Figure legends written
- [ ] All contrasts analyzed
- [ ] Pathway results validated

---

## 📚 Additional Resources

- **limma User's Guide**: https://bioconductor.org/packages/limma/
- **clusterProfiler book**: https://yulab-smu.top/biomedical-knowledge-mining-book/
- **RNA-seq best practices**: https://www.rna-seqblog.com/
- **MSigDB database**: https://www.gsea-msigdb.org/

---

**Version:** Enhanced 2024
**Author:** Updated pipeline with comprehensive pathway analysis
**License:** Same as original project
