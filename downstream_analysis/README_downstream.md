# KAT8 downstream analysis — TF activity & cell↔tissue concordance

Two analyses that address **"is the KAT8 effect direct or indirect, and is it
cell-autonomous?"** using **only data you already have**. Both are threshold-free
and aggregate weak signals across many genes, so the **low 3T3-L1 DEG count does
not limit them**.

- `part8_tf_activity.R` — which TF programs shift when KAT8 is lost (decoupleR ULM × CollecTRI).
- `part9_cell_tissue_concordance.R` — which tissue changes reproduce in the pure adipocyte cells.
- `config_downstream.R` — shared config; **edit the input paths here**.

## What to upload (this is all I need)

Just the three per-gene DESeq2 result CSVs your pipeline already writes. Put them
in `downstream_analysis/data/` (or edit `DE_PATHS` in `config_downstream.R`):

| Contrast | File your pipeline produces | Rename to |
|---|---|---|
| Cells | `DE_cells_KAT8KD_vs_CTL_<tag>.csv` | `data/DE_cells_KAT8KD_vs_CTL.csv` |
| iWAT | `DE_tissue_iWAT_KD_vs_CTL_<tag>.csv` | `data/DE_tissue_iWAT_KD_vs_CTL.csv` |
| gWAT | `DE_tissue_gWAT_KD_vs_CTL_<tag>.csv` | `data/DE_tissue_gWAT_KD_vs_CTL.csv` |

Each must contain columns `gene_name`, `log2FoldChange`, `stat`, `pvalue`, `padj`
(your current outputs already do — the first unnamed column is read as the gene name).

> If you also upload `counts.txt` + the cell DESeq2 script, I can additionally apply the
> DESeq2-level rigor fixes below (shrinkage, threshold-testing) and re-emit the DE tables.

## How to run

```r
# from the repo root
source("downstream_analysis/part8_tf_activity.R")            # TF activity
source("downstream_analysis/part9_cell_tissue_concordance.R") # concordance
```

Outputs land in `downstream_analysis/results/`. `part8` needs internet **once** to
fetch the CollecTRI network from OmniPath; it is then cached to
`downstream_analysis/collectri_mouse.rds`.

New packages used: `decoupleR`, `OmnipathR`, `fgsea` (Bioconductor); `pheatmap`,
`ggrepel`, optional `qvalue`. The scripts install anything missing.

## How to read the results

**Organizing hypothesis (from your slides):** KAT8 is an *activator*, so its *direct*
effect is loss of expression. Prediction — the **down** metabolic/adipocyte-identity
program is **direct & cell-autonomous** (reproduces in cells); the **up** inflammatory
program is **indirect & non-cell-autonomous** (infiltration, tissue-only).

- **part8 heatmap** — a TF active in *cells and* tissue is a cell-autonomous regulator;
  a TF active only in tissue (NF-κB/IRF/STAT) is likely infiltration-driven. Watch the
  adipogenic TFs (Pparg, Cebpa, Srebf1, nuclear receptors) in the **cells** column —
  their repression is your cell-autonomous mechanism.
- **part9 scatter + GSEA** — a positive, significant `NES(tissue_DOWN)` in cells means
  the tissue down-program is cell-autonomous; a non-significant `tissue_UP` in cells
  supports the inflammation being non-cell-autonomous. The `binom_p` directional tests
  quantify the same thing while borrowing power from the well-powered tissue DEGs.

---

## Rigor audit of the existing pipeline

Severity: **[H]** materially affects conclusions · **[M]** worth fixing · **[L]** minor.
"Applied" = handled in the new scripts. "Recommended" = a change to your DESeq2 scripts
(`KAT8_bulk_cells_comprehensive.R` / `part1_main_analysis.R`) — I can make these if you
upload counts + scripts.

1. **[H] Remove the literature DEG-count "consistency check."** *Recommended.*
   The cell script flags results "CONSISTENT" if the DEG count lands in 80–150 (vs Li 2012 = 89,
   Ravens 2014 = 120). DEG count depends on sequencing depth, n, and cell system, so matching an
   unrelated study is not a QC criterion — and it invites tuning the |log2FC| cutoff to hit a
   target count, which is circular. Report DEGs at pre-registered thresholds and drop the check.

2. **[H] Use `lfcShrink()` for effect sizes and ranking.** *Recommended.*
   At n = 4/group the raw MLE `log2FoldChange` is inflated for low-count genes. Add
   `lfcShrink(dds, coef="Genotype_KAT8KD_vs_CTL", type="apeglm")` and use the **shrunken** LFC
   for volcano, heatmaps, and the concordance scatter. Keep the **unshrunken Wald `stat`** for
   GSEA/TF ranking (correct as-is). This matters most for the low-power cells.

3. **[M] Replace the post-hoc `|log2FC| >` filter with threshold-testing.** *Recommended.*
   Filtering on a point-estimate LFC after an FDR cut does not control error at that fold-change.
   Prefer `results(dds, lfcThreshold=log2(1.5), altHypothesis="greaterAbs")`, which computes
   p-values against the threshold. For the cells especially, report FDR-only DEGs as the primary
   count and treat |log2FC| as a descriptor, not a gate.

4. **[M] For direct/indirect, do not lean on a hard DEG-overlap (Venn).** *Applied.*
   With few cell DEGs a Venn is underpowered and misleading. `part9` instead uses rank
   correlation, a signature-GSEA of the tissue program against the full cell ranking, and
   directional binomial tests that borrow power from the tissue side — all threshold-free on cells.

5. **[L] PCA on VST.** *Recommended.* The cell script runs `prcomp(t(vst), scale.=TRUE)` on all
   genes; unit-scaling every gene upweights noisy low-variance genes. DESeq2 convention is the
   top ~500 most-variable genes without unit scaling (`DESeq2::plotPCA`). Cosmetic, not a
   conclusion-changer.

6. **[L] Cross-dataset gene matching.** *Applied.* `part9` merges cells↔tissue on gene symbol and
   de-duplicates by keeping the most significant row, so `make.unique` suffixes don't silently
   drop genes.

7. **[L] Permutation calibration.** *Applied.* `part8` can permute gene labels (`n_permutations`)
   to confirm the analytic TF p-values aren't inflated at this sample size.

The validated tissue design (`~ 0 + GroupDepot`, sex handling) is sound and left unchanged.

## What still needs new data (out of scope here, but worth stating in the paper)
- **snRNA-seq** — the definitive test of *which cell type* drives each signature.
- **H4K16ac / KAT8 ChIP-seq** — the definitive *direct-target* map. Until then, integrating
  *public* MOF/H4K16ac ChIP (Cistrome/ENCODE) via BETA is the strongest no-new-experiment proxy
  (a possible part 10).
- KAT8 also has **non-histone substrates**, so some effects are post-translational and invisible
  to any RNA-based readout — state as a limitation.
