# KAT8 downstream analysis — TF activity & cell↔tissue concordance

Two analyses that address **"is the KAT8 effect direct or indirect, and is it
cell-autonomous?"** using **only data you already have**. Both are threshold-free
and aggregate weak signals across many genes, so the **low 3T3-L1 DEG count does
not limit them**.

- `part4_tf_activity.R` — **cells-focused** TF activity inference (decoupleR ULM × CollecTRI
  on the genome-wide Wald statistic). Centres the 3T3-L1 cells, uses iWAT/gWAT as comparison,
  and prints where the adipogenic vs inflammatory TFs land. **Needs OmniPath reachable once**
  to fetch CollecTRI (then cached to `collectri_mouse.rds`); run it where you have internet.
- `part4b_cell_tissue_concordance.R` — which tissue changes reproduce in the pure adipocyte cells.
  (Also has a runnable Python port: `run_concordance_python.py` + `run_program_level.py`, already executed —
  results in `results/`.)
- `part5_cells_mechanism.R` — **the mechanism analysis for the cells.** Modules:
  **A** secretome (what the KAT8-null adipocyte broadcasts + whether it replicates in vivo),
  **B/C/D** three independent TF layers (TF genes → regulon activity → ChIP-derived target
  sets), **E** PROGENy pathway activity, **F** KAT8 MSL/NSL complex output with mito-encoded
  and ribosomal specificity controls, **G** cross-layer TF convergence. Writes a
  `SANITY_part5_*.csv` ledger of every check. Uses GO annotations from `org.Mm.eg.db`
  (offline, non-circular) — no matrisome download, no hand-curated seed lists.
- `config_downstream.R` — shared config; **edit the input paths here**. Also holds the TF
  highlight sets (`tf_adipogenic`, `tf_inflammatory`) and the consensus toggle for part4.

**Prediction part4 (TF) tests (from the part4b result):** the metabolic-identity collapse is
non-cell-autonomous, so the adipogenic TFs (Pparg, Cebpa, Srebf1, …) should read strongly
repressed in **iWAT** but **flat in the cells**. Whatever *is* cell-autonomous (ECM/fibrosis
program) should surface as the TFs active in the cells column.

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
source("downstream_analysis/part4_tf_activity.R")            # TF activity
source("downstream_analysis/part4b_cell_tissue_concordance.R") # concordance
```

Outputs land in `downstream_analysis/results/`. `part4_tf_activity.R` needs internet **once** to
fetch the CollecTRI network from OmniPath; it is then cached to
`downstream_analysis/collectri_mouse.rds`.

New packages used: `decoupleR`, `OmnipathR`, `fgsea` (Bioconductor); `pheatmap`,
`ggrepel`, optional `qvalue`. The scripts install anything missing.

## How to read the results

**Organizing hypothesis (from your slides):** KAT8 is an *activator*, so its *direct*
effect is loss of expression. Prediction — the **down** metabolic/adipocyte-identity
program is **direct & cell-autonomous** (reproduces in cells); the **up** inflammatory
program is **indirect & non-cell-autonomous** (infiltration, tissue-only).

- **part4 heatmap** — a TF active in *cells and* tissue is a cell-autonomous regulator;
  a TF active only in tissue (NF-κB/IRF/STAT) is likely infiltration-driven. Watch the
  adipogenic TFs (Pparg, Cebpa, Srebf1, nuclear receptors) in the **cells** column —
  their repression is your cell-autonomous mechanism.
- **part4b scatter + GSEA** — a positive, significant `NES(tissue_DOWN)` in cells means
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
   With few cell DEGs a Venn is underpowered and misleading. `part4b` instead uses rank
   correlation, a signature-GSEA of the tissue program against the full cell ranking, and
   directional binomial tests that borrow power from the tissue side — all threshold-free on cells.

5. **[L] PCA on VST.** *Recommended.* The cell script runs `prcomp(t(vst), scale.=TRUE)` on all
   genes; unit-scaling every gene upweights noisy low-variance genes. DESeq2 convention is the
   top ~500 most-variable genes without unit scaling (`DESeq2::plotPCA`). Cosmetic, not a
   conclusion-changer.

6. **[L] Cross-dataset gene matching.** *Applied.* `part4b` merges cells↔tissue on gene symbol and
   de-duplicates by keeping the most significant row, so `make.unique` suffixes don't silently
   drop genes.

7. **[L] Permutation calibration.** *Applied.* `part4_tf_activity.R` can permute gene labels (`n_permutations`)
   to confirm the analytic TF p-values aren't inflated at this sample size.

The validated tissue design (`~ 0 + GroupDepot`, sex handling) is sound and left unchanged.

## What still needs new data (out of scope here, but worth stating in the paper)
- **snRNA-seq** — the definitive test of *which cell type* drives each signature.
- **H4K16ac / KAT8 ChIP-seq** — the definitive *direct-target* map. Until then, integrating
  *public* MOF/H4K16ac ChIP (Cistrome/ENCODE) via BETA is the strongest no-new-experiment proxy
  (a possible follow-up analysis).
- KAT8 also has **non-histone substrates**, so some effects are post-translational and invisible
  to any RNA-based readout — state as a limitation.


---

## Pipeline hygiene (added during the audit)

**Driver.** `run_all.R` at the project root runs everything in dependency order and
**stops at the first failure**, so a mid-chain error can never leave a half-populated
run folder for the next script to read:

```bash
Rscript run_all.R tissue      # parts 1-4
Rscript run_all.R cells       # cells DE/ORA/GSEA
Rscript run_all.R downstream  # TF activity, concordance, mechanism
Rscript run_all.R all         # everything
```
Each part runs in its own environment so variables cannot leak between scripts.

**Run-folder resolution.** Part 1 now writes `savepoints/LATEST_RUN.txt`. Parts 2/3/4 read
that pointer, falling back to a folder-**name** sort. They no longer sort by mtime —
OneDrive sync updates mtimes and could silently promote a stale run folder.

**Sample-annotation validation.** Part 1 validates the positional Depot/Sex/Genotype
assignment against the data itself: `Xist` vs Y-linked genes for sex, `Kat8` down in every
KD group for genotype, plus a design-balance table. Set `STOP_ON_ANNOTATION_FAIL <- TRUE`
to halt instead of warn.

**Cell DEG definition.** `parameters.R` now owns the cell thresholds (`cells_params`).
Both DEG definitions are computed and saved (`DEG_lists_cells_*.rds`); `primary_deg_rule`
selects which is authoritative. Default is `"fdr_only"` — the `|log2FC| > 0.5` gate
discarded 290 of 392 real cell DEGs (74%), because median |log2FC| in these cells is ~0.37.

**Bugs fixed in this pass**
| Where | Bug | Effect |
|---|---|---|
| `part3_fgsea.R` | `gsea_go <- gsea_go_sig` inside `error=function(e)` | assigned to the handler's local scope, so a failed `simplify()` silently saved the **full unfiltered** GO result |
| `part2/3/4` | run folder chosen by mtime | OneDrive sync could promote a stale run; wrong tables read *and* written |
| `part3_fgsea.R` | `gseKEGG(pvalueCutoff = 0.1)` vs `gseGO(1.0)` | `N_Pathways_Tested` counted only pathways already passing p<0.1 — not comparable to the GO row |
| cells script | `use_ora_universe <- FALSE` | liberal background inflated ORA significance and disagreed with `part2_ora.R` (TRUE) |
| `part5` (new) | `FDR` constant shadowed by the `FDR` column | `ifelse(FDR < FDR, …)` always FALSE — caught and fixed before release |
| `part5` (new) | leading-edge stat vector unnamed | reproduced the earlier silent-empty-table bug; now uses a named vector and asserts non-zero rows |
