# KAT8 downstream analysis

Four R scripts, **one question each**, all run from the three DESeq2 result tables
your main pipeline already writes. No new data, no Python, no counts matrix.

| Script | Question it answers |
|---|---|
| `part4b_cell_tissue_concordance.R` | Is the KAT8 effect **cell-autonomous**? (gene-level + program-level) |
| `part5a_secretome.R` | What does the KAT8-null adipocyte **broadcast** to its neighbours? |
| `part5b_tf_analysis.R` | Which **transcription factor** drives it? (3 independent layers) |
| `part5c_kat8_complex.R` | Is the canonical **KAT8/NSL** program engaged? |

`config_downstream.R` holds the shared settings (input paths, thresholds, colours,
plot theme, sanity-check helpers). **Edit the input paths there, nowhere else.**

## Setup

Put the three DE tables in `downstream_analysis/data/` (or point `DE_PATHS` at them):

| Contrast | Rename your pipeline's output to |
|---|---|
| Cells | `data/DE_cells_KAT8KD_vs_CTL.csv` |
| iWAT | `data/DE_tissue_iWAT_KD_vs_CTL.csv` |
| gWAT | `data/DE_tissue_gWAT_KD_vs_CTL.csv` |

Each needs `gene_name`, `log2FoldChange`, `stat`, `padj` — your current outputs already have them.

## Running

```r
# from the project root
source("downstream_analysis/part4b_cell_tissue_concordance.R")
source("downstream_analysis/part5a_secretome.R")
source("downstream_analysis/part5b_tf_analysis.R")
source("downstream_analysis/part5c_kat8_complex.R")
```
or let the driver do it:

```bash
Rscript run_all.R              # everything: tissue -> cells -> downstream
Rscript run_all.R downstream   # just these four
```

The driver **automatically copies** the newest `DE_tissue_*` / `DE_cells_*` tables out of
`savepoints/RUN_*/tables/` into `downstream_analysis/data/` before the downstream steps run,
so you don't have to move files by hand. (The tissue and cells scripts each create their own
timestamped run folder, which is why this staging step exists.) If a table can't be found the
downstream steps are **skipped with a message** rather than failing halfway.

Everything lands in `downstream_analysis/results/`, and each script writes a
`SANITY_*.csv` listing every check it performed.

**Offline by design.** The TF regulon is bundled (`regulon_dorothea_mm_ABC.csv`,
mouse DoRothEA confidence A–C) so nothing depends on the OmniPath server, which
has an intermittent bug that breaks `get_collectri()`. Set `tf_params$network_csv <- NULL`
in the config if you'd rather fetch CollecTRI live.

Packages: `decoupleR`, `fgsea`, `org.Mm.eg.db`, `GO.db`, `AnnotationDbi` (Bioconductor);
`msigdbr`, `pheatmap`, `ggrepel` (CRAN). Optional: `progeny`, `qvalue`. Scripts install what's missing.

## How to read the results, in order

**1. Is it cell-autonomous?** (`part4b`) — the program-level figure is the headline.
A program that is strong in tissue but **flat in the cells** is *not* cell-autonomous.
A program that reproduces in the cells is the direct-effect candidate.

**2. What does the cell broadcast?** (`part5a`) — secreted/surface genes changed in the
cells that move the **same direction in vivo**. These are the candidate signals linking an
apparently healthy KAT8-null adipocyte to a diseased tissue.

**3. Which TF?** (`part5b`) — read the **convergence table first**. Three deliberately
independent layers:
- *Layer 1* — TF genes that themselves change (no regulon, no inference)
- *Layer 2* — regulon activity (decoupleR); a re-weighted score on the *same* statistics,
  so a summary rather than independent evidence
- *Layer 3* — MSigDB C3:TFT, ChIP-derived target sets, independent of layer 2

A TF supported by **2–3 layers** is a real candidate. A layer-2-only hit is not.
`part5b_L2_leading_edge_*.csv` lets you audit any call gene by gene.

**4. Canonical KAT8 program?** (`part5c`) — only claim an NSL-specific effect if the
printed verdict says `SPECIFIC`, i.e. nuclear-encoded OXPHOS moved while **both** controls
(mito-encoded and ribosomal genes) stayed flat. The verdict is computed from the data,
never hardcoded.

## What still needs new data
- **snRNA-seq** — the definitive test of which cell type drives each signature.
- **H4K16ac / KAT8 ChIP-seq** — the definitive direct-target map. Everything here is
  inference; it tells you what to ChIP for.
- KAT8 also has **non-histone substrates**, so some effects are post-translational and
  invisible to any RNA-based readout. Worth stating as a limitation.
