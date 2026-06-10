# PROMPT: Compare our KAT8-KD adipose data to the Mof+/- het paper (GSE156463)

> Copy everything below the line into a new Claude Code session opened in the
> KAT8-Project repo. It is self-contained — you should not need to refer back
> to any prior chat.

---

## ROLE / GOAL
You are helping me (a metabolism/epigenetics researcher) compare my own bulk
RNA-seq dataset against a published dataset, to quantify where the two KAT8/MOF
loss-of-function models AGREE and DIFFER, at both the RNA-seq and ChIP-seq
level. Be rigorous, show effect sizes and CIs (not just p-values), and flag
confounds honestly. Produce runnable R code, run it where possible, save
figures/tables, and write a short interpretation memo at the end. Commit and
push results to my working branch.

## MY DATASET ("OURS")
- Model: **adipocyte-restricted KAT8 KNOCKDOWN** (Genotype levels: CTL vs
  KAT8KD), C57BL/6, BOTH sexes, n~4-5 per Sex x Genotype group.
- Two depots: **iWAT** (subcutaneous) and **gWAT** (visceral/gonadal).
- Data + pipeline already in this repo:
  - `counts.txt` — raw gene-count matrix; has a `gene_name` column, samples as
    columns. (If not present in the working dir, ask me for its path; it lives
    on my OneDrive analysis machine.)
  - `parameters.R` — thresholds (FDR cut, lfc cut, prefilter settings).
  - `part7_sex_analysis_consolidated.R` — THE reference script. It loads
    `counts.txt`, builds `sample_annot` (columns: Sample, Sex [F/M],
    Genotype [CTL/KAT8KD], depot), filters to `count_matrix`, runs per-sex
    DESeq2 (`~ Genotype`) giving `de_m` and `de_f` (data.frames with
    `log2FoldChange`, `padj`, `gene`), and the interaction model
    `~ Sex + Genotype + Sex:Genotype` (LRT vs `~ Sex + Genotype`).
  - Reuse that script's data-loading and DESeq2 conventions so results are
    comparable. Do NOT re-invent the loading logic.
- Prior findings I have already established (for context, don't re-derive):
  - M-vs-F concordance is ~99% (same direction), Pearson r ~0.6-0.73.
  - Genotype explains ~17-33% median per-gene variance; Sex ~7-11%;
    Sex x Genotype only ~2% (real but minority signal).
  - The Sex x Genotype interaction genes (630 iWAT / 374 gWAT) concentrate on
    the PPAR / fatty-acid / lipid-storage / insulin-response axis.

## THE PUBLISHED DATASET ("THEIRS") — GSE156463
- Paper: Pessoa Rodrigues C, Chatterjee A, et al., Akhtar A. "Histone H4 lysine
  16 acetylation controls central carbon metabolism and diet-induced obesity in
  mice." Nat Commun 2021;12:6212. doi:10.1038/s41467-021-26277-w.
- Model: **whole-body Mof (Kat8) HETEROZYGOTE (Mof+/-)** — ~50% global,
  constitutive haploinsufficiency. NOT adipocyte-restricted, NOT a knockdown.
- GEO **GSE156463** contains BOTH:
  1. **Bulk WAT mRNA-seq**: visceral WAT (≈ my gWAT), **MALE ONLY**, Mof+/- vs
     Mof+/+, on standard diet (SD; triplicates) and 26-wk high-fat diet (HFD;
     duplicates). STAR -> GRCm38 -> featureCounts -> DESeq2 in the original.
  2. **ChIP-seq** in WAT (male, SD): **MOF** and **H4K16ac** (normalized as
     log2 over H3 IP; peaks called with MACS2; mm10/GRCm38).
- Their central mechanism (my anchor): MOF binds the **Pparg** promoter ->
  deposits H4K16ac -> drives **Pparg** -> **Pgc1a (Ppargc1a) / Mef2c** ->
  **Glut4 (Slc2a4)** -> glucose uptake -> neutral lipid storage. Loss of MOF
  collapses this -> obesity-resistant but glucose-intolerant. Rbp4 goes UP
  (insulin-resistance marker). Rescue by TZD (PPARg agonist), ectopic Glut4,
  or SIRT1 inhibition. NOTE: in their ChIP, MOF bound Pparg promoter (positive
  control) but NOT Glut4/Pgc1a/Mef2c promoters (negatives).

## CRITICAL CAVEATS (state these in the memo)
- Their RNA-seq is **MALE-ONLY and underpowered (n=2-3)** — I cannot get a
  female comparison from their data. Use their MALE Mof+/- signature only to
  (a) validate that my KD reproduces the canonical KAT8 program and (b) anchor
  the PPARg-Glut4 module. The female/sex story stays in MY data.
- **Het (50% global, constitutive)** vs **adipocyte KD** are different doses
  and compartments. Differences are expected and interesting — frame as
  dose/compartment, NOT as contradiction.
- Match conditions: compare THEIR **SD** Mof+/- vs +/+ to MY **baseline KD**
  (no extra diet stress). Use my **gWAT MALE** arm (`de_m` for gWAT) as the
  primary comparator since their tissue is visceral and male.
- Their "no sex bias" claim is **weight gain only**; it does NOT contradict a
  sex-resolved transcriptome (nobody did that before — that's my novelty).

## TASKS

### TASK 1 — Acquire and process THEIRS (GSE156463)
- Download via `GEOquery::getGEO("GSE156463")` + `getGEOSuppFiles("GSE156463")`.
  Prefer the processed count matrix / DEG table from the supplementary files;
  if only raw counts are available, run DESeq2 (`~ Genotype`, Mof+/- vs +/+) on
  the SD samples. Also grab Supplementary Data 3 from the paper if present.
- Map their gene IDs to symbols (GRCm38) so they merge with my `gene_name`s.
- Report their n per group, what files you found, and exactly how you derived
  their logFC/padj. If download fails (no internet), tell me precisely which
  files to download manually and where to put them, then proceed.

### TASK 2 — RNA-seq CONCORDANCE (theirs male WAT vs my male gWAT KD)
Primary comparison: their SD Mof+/- vs +/+ logFC  ↔  my `de_m` (gWAT) logFC.
Also run my male iWAT as a secondary. Compute and visualize:
  a. **Pearson AND Spearman** of logFC on shared genes (all genes; and
     restricted to genes significant in either dataset). Scatter with the
     identity line, colored by concordant/discordant quadrant.
  b. **RRHO2** (rank-rank hypergeometric overlap) — threshold-free
     directional concordance; show the heatmap; report the strongest
     concordance quadrants (expect down-down and up-up).
  c. **GSEA / fgsea**: build "their UP" and "their DOWN" SD DEG sets, test for
     enrichment in MY ranked male gWAT logFC list. Report NES + padj.
  d. **Overlap stats**: Fisher's exact / odds ratio for shared
     directional DEGs; a 2x2 (up/down) contingency and a signed Venn.
  e. **Where do they DIFFER**: list genes strongly DE in one model but not the
     other (and opposite-direction genes). ORA (GO:BP + KEGG, clusterProfiler,
     org.Mm.eg.db) on the "ours-only", "theirs-only", and "discordant" sets to
     see which pathways are model-specific.

### TASK 3 — The PPARg-GLUT4 axis, head to head
- Module gene set (their established effector axis):
  `Pparg, Ppargc1a, Mef2c, Slc2a4, Cd36, Fabp4, Lpl, Plin1, Plin4, Cidec,
   Adipoq, Scd1, Lipe, Pnpla2, Rbp4`.
- Tabulate logFC + padj for every module gene in: theirs (SD), my male gWAT,
  my female gWAT, my male iWAT, my female iWAT. One heatmap, rows = genes,
  columns = those 5 contrasts.
- Test whether the module moves the SAME direction in theirs vs mine, and
  whether MY data shows a Sex x Genotype interaction on the module score
  (per-sample mean z-score; `lm(score ~ Sex*Genotype)`), per depot.
- Specifically check **Rbp4** (their insulin-resistance marker): does it go up,
  and more in males?

### TASK 4 — ChIP-seq (theirs) integrated with my DE
- From GSE156463, get the **H4K16ac** and **MOF** ChIP signal/peaks in WAT.
  Use bigWig + MACS2 peaks (mm10/GRCm38). If processed peaks aren't in GEO,
  tell me what to download and proceed with what's available.
- Build a **target/control panel**: Pparg main-isoform promoter (their Fig 6h
  positive control), Slc2a4/Glut4, Ppargc1a, Mef2c (their NEGATIVE controls —
  MOF did not bind), plus MY interaction-gene promoters (Elovl3, Scd1, Cd36,
  and the top female-amplified iWAT genes).
- Key integrative question: **are MY KD-downregulated genes enriched for
  MOF/H4K16ac binding in their WAT ChIP?** Test: for genes bound by
  MOF/H4K16ac (TSS +/- promoter window) vs not, is the KD logFC more negative?
  (boxplot + Wilcoxon; and a binding-vs-|logFC| trend). Do this for my male
  and female arms separately — if female KD genes are LESS dependent on
  MOF/H4K16ac binding, that hints at the buffering mechanism.
- Output a **pre-registration table** of loci + predicted direction for my own
  planned ChIP/CUT&RUN, anchored to their coordinates.

### TASK 5 — Memo + deliverables
- Write `GSE156463_comparison_memo.md`: (1) how concordant our KD is with the
  Mof+/- model (one headline number + RRHO + GSEA), (2) the biggest
  model-specific differences and likely reasons (dose/compartment/diet/sex),
  (3) what the ChIP integration says about which of my DE genes are direct
  MOF/H4K16ac targets, (4) a bullet list of testable hypotheses with
  predictions, (5) explicit limitations.
- Save all figures (PNG, 150 dpi) and tables (CSV) under a clearly named
  results subfolder with a run timestamp. List every output file.
- Reuse my existing helper conventions (DESeq2 loading, `run_ora_block`-style
  ORA, color schemes) from `part7_sex_analysis_consolidated.R` where useful.

## GIT
- Develop on branch `claude/debug-concordance-plots-MY9Dc` (create if needed).
- Commit logically; push with `git push -u origin claude/debug-concordance-plots-MY9Dc`.
- Do not commit large downloaded raw data — add it to `.gitignore`; commit only
  scripts, figures, tables, and the memo.

## STYLE
- Show me the effect sizes and CIs, not just significance.
- When their data is male-only/underpowered, say so and don't over-interpret.
- If something can't be computed (missing file, no internet), tell me exactly
  what to provide rather than silently skipping.
