# KAT8-KD RNA-seq: what each analysis does and why

A working explanation of the pipeline — written to be understood and re-explained,
not just cited. Every number and setting below comes from the scripts themselves
(`parameters.R` and the part scripts), not from memory.

Each section follows the same shape:

- **The question** — what we are actually asking
- **How it works** — the mechanics, plainly
- **Why this way** — what the obvious alternative would have got wrong
- **Reading the output** — what a good result looks like, and how to be wrong

A journal-ready Methods paragraph is at the end.

---

## 0. The experiment

KAT8 knockdown by siRNA in mouse adipose tissue, plus a parallel knockdown in
3T3-L1 adipocytes.

- **Two depots:** iWAT (inguinal, subcutaneous) and gWAT (gonadal, visceral)
- **Design:** 2 × 2 per depot — Sex (F, M) × Genotype (CTL, KAT8-KD)
- **Libraries:** 40 sequenced, 2 removed as outliers (JS_08, JS_28), **38 analysed**
- **Per depot:** 19 samples — CTL-F 5, KD-F 4, CTL-M 5, KD-M 5

KAT8 is the catalytic subunit that writes **H4K16 acetylation**. It is not a
transcription factor; it is a co-activator. That single fact drives most of the
analysis design: you cannot find KAT8's "targets" by looking for a motif, and
because it *activates*, its direct consequence should be genes going **down**.

---

## 1. Differential expression (Part 1)

### The question
Which genes change when KAT8 is knocked down, in each depot?

### How it works

**Prefiltering.** Genes are kept if they have **≥ 10 counts in ≥ 4 samples**,
applied across the whole experiment rather than per group. Genes below that are
noise — they cost statistical power (every gene tested makes the multiple-testing
correction harsher) and return nothing.

**Normalisation.** DESeq2's median-of-ratios. Libraries differ in sequencing
depth; you cannot compare raw counts. Rather than dividing by total reads — which
a handful of enormously expressed genes would distort — DESeq2 computes, for each
sample, the median ratio of each gene to that gene's average across all samples.
That median is the size factor. It is robust because a few extreme genes cannot
move a median.

**The model.** `~ 0 + GroupDepot`, fitted once across all 38 samples.

`GroupDepot` is a single factor combining depot and genotype, so its levels are
`iWAT_CTL`, `iWAT_KAT8KD`, `gWAT_CTL`, `gWAT_KAT8KD`. The `0 +` means "no
intercept" — instead of one baseline plus offsets, each group gets its own mean
estimated directly. This is a **cell-means** parameterisation. Depot-specific
KAT8 effects are then pulled out as contrasts: `iWAT_KAT8KD − iWAT_CTL`, and
likewise for gWAT.

**Dispersion.** RNA-seq counts are over-dispersed — the variance exceeds the mean,
because biological replicates genuinely differ. DESeq2 estimates a dispersion per
gene, then shrinks each gene's estimate toward a fitted trend across all genes.
With n ≈ 5 per group, a per-gene variance estimate on its own is unreliable;
borrowing strength across thousands of genes is what makes the test work at this
sample size.

**Testing.** Wald test per gene, Benjamini–Hochberg FDR correction.

**Calling a DEG.** `FDR < 0.05` **AND** `|log2FC| > 1` (i.e. a two-fold change),
for both depots.

**Shrunken fold changes.** `lfcShrink` with `ashr` is also computed and exported
as an extra column. Genes with low counts produce wild fold-change estimates —
2 counts vs 8 counts is "4-fold" but means nothing. Shrinkage pulls poorly
supported estimates toward zero while leaving well-supported ones alone. It is
exported for ranking and plotting but **DEG calling and GSEA ranking are
unchanged**, so the validated counts are unaffected.

**VST.** A variance-stabilising transform produces a matrix on a log-like scale
where variance no longer depends on the mean. Used for PCA, clustering,
correlation and heatmaps — never for the statistical testing, which always runs
on raw counts.

### Why this way
The alternative is fitting each depot separately. The unified model was chosen
because dispersion is estimated once across all 38 samples rather than twice
across 19, which makes the variance estimates more stable. The trade-off is real
and shows up later: a gene contaminated in one depot gets an inflated dispersion
that propagates into the *other* depot's contrast (see §7).

### Reading the output
`DE_tissue_<contrast>_<tag>.csv`, one row per gene: `baseMean` (average
normalised expression), `logFC`, `lfcSE` (the standard error — the uncertainty on
that fold change), `stat` (Wald statistic), `pvalue`, `padj`.

**Result:** iWAT 3,022 DEGs, gWAT 1,159. The response is strongly depot-asymmetric.

---

## 2. Should males and females be pooled? (the sex analysis)

### The question
Does KAT8 knockdown do something *different* in males than in females? If yes,
the sexes must be analysed separately — which halves the sample size. If no, sex
can be included as a covariate, or omitted entirely if the design is balanced.

### How it works — five independent lines of evidence

**(a) Likelihood ratio test for interaction.** Run per depot. Two models are
fitted to every gene:

- Full: `~ Sex + Genotype + Sex:Genotype`
- Reduced: `~ Sex + Genotype`

The LRT asks whether adding the interaction term improves the fit enough to
justify the extra parameter. `Sex:Genotype` captures genes where the KD effect
*differs* by sex — up in female KD but down in male KD would be a strong
interaction; up equally in both is none.

Result: iWAT 630/16,466 significant (3.83%); gWAT 374/16,827 (2.22%). At
FDR < 0.05 you expect ~5% false positives anyway, so this is a low rate.

**(b) π₀ estimation (Storey & Tibshirani).** If no gene truly had an interaction,
the p-values would be uniformly distributed. Real signal produces a spike near
zero on top of that flat background. π₀ measures how much of the distribution is
flat, i.e. what proportion of genes are true nulls.

Result: iWAT π₀ = 0.87, gWAT 0.85 — 85–87% of genes genuinely have no
sex-dependent KD response.

**(c) Variance decomposition.** ANOVA per gene on VST values, asking what fraction
of each gene's variance is attributable to Sex, Genotype, and their interaction.

| Term | iWAT (median) | gWAT (median) |
|---|---|---|
| Sex | 6.5% | 11.5% |
| Genotype | 33.2% | 16.7% |
| Sex × Genotype | 2.1% | 2.4% |
| Residual | 46.3% | 50.1% |

The interaction gets ~2% against genotype's 17–33% — eight to sixteen times less.

**(d) Sex-stratified concordance.** DESeq2 run on females only, then males only,
and the two sets of fold changes correlated.

Result: iWAT r = 0.727, gWAT r = 0.605; **directional** concordance 99.0% and
99.5%. The correlation is not higher because per-sex analyses at n = 4–5 are
noisy, but essentially every significant gene moves the *same way* in both sexes.
The lower gWAT r reflects magnitude differences, and visceral fat is known to be
more sexually dimorphic at baseline than subcutaneous.

**(e) Model comparison.** Fit with and without Sex, count DEGs.

- Per depot (n = 19): adding Sex **gains** DEGs (+14.1% iWAT, +5.0% gWAT)
- Unified (n = 38): adding Sex **loses** DEGs (−8.1% iWAT, −6.0% gWAT)

The reversal is a degrees-of-freedom argument. At n = 19, spending one df on Sex
costs 5.6% of your residual evidence, but Sex absorbs enough variance to more than
repay it. At n = 38 the same df costs only 2.9%, *and* the balanced design already
prevents sex from confounding genotype — so the cost now outweighs the benefit.

### The decision
`~ 0 + GroupDepot`, unified, no Sex covariate. All five lines agree, logFC
correlation between the ± Sex models is r = 0.94–0.97, so the biology is
unchanged either way.

### How to explain it in one line
> Sex changes the *baseline* in fat, but it does not change *how fat responds to
> losing KAT8* — and because the design is balanced, sex cannot masquerade as a
> genotype effect.

---

## 3. Pathway analysis: ORA (Part 2) and GSEA (Part 3)

These answer different questions and must not be conflated.

### ORA — over-representation analysis

**The question.** Among the genes that passed my threshold, are any pathways
present more often than chance?

**How it works.** Take the DEG list (FDR < 0.05 and |log2FC| > 1), split by
direction, and for each gene set ask: of my 800 up-genes, 40 are in "fatty acid
oxidation" — given that pathway's size and the size of the background, how
surprising is 40? A hypergeometric test answers it. Think of drawing balls from
an urn without replacement.

**The background matters more than people realise.** The universe here is **all
genes tested by DESeq2**, not all genes in the genome. If you used the whole
genome as background, every result would be inflated by the fact that you only
sequenced expressed genes in the first place — adipose pathways would look
enriched simply because you were looking at adipose.

**Settings.** `clusterProfiler` against GO:BP and KEGG; gene sets restricted to
5–500 genes (smaller are unstable, larger are too vague to interpret);
`pvalueCutoff = 0.1`, `qvalueCutoff = 0.2` at the enrichment step; GO:BP
redundancy reduced with `simplify(cutoff = 0.7)`.

**Why simplify.** GO is a hierarchy, so "cellular respiration", "aerobic
respiration" and "oxidative phosphorylation" all fire together and are largely
the same finding. `simplify` measures semantic similarity between terms and keeps
one representative per cluster. Without it a results table is twenty restatements
of three findings.

### GSEA — gene set enrichment analysis

**The question.** Ignoring thresholds entirely: is any pathway systematically
shifted toward one end of my ranked gene list?

**How it works.** Rank **every** tested gene by the **DESeq2 Wald statistic** —
not log2FC. Walk down the ranked list keeping a running score that rises on
pathway members and falls on non-members. The maximum deviation from zero is the
enrichment score; significance comes from permutation.

**Why rank on the Wald statistic.** The Wald stat is `logFC / lfcSE` — effect size
divided by its uncertainty. Ranking on logFC alone puts noisy low-count genes with
huge fold changes at the top. The Wald stat demotes them automatically, because a
large fold change with a large standard error scores low.

**Why GSEA as well as ORA.** ORA throws away everything below threshold. If a
pathway's fifty genes all move 1.5-fold in the same direction, ORA sees nothing —
none passed |log2FC| > 1 — while GSEA sees a coordinated shift. Coordinated
subthreshold movement is exactly what a chromatin modifier is expected to produce.

**Settings.** `gseGO` / `gseKEGG`, gene sets 5–500, FDR < 0.05, seed 12345,
GO:BP simplification at 0.5 applied to the top 300 terms by p-value (the
similarity computation is O(n²) and hangs on larger inputs).

### Reading them together
Agreement between ORA and GSEA is strong evidence. ORA-only suggests a few large
effects; GSEA-only suggests a broad coordinated shift too subtle to pass a
threshold. Neither is wrong — they are different-shaped findings.

---

## 4. Is the effect cell-autonomous? (Part 4b)

### The question
The tissue result could mean two completely different things:

1. KAT8 loss changed what adipocytes express — **cell-autonomous**
2. KAT8 loss changed the *cell composition* of the tissue — macrophages moved in,
   fibrosis developed — and we are reading the new cell mixture

A bulk RNA-seq tissue sample cannot tell these apart on its own. Both look like
"gene expression changed".

### How it works
3T3-L1 adipocytes are a pure adipocyte system: no immune cells, no vasculature,
no infiltrate. Overlaying the cell result on the tissue result splits every change
into three classes:

| Pattern | Interpretation |
|---|---|
| Changed in **both**, same direction | **Cell-autonomous** — adipocyte-intrinsic |
| **Tissue only** (especially inflammatory, up) | **Non-cell-autonomous** — composition or signalling |
| **Cell only** | In-vitro-specific |

Because KAT8 is an **activator**, its direct consequence should be
down-regulation — so the **concordant-down** set is the prime direct-target
candidate list.

**The cell side is threshold-free.** Rank correlation, gene-set enrichment and
directional tests that borrow power from the well-powered tissue DEGs. This
matters: if the comparison required the cells to produce their own long DEG list,
a modest cell experiment would limit everything. Asking instead "do the tissue's
DEGs move consistently in the cells?" is a far more powerful question.

### The critical interpretive point
In KAT8-KD 3T3-L1, the adipocyte-marker, lipid and thermogenic programs are
**flat, with tight confidence intervals** — even though the knockdown is
*stronger* than in tissue.

This is **not** a failed experiment, and reading it that way is the single easiest
mistake to make with this dataset. The 3T3-L1 work is a **maintenance** knockdown:
KAT8 is removed *after* the cells have differentiated. Published work shows KAT8
is required when lost **before** differentiation. So a flat result is a real,
informative finding — *the mature adipocyte program does not need ongoing KAT8*.
Tight confidence intervals around zero mean "genuinely unchanged", not
"underpowered".

That result is what makes Part 5a necessary.

---

## 5. What does the KAT8-null adipocyte broadcast? (Part 5a)

### The question
Part 4b left a puzzle: if the KAT8-null adipocyte still looks and behaves like a
healthy adipocyte, why is the tissue inflamed and fibrotic?

The obvious candidate is what the cell **sends** to its neighbours.

### How it works
Take the genes that change in the cells, restrict to those annotated as
**secreted or cell-surface**, then test whether each one reproduces *in vivo* —
does it move the same direction in a depot?

A gene changed in cells **and** moving the same way in tissue is a candidate
non-cell-autonomous signal: the mechanistic bridge between an apparently healthy
KAT8-null adipocyte and a diseased tissue.

### Why this is the right follow-up
It converts a negative result into a directed hypothesis. "The cells look fine"
is only a dead end if you stop at the cell's own phenotype; the secretome asks what
the cell does *to things around it*, which is exactly the gap between a normal
cell and an abnormal tissue.

---

## 6. Which transcription factor drives it? (Part 5b)

### The question
KAT8 is a co-activator, not a TF. It has no motif to search for. So: which TF
*programs* shift when it is lost?

### How it works — three deliberately independent layers

**Layer 1 — TF genes that themselves change.** Just look at whether TF-encoding
genes are differentially expressed. Zero assumptions, zero inference. Also the
weakest, because a TF's *activity* is often regulated without its mRNA changing at
all.

**Layer 2 — regulon activity (decoupleR, univariate linear model).** A regulon is
a curated set of a TF's known target genes. If most of PPARγ's targets went down,
PPARγ activity probably went down, whether or not `Pparg` mRNA moved. ULM
regresses the DE statistics on regulon membership to score each TF.

**Layer 3 — MSigDB C3:TFT target enrichment.** Gene sets built from **ChIP-derived
binding data** — an entirely different evidence base from Layer 2's curated
network.

### Why three layers, and why Layer 2 must not be read alone
Layer 2 is the one people quote, and on its own it is the weakest of the three:

- It is a **re-weighted score computed from the same statistics as the DE result**,
  so in part it restates what you already know
- Its curated regulon is **largely non-adipose** and thin for adipocyte TFs
- It is **blind to co-activators** — the very class KAT8 belongs to

**Agreement across layers is the evidence. No single layer is.** Section 6 of the
script tabulates that convergence explicitly, and that table is what should be
reported.

**Also included:** PROGENy, which scores signalling pathway activity one level
*above* TFs; and a cells-vs-tissue comparison, so a TF active in **both** can be
distinguished from one that is tissue-only and therefore likely driven by
infiltrating cells rather than by adipocytes.

---

## 7. Is the canonical KAT8 program engaged? (Part 5c)

### The question
KAT8 is the catalytic subunit of **two** complexes:

- **MSL** (Msl1/2/3) — broad H4K16ac, general transcriptional output
- **NSL/KANSL** (Kansl1–3, Mcrs1, Phf20, Ogt, Wdr5, Hcfc1) — promoter-proximal,
  and the best-documented direct regulator of **nuclear-encoded OXPHOS genes**

Did knocking down KAT8 actually engage its canonical program?

### The problem this analysis solves
You cannot answer it by looking at whether the complex members changed. **siRNA
depletes KAT8 protein, not its partners' transcripts** — so partner mRNA stays
flat *by design*. The complex effect is invisible in membership expression. It has
to be read from the complex's regulatory **output**: its target genes.

### How it works — and why the controls are the whole point

| Gene set | Expected if NSL is genuinely hit |
|---|---|
| Nuclear-encoded OXPHOS (NSL targets) | **moves** |
| Mitochondrially encoded (`mt-`) genes | **flat** — transcribed by mtDNA, outside NSL control |
| Ribosomal proteins | **flat** — comparable housekeeping, generic-drift control |

Only if nuclear OXPHOS moves **while both controls stay flat** is the shift
specifically NSL, rather than a normalisation artefact or a generic stress
response. If a control also moves, complex specificity cannot be claimed.

**Thresholds.** `SPECIFIC` requires a **≥ 2.0×** margin over every control;
**1.5–2.0×** is reported as `BORDERLINE`. The verdict is computed from the data,
never written into a caption.

### The result, and how to report it honestly
Using the curated 50-gene nuclear-OXPHOS list, iWAT is **SPECIFIC at 2.54×**.
Using GO:0006119 (oxidative phosphorylation) minus mito-encoded genes, the same
data give **BORDERLINE at 1.96×**.

That is **threshold-adjacency between two defensible gene sets, not a
contradiction.** Quoting 2.54× as though it were robust would be overstating it.
The honest report is borderline-to-specific, naming which list produced which
number. The sensitivity check exists precisely so this cannot be quietly ignored.

---

## 8. Quality control (Part 1b)

Six checks. What each one is actually for:

**A. p-value histograms.** The single best model diagnostic. A healthy contrast
gives a flat histogram with a spike near zero. A *hill in the middle* means the
test is mis-calibrated or the model is misspecified. A spike at 1 means filtering
is over-conservative. Anything but the first pattern makes the results suspect no
matter how good the biology looks.

**B. Independent filtering.** DESeq2 sets `padj = NA` for genes it filters out to
improve power. If a large fraction is NA, the filter is eating your data.

**C. Unstable log2FC estimates.** Flags significant genes with `lfcSE > 2` — a
fold change nobody can stand behind, usually driven by a few libraries.

**D. Contamination screen.** Tissue-restricted marker panels (sperm/testis, blood,
muscle, liver, lymph node), scored per sample. Three questions per panel: is any
library an outlier; does the score track genotype *within depot*; and — the test
that actually decides it — **does the association survive dropping its strongest
single gene?** Carry-over transfers a whole tissue, so it survives losing any one
marker. Regulation of one gene does not. Survival only counts if something
tissue-**restricted** still carries it, since shared markers are expressed in fat
in their own right.

**E. PCA and covariate association.** Which experimental variables explain which
principal components.

**F. Sample–sample correlation.** Libraries that correlate poorly with everything
else — with their annotation printed, because a whole group differing from the
rest is a treatment effect, not a bad library.

### What this dataset's QC found — all resolved

**Alb is down in both depots — regulation, not liver contamination.** No outlier
library; every other hepatocyte-restricted marker sits at background (VST ~6);
the association collapses to a single gene when Alb is dropped; and Alb **falls**,
whereas carry-over is additive and can only ever *raise* a foreign marker.
Cyp2e1, the other mover, goes *down* in iWAT but *up* in gWAT — contamination
cannot do that from one dissection.

**Sperm transcripts in 8 of 10 male gWAT libraries — epididymal carry-over.**
Gonadal fat sits against the epididymis. Balanced across genotype (4 CTL, 4 KD;
p = 0.57), and absent from both sex-interaction lists, so it cannot bias anything.

**Six unstable-estimate genes in the iWAT contrast are all sperm genes.** This is
the unified model's one real cost. DESeq2 estimates one dispersion per gene across
the whole dataset, so a gene that is zero in most libraries and enormous in eight
gWAT males gets a huge dispersion — which inflates its standard error in *every*
contrast, including iWAT where the gene is flat zero. They were never iWAT
findings. They are excluded from figures and retained in all tables.

**JS_06/07/09/10 correlate poorly with the rest.** All four are iWAT KAT8KD
females — the largest treatment-effect group. Correlation is computed against
everything else, mostly controls and the other depot, so the strongest-responding
group is *expected* to correlate least. Biology, not instrument failure.

### Two exclusion rules, and what they deliberately do not touch
Contaminant genes are removed from **volcano labels and heatmap selection only**,
by name (`CONTAMINATION_EXCLUDE`) and by statistic (`lfcSE > 2`). DEG counts, ORA
input and GSEA ranking are **untouched**. Deleting a real measurement to make a
figure read better is not defensible; stopping it from occupying a label slot on
the strength of a fold change driven by four libraries is.

---

## Methods paragraph (publication prose)

> Reads were quantified and analysed with DESeq2 (R/Bioconductor). Genes with at
> least 10 counts in at least 4 samples were retained. Two libraries (JS_08,
> JS_28) were excluded as outliers, leaving 38 samples. Counts were normalised by
> the median-of-ratios method and fitted with the cell-means design
> `~ 0 + GroupDepot`, where GroupDepot combines depot and genotype; depot-specific
> knockdown effects were extracted as contrasts. Differentially expressed genes
> were defined as FDR < 0.05 and |log2 fold change| > 1. Shrunken fold changes
> (ashr) were computed for ranking and visualisation only; differential
> expression calls and gene set ranking used the unshrunken estimates and Wald
> statistics respectively. Sexes were pooled after confirming sex-independence of
> the knockdown response by five criteria: a per-depot likelihood ratio test of
> Sex × Genotype (3.83% of genes in iWAT, 2.22% in gWAT at FDR < 0.05), true-null
> proportion π₀ = 0.87 and 0.85, a median interaction variance component of ~2%
> against 17–33% for genotype, sex-stratified directional concordance of 99.0%
> and 99.5%, and a model comparison showing that adding sex to the unified model
> reduced power. Over-representation analysis (clusterProfiler; GO:BP and KEGG)
> used all tested genes as background, gene sets of 5–500 genes, and GO:BP
> redundancy reduction at a semantic similarity of 0.7. Gene set enrichment
> analysis ranked all tested genes by the DESeq2 Wald statistic (gene sets
> 5–500 genes, FDR < 0.05, seed 12345). Contamination screening with
> tissue-restricted marker panels identified sperm/testis transcripts in eight of
> ten male gonadal libraries, consistent with epididymal carry-over during
> dissection; the signal was balanced across genotype (p = 0.57) and absent from
> the interaction results, and affected transcripts were retained in all tables
> but excluded from figure labelling, as was any gene with a standard error
> exceeding 2 on the log2 scale. All analyses used a fixed seed (12345).

---

## Quick reference

| Setting | Value |
|---|---|
| Prefilter | ≥ 10 counts in ≥ 4 samples |
| Samples | 40 sequenced, 38 analysed (JS_08, JS_28 removed) |
| Design | `~ 0 + GroupDepot` (unified, 38 samples) |
| DEG threshold | FDR < 0.05 and \|log2FC\| > 1 |
| GSEA ranking | DESeq2 Wald statistic |
| ORA background | all tested genes |
| Gene set size | 5–500 |
| GO:BP simplify | 0.7 (ORA), 0.5 on top 300 terms (GSEA) |
| NSL specificity | ≥ 2.0× = SPECIFIC, 1.5–2.0× = BORDERLINE |
| Figure exclusions | `CONTAMINATION_EXCLUDE` + `lfcSE > 2` |
| Seed | 12345 |
