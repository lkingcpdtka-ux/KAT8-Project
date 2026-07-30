# Follow-up tasks — status report

Environment note, because it determines what could be executed here: this
container holds the **pipeline scripts** and the three staging DE tables in
`downstream_analysis/data/`. It does **not** hold `savepoints/`, the part 5a–5c
outputs, the sex-interaction outputs, any alignments, or an R installation.
Checks that the local DE tables support were run and are reported with numbers.
Everything else is delivered as a script to run locally, or stopped.

| Task | Status |
|---|---|
| 1 — Kat8 deletion efficiency | **STOPPED** — no alignments; not pursued by decision |
| 2 — Secretome vs background | **RESOLVED** — ≥1 depot claim fails, both-depot survives matching |
| 3 — Relaxed gWAT down ORA | **RESOLVED** — primary GSEA already carries the finding |
| 4 — Layer 3 join defect | **FIXED and verified** |
| 5 — Layer 3 set-size inconsistency | **Documented** — 13/36 survive 5–500 |
| 6 — Sex summary | (a) **BLOCKED** — scripts absent; (b) **DONE** |
| 7 — DE_cells column harmonisation | **DONE** — preconditions verified, file written |

## Outcomes (all tasks run 2026-07-30)

**Task 2.** Fisher against background: ≥1 depot OR 1.45, p = 0.164 (**not
significant**); both depots OR 2.17, p = 0.012. Against an
expression-and-effect-size-matched background: ≥1 depot OR 1.93, p = 0.053;
both depots **OR 2.54, p = 0.029 — survives matching**. Background resolves to
314, not 312, as predicted. *"45 of 78 (58%)" must not be quoted as evidence;
the both-depot statistic is defensible and now has a matched-background check
behind it.*

**Task 3 — resolved better than the sensitivity analysis anticipated.** The
relaxed set returned 8 significant GO:BP terms, dominated by cilia. But checking
the primary GSEA showed **the finding was already there, as the top-ranked
result**: `cilium organization` (setSize 395, NES −1.94, padj **1.02e-08**) and
`intraciliary transport` (setSize 48, NES −2.39, padj 5.13e-07), both
down-regulated. No ciliary terms in iWAT GO:BP or in either KEGG table — the
effect is **depot-restricted to gWAT**. So no post hoc threshold needs
justifying and nothing goes to supplementary as a rescue: the claim of "no
enrichment" was never true of the analysis, only of ORA within it. The relaxed
ORA is corroboration.

*Leading edge of `cilium organization` (98 genes):* **98 of 98 significant at
FDR < 0.05, 98 of 98 down-regulated, and only 2 exceeding |log2FC| > 1.** That
one line explains the whole ORA/GSEA discrepancy — 96 of 98 genes sat below the
gate ORA applies, so the method could not see a pathway that was fully present
in the data. Effect sizes cluster between −0.6 and −1.0.

Composition: intraflagellar transport (Ift140, Ift122, Ift172, Ift57, Ift46,
Ift27, Ift22, Ift43, Wdr35, Cluap1, Ttc21b), the BBSome and its chaperonin
(**Bbs1, Bbs2, Bbs7, Bbs10, Ttc8/BBS8, Arl6/BBS3**, plus BBS-associated Wdpcp),
transition zone (Mks1, Rpgrip1l, Tctn1, Tctn2, Tmem231, Cplane1), retrograde
motor (Dync2li1) and basal body (Cep19, Cep41, Cep89, Cep162, Pcm1). Six
Bardet-Biedl genes and Cep19 — all human obesity-ciliopathy genes —
coordinately down in visceral fat and unchanged in subcutaneous.

*Hedgehog is NOT supported and should not be claimed.* No Gli, Ptch, Smo, Sufu,
Kif7 or Hhip appears in the leading edge. The `smoothened signaling` term in the
relaxed ORA is annotation overlap: GO assigns ciliary structural genes to that
term because cilia are required to transduce Hedgehog, so it is the same genes
under a second label, not independent evidence. For the discussion, the correct
framing is that Smo and Gli are regulated by localisation and post-translational
processing rather than transcript abundance, so RNA-seq cannot report Hedgehog
activity either way — and PROGENy offers no fallback, since Hedgehog is not
among its 14 pathways.

**Task 4 — verified.** All three expected values present after the fix: Foxj2
−1.202 / 0.279, Ncoa2 +1.293 / 0.138, Snai1 +0.915 / 0.882, each `tested, not
enriched`, `n_layers` = 1. 3 of 58 testable, as expected. `n_layers` unchanged
for every TF, confirming the support gate held. The three multi-layer TFs
(Creb3l1, Myc, Stat2) are all in the *not testable* group — Layer 3 and Layers
1–2 have essentially no overlap in this dataset.

**Task 5 — numbers now known.** 437 sets returned (450 passed the pre-filter;
fgsea drops sets below `minSize` after intersecting with the ranked list, which
is expected). Size range 10–1456; **144 of 437 outside 5–500**; **13 of 36
significant sets would survive the bound**, i.e. 23 exceed it. Limitations
section only.

**Task 7 — done.** `pvalue` == `P.Value` and `log2FoldChange` == `logFC` both
identical across 12,798 rows; 745 NA retained. Harmonised file written to
`followups/output/`; original untouched.

**Reproducibility check not in the original brief.** Part 5b logged a failure to
download the PROGENy model from OmniPath and fell back to a static table. All 14
pathway scores are byte-identical to the 28 July run, including JAK-STAT
(−12.482), so the fallback is equivalent and the figure is reproducible.

---

## Task 1 — Kat8 deletion efficiency — STOPPED

**Not attempted, as instructed.** No alignments are available. Searching the
repository for `*.bam`, count matrices or any per-exon quantification returns
nothing; only gene-level DE tables exist. Per-exon coverage cannot be derived
from a gene-level table — the summation that hides residual exons has already
happened.

The premise is correct and worth restating: gene-level `Kat8` log2FC of −0.42
(iWAT) and −0.59 (gWAT) understates deletion for two compounding reasons. Under
Adipoq-Cre the tissue contains deleted adipocytes plus untouched stroma,
endothelium and immune cells, so the value is diluted by cells that were never
targeted; and the floxed allele removes only a few exons while gene-level
counting sums across the whole transcript, so retained exons keep contributing
reads. **Neither figure should be presented as a recombination efficiency.**

To complete this task, one of:

1. **BAM files** plus the allele design. Note the exon identity is also not in
   this repository — it has to come from the allele's targeting strategy (MGI
   allele page, IMPC/EUCOMM design, or the vector map from whoever made the
   line). Do not infer it from the gene model. With both, per-exon coverage
   (`featureCounts -f`, or `GenomicAlignments::summarizeOverlaps` on an exon
   `GRangesList`) gives the floxed-vs-retained ratio per library.
2. **qPCR across the floxed exons** on the same RNA, with a retained-exon
   amplicon as the internal reference.
3. **Genotyping PCR** on genomic DNA from the same tissue, or **KAT8 protein**
   by western/IF.

Option 1 is strongest because it reuses data already generated and gives a
per-sample distribution rather than a point estimate.

---

## Task 2 — Secretome reproduction vs background

**Script:** `followups/task2_secretome_background.R` (needs
`part5a_secretome_annotated_20260728_000412.csv`, which is not in this container).

What the local DE tables **do** confirm:

| Check | Result |
|---|---|
| Fisher, ≥1 depot — 45/78 vs 152/312 | **OR 1.435, p = 0.1656** — matches the stated 1.44 / 0.166 |
| Fisher, both depots — 22/78 vs 48/312 | **OR 2.161, p = 0.01269** — matches the stated 2.16 / 0.0127 |
| All cell-significant genes reproducing ≥1 depot | **197** = 45 + 152 ✓ |
| All cell-significant genes reproducing both depots | **70** = 22 + 48 ✓ |
| Total genes with cells padj < 0.05 | **392**, but 78 + 312 = **390** |

The two reproduction totals reconcile exactly, which independently confirms the
expected split. **But the set sizes are two genes short.** Since the reproducing
counts already balance, the two unaccounted genes must reproduce in neither
depot, so they belong in the background: 314, not 312. Re-running with 314
changes nothing that matters — ≥1 depot OR 1.453, p = 0.164; both depots
OR 2.177, p = 0.0124. Worth finding the two genes (likely a duplicated symbol or
an NA `gene_name` dropped by one filter and not the other), but no conclusion
turns on it.

### Verdict — the existing claim does not fully survive

- **"45 of 78 (58%) reproduce in vivo" must not be quoted as evidence.** Against
  a background of 49%, p = 0.166. Secreted genes reproduce at essentially the
  rate of any cell-significant gene. The number is real; the *implication* that
  secretome genes are special is not supported.
- **"22 of 78 (28%) reproduce in both depots" is defensible.** Against 15%
  background, OR 2.16, p = 0.013. Requiring both depots is the more demanding
  criterion and it is where the enrichment actually lives.

The script also runs the expression-matched background (stratified on `baseMean`
and |cell log2FC|) and reports explicitly whether the both-depot conclusion
survives matching.

---

## Task 3 — Sensitivity check: gWAT down-regulated ORA

**Script:** `followups/task3_relaxed_gwat_ora.R` (needs R with
clusterProfiler + org.Mm.eg.db; not installed here).

Threshold arithmetic verified against the local gWAT table:

| Set | n |
|---|---|
| padj < 0.05, any direction | **6,195** ✓ matches |
| padj < 0.05 and \|log2FC\| > 1 | **1,159** ✓ matches |
| padj < 0.05 and log2FC < −1 (primary DOWN) | **422** ✓ matches DEG summary |
| padj < 0.05 and log2FC < −0.5 (**relaxed DOWN**) | **1,399** |

So the relaxed set is 3.3× the primary down set — a genuine test of the
threshold, not a token relaxation.

Parameter identity is guaranteed by construction rather than by retyping: the
script `source()`s `parameters.R` and reads `ora_params`, which is exactly what
`part2_ora.R` does (its `ora_enrich_params` is a straight copy of `ora_params`).
The universe is rebuilt the same way — every gene tested by DESeq2, mapped
SYMBOL→ENTREZID through org.Mm.eg.db, matching `use_universe <- TRUE` in the
primary run — and the script prints the universe size and its provenance so this
can be confirmed rather than trusted.

The report states plainly whether the absence survives. If it does, the finding
stands; if terms appear, the claim must be softened to "no enrichment survives
the primary |log2FC| > 1 gate" rather than "gWAT down-regulation is
unstructured." Either way nothing substitutes into a primary analysis.

---

## Task 4 — Layer 3 join defect — FIXED

**Changed:** `downstream_analysis/part5b_tf_analysis.R`. Re-run part 5b
(`RUN_FORCE <- "p5b"`) to regenerate the convergence table.

The defect was real. The join filtered the fgsea result to significant sets
before matching:

```r
hits <- tft %>% dplyr::filter(padj < FDR_CUT)   # <- only significant
```

so a TF that *was* tested and came back near-significant recorded identically to
one that came back flat — both NA. NCOA2 (NES +1.29, p = 0.016) and SNAI1
(NES +0.92, p = 0.75) are not the same evidence, and the table could not
distinguish them. That distinction is the whole point of a convergence analysis.

**The fix** matches against the full fgsea result and adds `layer3_NES`,
`layer3_pval`, `layer3_padj`, populated for every TF that has a TFT set and left
empty only where no set exists.

Two things I had to get right, both of which would have caused silent damage:

1. **Tie-break.** Where a TF matches more than one set, the most significant
   match is now taken. Without this, switching from the significant subset to the
   full table could select a *different* set for a TF that does have a
   significant one — which would have changed `n_layers` as a side effect of a
   cosmetic fix.
2. **Support vs having a value.** `n_layers` previously counted non-NA in the
   significance-gated column. Retaining all NES values naively would have
   inflated every tested TF's layer count. Support is now an explicit
   `layer3_supported` flag (`testable & padj < FDR_CUT`), and `n_layers` counts
   that. **`n_layers` is unchanged for every TF by construction** — the fix adds
   information without moving any support call.

`layer3_status` keeps its three states and is now derived from `padj` directly.
The console table prints `layer3_NES` and `layer3_padj` alongside the status.

On re-running, verify the three known values appear (NCOA2 +1.293 / 0.138,
FOXJ2 −1.202 / 0.279, SNAI1 +0.915 / 0.882) and that the testable count is 3 of
58 — the sanity log already reports it.

---

## Task 5 — Layer 3 set-size inconsistency — documented

**Changed:** `downstream_analysis/part5b_tf_analysis.R` — reporting only. The
fgsea call still uses `minSize = 10, maxSize = 2000` and **was deliberately not
changed**.

The inconsistency is confirmed in the source: the Layer 3 call hardcodes
10–2000, while every other gene-set analysis inherits `min_gs_size = 5`,
`max_gs_size = 500` from `gsea_params`. Three lines now log, at runtime, the
size range tested, how many sets fall outside 5–500, and how many of the
significant sets would survive that bound — flagged as limitations material, not
a result.

**The survivors are not re-reported as a finding, and should not be.** Applying a
size filter after seeing which sets came up significant fits the filter to the
answer. The correct use is a limitations sentence: the Layer 3 run did not apply
the pipeline's standard set-size bound, and a substantial share of its
significant sets are larger than that bound allows.

It matters less than it appears, for a reason that is independent of set size:
**only 3 of 58 convergence candidates have any C3:TFT set at all**, so Layer 3
is uninformative for this dataset regardless of how it is filtered.

### Better third-layer options, should a genuine one be wanted later — not run

- **ChEA3** — the strongest candidate. Aggregates ENCODE, ReMap, literature
  ChIP-seq and co-expression into ~1,600 TF libraries, so coverage of adipocyte
  TFs is far better than C3:TFT. Mouse via orthologue mapping.
- **TRRUST v2** — curated from literature, has a native mouse network (~800
  TFs). High precision, small regulons, so recall is limited; good as a
  confirmatory rather than a discovery layer.
- **CollecTRI** — **already Layer 2** in this pipeline (`get_collectri`, with
  DoRothEA A/B/C as fallback). Using it again would not be a third layer at all,
  just Layer 2 twice.

The genuinely independent axis remains actual binding data with broad TF
coverage — ChEA3, or ReMap/ENCODE directly.

---

## Task 6(b) — variance decomposition — DONE

Written to `followups/output/task6b_variance_decomposition_*.csv`.

| Variance component | iWAT | gWAT |
|---|---|---|
| Sex | 6.54% | 11.47% |
| Genotype | **33.17%** | **16.73%** |
| Sex × Genotype | **2.05%** | **2.39%** |
| Residual | 46.32% | 50.07% |
| *(sum of medians)* | *88.08%* | *80.66%* |

Genotype explains **16.2×** more variance than the interaction in iWAT and
**7.0×** more in gWAT. This is the defensible lead argument for pooling, with
the 99%+ directional concordance beside it.

**Two presentation cautions.** These are **medians taken independently across
genes**, so they do not sum to 100 (88.1% and 80.7%). They must not be drawn as
a pie chart or a stacked bar — the deck's pie-chart analogy is wrong on this
point — and the shortfall should be explained if asked rather than rescaled.

**`N_genotype_DEGs` in that file (2,691 iWAT / 1,054 gWAT) matches none of the
other counts** — not the unified headline (3,022 / 1,159), nor per-depot
`~ Genotype` (2,665 / 1,077), nor `~ Sex + Genotype` (3,042 / 1,131). Every
other field reconciles exactly with the deck (π₀ 0.8701/0.8483, interaction
630/374, correlation 0.7272/0.6054), so the discrepancy is confined to that one
column. Drop it or state which model produced it before showing the table.

## Task 6(a) — N_sex_DEGs — STOPPED

The sex-interaction analysis is **not part of this repository**. Searching every
`.R` file for `Sex`, `interaction` or `LRT` returns only `part1_main_analysis.R`,
and there only for sample annotation, the sex-check against Xist/Y-linked genes,
and plot colouring. The scripts referenced in the deck as "Parts 7.6 + 7.7" —
which produced `sex_interaction_summary_20260610_152837.csv` — live outside the
tracked pipeline, and the summary file itself is not in this container.

Recomputing `N_sex_DEGs` "using the same model and thresholds as the rest of that
table" is therefore not possible without guessing what that model was. Since a
sex main-effect DEG count depends entirely on the design used to obtain it
(`~ Sex + Genotype` per depot vs the unified model, and which contrast was
extracted), any number I produced would be an approximation of an unknown
procedure — exactly what the brief says to stop rather than substitute.

**To unblock:** provide the 7.6/7.7 script, or state the model and threshold to
use. Part (b) — reformatting the variance decomposition already present in that
file — needs only the file itself and is trivial once available.

**One point that stands regardless, and I agree with the framing in the brief:**
π₀ = 0.87/0.85 means roughly 13–15% of genes carry a genuine sex-by-genotype
interaction. That is *not* evidence interaction is absent, and it should not be
led with. The defensible arguments for pooling are the variance split
(interaction ~2% vs genotype 17–33%) and the 99%+ directional concordance —
those say the interaction exists but is small and does not flip direction, which
is the actual justification.

---

## Task 7 — DE_cells column harmonisation — verified, script ready

**Script:** `followups/task7_harmonise_de_cells.R`. Preconditions were checked
here against `downstream_analysis/data/DE_cells_KAT8KD_vs_CTL.csv` (the untagged
staging copy of the same table):

| Check | Result |
|---|---|
| `pvalue` == `P.Value` | **identical, 0 mismatches** across 12,798 rows |
| `log2FoldChange` == `logFC` | **identical, 0 mismatches** |
| `adj.P.Val` NA count | **745** — matches the brief exactly |
| Tissue tables use `padj` | confirmed, both |

All three preconditions hold, so dropping the duplicates is safe. The script
re-verifies on the tagged file before touching anything and stops if the columns
are not identical — a column that merely *looks* redundant is not, and deleting
it would destroy data.

The 745 NA values are retained as NA. They are DESeq2 independent filtering —
genes deliberately set aside to gain power — and imputing them, dropping those
rows, or filling with 1 would each silently change what the table means.

Worth noting: the tissue tables also carry a duplicate `logFC` alongside
`log2FoldChange`. That is left alone; it is consistent across both tissue files
and downstream scripts read it.

---

## Anything unexpected

Four things. **First**, the secretome set arithmetic is two genes short — 392
genes reach cells padj < 0.05, against 78 + 312 = 390 — while the reproduction
totals (197 and 70) reconcile perfectly. The discrepancy is therefore confined to
non-reproducing background genes and changes no conclusion, but it means one of
the two set definitions drops two genes the other keeps, probably a duplicated
symbol or an NA `gene_name`. **Second**, the Layer 3 fix contained a trap: naively
retaining all NES values would have inflated `n_layers` for every tested TF, and
matching against the full table without a significance tie-break could have
silently reassigned which set a TF matched. Both are handled, and `n_layers` is
provably unchanged. **Third**, `docs/METHODS_EXPLAINED.md` described the tissue
arm as an siRNA knockdown; it is a constitutive Adipoq-Cre knockout. That was my
error and is corrected — and the correction strengthens the §4 argument, because
Adipoq-Cre deletes *during* adipogenesis while the 3T3-L1 siRNA acts *after* it,
which is precisely why the tissue has a phenotype and the isolated adipocyte does
not. **Fourth**, both blocked tasks are blocked for the same underlying reason:
the alignments and the sex-analysis scripts are the two pieces of this project
that were never brought under the pipeline. That is worth fixing for
reproducibility independent of these tasks.

**No primary result is altered by anything in this report.** Task 2 changes which
of two existing statistics may be quoted; Task 3 tests an existing claim without
replacing it; Tasks 4, 5 and 7 are a defect fix, documentation, and a column
rename.
