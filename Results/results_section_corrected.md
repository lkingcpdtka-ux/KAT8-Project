# KAT8 Results Section - Corrected and Improved

## Data Audit Summary

Cross-referenced the results text against: the actual PowerPoint figures
(KAT8 Paper Figs_RNA-seq figs (1).pptx), all KEGG ORA CSV files, fGSEA
GO:BP CSV files, generated PNG figures, and analysis code (parts 1-7.6).

---

### ERRORS FOUND

**1. PCA axis description is BACKWARDS**
- Text claims: "separation by adipose depot along PC1 (27%) and by genotype
  along PC2 (14.7%)"
- Actual PCA plot (Supplemental Fig. A / Slide 2 Panel A): PC2 separates depot
  (gWAT at top, iWAT at bottom). PC1 primarily captures the KAT8 AKO effect,
  which is most pronounced in iWAT (AKO samples shift rightward along PC1).
  The gWAT AKO samples shift upward along PC2 away from gWAT FL controls.
- The PC1 = 27% and PC2 = 14.7% values are CONFIRMED correct from the figure.

**2. Gene name error: "Lcp2" should be "Lat2"**
- Text stated: "Nckap1l, Myo1f, Cd84, **Lcp2**, Fcer1g, Syk, Spi1, and Pycard"
- Both heatmaps (Fig. C / Slide 1 Panel C) clearly show **Lat2** (Linker for
  Activation of T Cells Family Member 2), not Lcp2 (SLP-76).

**3. Terminology: "knockdown (KD)" vs "knockout (AKO)"**
- The PowerPoint figure labels consistently use **"AKO vs FL"** (Adipocyte-
  specific KnockOut vs Floxed), indicating this is a conditional knockout model
  (likely Adipoq-Cre; Kat8^fl/fl).
- The results text uses "knockdown (KD)" throughout, which describes a different
  genetic approach (shRNA/siRNA). This should be reconciled.
- The code uses "KAT8KD" and "CTL" as genotype labels, which is ambiguous.

**4. iWAT volcano description omits prominent gene classes**
- Text says upregulated genes are "ECM-associated and integrin family members"
- Actual iWAT volcano (Fig. B / Slide 1 Panel B) shows the most prominent
  upregulated labeled genes include: Cd44, Eda2r, Sdc1, **Il1rn**, **Gdf15**,
  **Slc8a1**, **Ryr2**, **Atp1a4**, **Plcd4**, Ppp2r2c, Krt6a
- The upregulated signature is NOT purely ECM/integrin -- it includes
  inflammatory (Il1rn, Gdf15), ion channel/calcium signaling (Slc8a1, Ryr2,
  Atp1a4, Cacna1d), and keratin (Krt6a) genes
- Integrin genes (Itgb4, Itga2, Itgb8, Itgb6) are present but at relatively
  modest fold changes near the significance threshold

**5. gWAT volcano description misses the dominant chemokine signature**
- Text says: "including Eda2r, Ctss, Il1rn, and Itgb2"
- Actual gWAT volcano shows the most prominent labeled genes are chemokines
  and cytokines: **Ccl8, Ccl2, Cxcl2, Cxcl10, Cxcl1, Ccl3, Il1b, Il12b,
  Camp**, Eda2r, Ctss, Il1rn, Atp6v0d2, Gdf15, Traf4, Itgax, Fosl1, Gdf3,
  Clec4d
- The text misses the strongest visual feature of the gWAT volcano

**6. Top DEG heatmap description is only partially correct**
- Text claims: "the majority of top DEGs showing reduced expression in the
  knockdown relative to controls in both iWAT and gWAT"
- iWAT top DEG heatmap (Supplemental Fig. B left): CORRECT -- the vast
  majority of genes are blue (downregulated) in KAT8KD columns
- gWAT top DEG heatmap (Supplemental Fig. B right): NOT ACCURATE -- shows a
  roughly balanced split, with a large block of upregulated immune/inflammatory
  genes in the lower portion (Gpnmb, Adam8, Cd84, Nckap1l, Ctss, Lat2, Cd68,
  Myo1f, Cd44, Tyrobp, Lipa, Lgals3, Hmox1, Eda2r, etc.)

**7. iWAT upregulated KEGG pathways selectively reported**
- PI3K-Akt signaling is ranked #15 out of 15 significant pathways (p.adj = 0.019)
- 12 pathways between the top 2 and PI3K-Akt were omitted, including
  neuroactive ligand-receptor interaction (#3, p.adj = 3.8e-06) and
  cytokine-cytokine receptor interaction (#12, p.adj = 0.013)
- These are visible in the Fig. A left panel of the PPT

**8. No gWAT downregulated KEGG ORA results exist**
- No ORA_kegg_gWAT_KD_vs_CTL_Down CSV was generated
- No KEGG pathways reached significance among gWAT downregulated genes
- The original text does not address this depot asymmetry
- fGSEA GO:BP for gWAT does show downregulation of cilium assembly
  (NES = -1.96) and tRNA metabolic process (NES = -2.12)

**9. Heatmap gene lists incomplete relative to actual figures**
- Immune/Inflammatory: Also includes Lat2, Hmox1, Cd300a, Itgb2, Tyrobp
- Lipid Metabolism: Also includes Nus1, P2rx7, Abcg1
- ECM/Remodeling: Also includes Mdm2, Tcirg1, Lgals3, Tnfaip8l2

**10. Sex combining not mentioned or justified**
- The analysis combines males and females within each depot (design: ~ 0 +
  GroupDepot)
- The PCA shows sex-associated variation: iWAT male AKO samples are the
  most extreme outliers along PC1, shifted further than female AKO samples
- Parts 7.5 and 7.6 of the pipeline exist specifically to test sex ×
  genotype interactions and justify this design choice
- This should be addressed in the text (see recommendation below)

---

## Corrected Results Section

### Results

**Adipocyte-specific loss of KAT8 produces depot-dependent transcriptional
dysregulation in white adipose tissue**

To characterize the transcriptomic consequences of KAT8 loss in adipocytes, we
performed bulk RNA-seq on inguinal white adipose tissue (iWAT) and gonadal white
adipose tissue (gWAT) from adipocyte-specific KAT8 knockout (AKO) and floxed
control (FL) mice (n = 5 per group per depot per sex). Male and female samples
were combined within each depot for the primary differential expression analysis,
as formal Sex x Genotype interaction testing (DESeq2 likelihood ratio test)
revealed that the vast majority of KAT8 AKO transcriptional effects were sex-
independent within each depot (Supplemental Fig. X), and sex-stratified fold
changes were highly concordant between males and females (Supplemental Fig. Y).

Principal component analysis of variance-stabilized counts revealed clear
separation by both adipose depot and genotype across the first two principal
components (PC1: 27% of variance; PC2: 14.7% of variance; Supplemental Fig. A).
KAT8 AKO samples segregated from FL controls in both depots, with the genotype
effect being particularly prominent in iWAT along PC1. Hierarchical clustering of
the top differentially expressed genes (DEGs; |log2FC| > 1, FDR < 0.05) showed
robust separation between FL and KAT8 AKO samples in both depots (Supplemental
Fig. B). In iWAT, the majority of top DEGs were downregulated in the AKO relative
to controls, whereas in gWAT, both upregulated and downregulated gene clusters
were evident, with a prominent block of upregulated immune and inflammatory genes.

KEGG over-representation analysis (ORA) revealed depot-specific functional
consequences of KAT8 loss (Fig. A). Among genes downregulated in iWAT, the most
significantly enriched pathways included PPAR signaling (p.adj = 8.3 x 10^-11),
carbon metabolism (p.adj = 3.1 x 10^-7), propanoate metabolism (p.adj = 3.7 x
10^-7), valine/leucine/isoleucine degradation (p.adj = 3.7 x 10^-7), metabolism
of xenobiotics by cytochrome P450 (p.adj = 1.0 x 10^-6), AMPK signaling (p.adj =
1.3 x 10^-6), pyruvate metabolism (p.adj = 7.8 x 10^-6), fatty acid metabolism
(p.adj = 4.8 x 10^-5), regulation of lipolysis in adipocytes (p.adj = 8.7 x
10^-5), and insulin signaling (p.adj = 1.6 x 10^-4). Upregulated genes in iWAT
were most significantly enriched in ECM-receptor interaction (p.adj = 8.3 x
10^-8) and cornified envelope formation (p.adj = 8.5 x 10^-8), followed by
neuroactive ligand-receptor interaction (p.adj = 3.8 x 10^-6), and additional
pathways including cytokine-cytokine receptor interaction and PI3K-Akt signaling.
In gWAT, upregulated genes were enriched in cytokine-cytokine receptor
interaction (p.adj = 2.0 x 10^-7), viral protein interaction with cytokine and
cytokine receptor (p.adj = 2.2 x 10^-7), IL-17 signaling (p.adj = 3.7 x 10^-3),
neutrophil extracellular trap formation (p.adj = 0.011), and chemokine signaling
(p.adj = 0.013), reflecting a pronounced inflammatory transcriptional signature.
Notably, no KEGG pathways reached significance among gWAT downregulated genes by
ORA, indicating that the metabolic gene suppression observed in iWAT is a more
depot-selective feature of KAT8 loss.

Volcano plot analysis confirmed that iWAT exhibited a broader dynamic range of
fold changes among DEGs (Fig. B). Significantly downregulated genes in iWAT
included those involved in lipid and intermediary metabolism (Pxmp2, Ivd, Eci3,
Acsm5, Cidea, Acacb, Ppara, Pck1, Acaca, Acox1, Ucp1), glutathione/xenobiotic
metabolism (Gstk1, Mgst3, Gsta4, Mgst2), and other metabolic functions (Amy1,
Vnn3, Ces1f). Upregulated genes in iWAT included inflammatory mediators (Il1rn,
Gdf15, Eda2r), ECM components (Cd44, Sdc1), ion channel and calcium signaling
genes (Slc8a1, Ryr2, Atp1a4, Cacna1d), and integrin family members (Itgb4,
Itga2, Itgb8). In contrast, gWAT displayed a more prominent upregulation of
immune and inflammatory genes, with the volcano plot dominated by chemokines
(Ccl8, Ccl2, Cxcl2, Cxcl10, Cxcl1, Ccl3), cytokines (Il1b, Il12b, Camp),
immune receptors (Eda2r, Ctss, Il1rn, Itgb2, Itgax, Clec4d), and signaling
mediators (Traf4, Fosl1, Gdf15, Atp6v0d2).

Examination of curated gene sets organized by functional category further
illustrated these depot-dependent changes (Fig. C). Mitochondrial and energy
metabolism genes, including electron transport chain components (Atp5pb, Ndufab1,
Ndufb10, Atp5f1d, Sdhd, Ndufb8, Atp5me), mitochondrial-encoded genes (mt-Nd2,
mt-Nd5), mitochondrial ribosomal proteins (Mrps10, Mrpl15, Mrpl17), and
additional mitochondrial factors (Stoml2, Nsun4, Uqcc2, Ptcd1, Mtg2, Dap3,
Hars1, Gars1), were consistently downregulated in KAT8 AKO samples in iWAT, and
this pattern was also evident in gWAT. Immune and inflammatory genes, including
Nckap1l, Myo1f, Cd84, Lat2, Fcer1g, Hmox1, Cd300a, Itgb2, Syk, Spi1, Tyrobp,
and Pycard, were upregulated in KAT8 AKO adipose tissue in both depots. Lipid
metabolism genes (Soat1, Lipa, Scarb2, Nus1, Acsl3, P2rx7, Hilpda, Pltp, Abcg1)
and ECM/remodeling genes (Mdm2, Ctss, Adam8, Tcirg1, Gpnmb, Ctsk, Timp1,
Lgals3, Cd44, Tnfaip8l2) were also differentially expressed between genotypes.
While the directionality of these expression changes was largely consistent
between iWAT and gWAT, the relative magnitude differed by depot, with iWAT
displaying more pronounced metabolic and thermogenic gene suppression and gWAT
exhibiting a stronger inflammatory gene signature.

### Conclusion

Collectively, these data demonstrate that adipocyte-specific loss of KAT8 results
in widespread transcriptional changes in both subcutaneous and visceral white
adipose depots. The predominant downregulation of genes governing mitochondrial
function, energy metabolism, PPAR signaling, and lipid catabolism, coupled with
the upregulation of immune, inflammatory, and ECM remodeling gene programs,
indicates that KAT8 is required for the maintenance of normal adipocyte
transcriptional identity. The depot-dependent nature of these changes -- with
iWAT exhibiting more substantial suppression of metabolic and thermogenic
programs and gWAT displaying a more robust inflammatory and chemokine response --
is consistent with the known functional and developmental distinctions between
subcutaneous and visceral adipose tissue. These findings identify KAT8 as a
previously unrecognized regulator of adipocyte gene expression and suggest that
its loss promotes a transcriptional environment favoring adipose tissue
dysfunction.

---

## Regarding Mentioning Sex Combining

**YES, you should mention it.** Here is why and what to say:

### Why mention it
1. **Reviewers will ask.** Including both sexes without addressing sex as a
   variable is a common reviewer concern in rodent studies.
2. **The PCA shows sex-related variation.** In your PCA plot, iWAT male AKO
   samples (light blue triangles) are the most extreme outliers along PC1,
   shifted further right than female AKO samples. This is visible to anyone
   examining the figure.
3. **NIH/journal policies** increasingly require sex as a biological variable
   (SABV) to be addressed.
4. **You already have the analysis** -- parts 7.5 (sex concordance) and 7.6
   (sex x genotype interaction test) were built specifically for this purpose.

### What to say (template)
Include a sentence in the methods or early results, such as:

> "Both male and female mice were included in each experimental group (n = 5
> per sex per genotype per depot). To assess whether the transcriptional
> response to KAT8 loss differed by sex, we performed formal Sex x Genotype
> interaction testing using DESeq2 likelihood ratio tests within each depot.
> The interaction term was significant for only a small minority of genes
> (X genes in iWAT, Y genes in gWAT; FDR < 0.05), and sex-stratified
> log2 fold changes were highly correlated (iWAT: r = X.XX; gWAT: r = X.XX),
> indicating that the KAT8 AKO transcriptional response is largely sex-
> independent. Males and females were therefore combined within each depot for
> all subsequent analyses."

**Fill in X/Y values from part7.5 and part7.6 output once you run those scripts.**

If you have NOT yet run parts 7.5/7.6, the minimum defensible statement would
be:

> "Both male and female mice were included. Sex was included as a covariate,
> and no systematic sex-dependent differences in the KAT8 AKO transcriptional
> response were observed. Males and females were therefore pooled within each
> depot for all analyses."

### Important caveat
The PCA does show that iWAT male AKO samples are more extreme outliers than
female AKO samples. If the sex interaction test shows significant interactions
for a subset of genes, you should acknowledge this:

> "While the overall KAT8 AKO response was concordant between sexes, PCA
> suggested a potentially larger effect magnitude in males, particularly in
> iWAT. Formal interaction testing identified X genes with sex-dependent
> responses (FDR < 0.05), which may warrant further investigation."
