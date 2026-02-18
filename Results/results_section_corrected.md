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

To investigate the transcriptomic consequences of KAT8 loss in adipocytes, we
performed bulk RNA-seq on iWAT and gWAT from adipocyte-specific KAT8 knockout
(AKO; Adipoq-Cre; Kat8^fl/fl) and floxed control (FL) mice (n = 4-5 per group
per depot). Both male and female mice were included with a balanced design, and
formal interaction testing confirmed that the KAT8 AKO transcriptional response
was largely sex-independent in both depots (<4% of genes with significant
Sex x Genotype interactions; >99% directional concordance; Supplemental Table X).

Principal component analysis revealed clear separation by both adipose depot and
genotype (PC1: 27%; PC2: 14.7%; Supplemental Fig. A). KAT8 AKO samples
segregated from FL controls in both depots, with the genotype effect being
particularly prominent in iWAT. Hierarchical clustering of the top DEGs
(|log2FC| > 1, FDR < 0.05) confirmed robust separation between genotypes in both
depots (Supplemental Fig. B). In iWAT, the majority of top DEGs were
downregulated in AKO relative to controls, whereas gWAT displayed both
upregulated and downregulated gene clusters, with a prominent block of
upregulated immune and inflammatory genes.

KEGG pathway analysis revealed depot-specific functional consequences of KAT8
loss (Fig. A). In iWAT, downregulated genes were strongly enriched in metabolic
pathways including PPAR signaling (p.adj = 8.3 x 10^-11), carbon metabolism,
branched-chain amino acid degradation, AMPK signaling, pyruvate metabolism, fatty
acid metabolism, regulation of lipolysis in adipocytes, and insulin signaling.
Upregulated genes in iWAT were enriched in ECM-receptor interaction (p.adj = 8.3
x 10^-8), neuroactive ligand-receptor interaction, cytokine-cytokine receptor
interaction, and PI3K-Akt signaling. In gWAT, upregulated genes were enriched in
cytokine-cytokine receptor interaction (p.adj = 2.0 x 10^-7), IL-17 signaling,
neutrophil extracellular trap formation, and chemokine signaling, reflecting a
pronounced inflammatory signature. Notably, no KEGG pathways reached significance
among gWAT downregulated genes, indicating that the metabolic gene suppression
observed in iWAT is depot-selective.

Volcano plot analysis confirmed that iWAT exhibited a broader range of fold
changes among DEGs (Fig. B). Downregulated genes in iWAT included those involved
in lipid metabolism (Cidea, Acacb, Ppara, Pck1, Acaca, Acox1, Ucp1),
glutathione metabolism (Gstk1, Mgst3, Gsta4), and peroxisomal function (Pxmp2,
Eci3). Upregulated genes included inflammatory mediators (Il1rn, Gdf15, Eda2r),
ECM components (Cd44, Sdc1), and calcium signaling genes (Slc8a1, Ryr2). In
contrast, gWAT displayed a more prominent upregulation of immune and
inflammatory genes, including chemokines (Ccl8, Ccl2, Cxcl2, Cxcl10, Ccl3),
cytokines (Il1b, Il12b), and immune receptors (Ctss, Itgax, Clec4d).

Examination of curated gene sets organized by functional category further
illustrated these depot-dependent changes (Fig. C). Mitochondrial and energy
metabolism genes, including electron transport chain components (Atp5pb, Ndufab1,
Ndufb10, Sdhd), mitochondrial-encoded genes (mt-Nd2, mt-Nd5), and mitochondrial
ribosomal proteins (Mrps10, Mrpl15), were consistently downregulated in KAT8 AKO
iWAT. A subset of these genes also trended downward in gWAT, although the gWAT
response was dominated by inflammatory and immune programs. Immune-related genes,
including Nckap1l, Myo1f, Cd84, Lat2, Fcer1g, Syk, Spi1, and Pycard, were
upregulated in KAT8 AKO adipose tissue in both depots. Lipid metabolism genes
(Lipa, Scarb2, Acsl3, Hilpda, Pltp) and ECM/remodeling genes (Ctss, Adam8,
Gpnmb, Ctsk, Timp1, Lgals3) were also differentially expressed between
genotypes. While the directionality of these changes was largely consistent
between depots, iWAT displayed more pronounced metabolic and thermogenic gene
suppression whereas gWAT exhibited a stronger inflammatory gene signature.

### Conclusion

Collectively, these data demonstrate that adipocyte-specific loss of KAT8 results
in widespread transcriptional dysregulation in both subcutaneous and visceral
white adipose depots. The predominant downregulation of genes governing
mitochondrial function, energy metabolism, and lipid catabolism, coupled with the
upregulation of immune and inflammatory gene programs, indicates that KAT8 is
required for the maintenance of normal adipocyte transcriptional identity. The
depot-dependent nature of these changes -- with iWAT exhibiting more pronounced
metabolic suppression and gWAT displaying a stronger inflammatory and chemokine
response -- is consistent with the known functional distinctions between
subcutaneous and visceral adipose tissue.

---

## Methods: RNA-seq Analysis

Sequencing reads were aligned and quantified to generate a gene-level count
matrix. Following quality assessment, two outlier samples were removed, yielding
38 samples for downstream analysis. Differential expression analysis was
performed using DESeq2 with a unified model across both depots, and
depot-specific KAT8 AKO effects were extracted via pairwise contrasts. Genes
with |log2FC| > 1 and FDR < 0.05 were considered differentially expressed. KEGG
pathway over-representation analysis and gene set enrichment analysis (fGSEA)
were performed on DEG lists and ranked gene lists, respectively. Both male and
female mice were included in each experimental group with a balanced design
(equal numbers of males and females per genotype per depot). To validate that
sexes could be combined, formal Sex x Genotype interaction testing was performed
using DESeq2 likelihood ratio tests within each depot. Significant interactions
were limited to 630 genes in iWAT (3.83%) and 374 genes in gWAT (2.22%;
FDR < 0.05), with >99% directional concordance between sex-stratified fold
changes (Supplemental Table X).

### Part 7.7 Validation Numbers (for Supplemental Table X)

| Metric                          | iWAT          | gWAT          |
|---------------------------------|---------------|---------------|
| Genes tested                    | 16,466        | 16,827        |
| Significant interactions (LRT)  | 630 (3.83%)   | 374 (2.22%)   |
| Pi0 (true null proportion)      | 0.87          | 0.85          |
| Sex-stratified Pearson r        | 0.727         | 0.605         |
| Directional concordance         | 99.0%         | 99.5%         |
| DEGs: ~ 0 + GroupDepot          | 3,022         | 1,159         |
| DEGs: ~ Sex + GroupDepot        | 2,777         | 1,089         |
| logFC correlation (+/- Sex)     | 0.939         | 0.967         |

### Discussion caveat (if reviewers ask about sex)

> While the overall KAT8 AKO response was concordant between sexes, PCA
> suggested a potentially larger effect magnitude in males, particularly in
> iWAT. A small subset of genes (2-4%) exhibited significant Sex x Genotype
> interactions that may warrant further investigation in follow-up studies.
