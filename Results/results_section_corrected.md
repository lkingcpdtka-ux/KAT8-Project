# KAT8 Results Section - Corrected and Improved (v2)

## Change Log (v2)

All changes below address the following review concerns:

1. **Merged sex numbers**: Verified correct (<4% interactions, >99% concordance).
2. **gWAT downregulated overstatement**: Removed the claim "indicating that the
   metabolic gene suppression observed in iWAT is depot-selective." Absence of
   significant KEGG enrichment among gWAT downregulated genes may reflect fewer
   downregulated DEGs in gWAT (reducing ORA power), not necessarily a true absence
   of metabolic suppression. Replaced with neutral language.
3. **Gene miscategorizations fixed**:
   - "inflammatory mediators (Il1rn, Gdf15, Eda2r)" → "immune signaling and
     stress-response genes (Il1rn, Gdf15, Eda2r)". Il1rn is anti-inflammatory
     (IL-1 receptor antagonist), Gdf15 is a stress-response cytokine, and Eda2r
     is a TNF receptor superfamily member / NF-kB signaling target.
   - "immune receptors (Ctss, Itgax, Clec4d)" → "immune effectors and receptors
     (Ctss, Itgax, Clec4d)". Ctss (Cathepsin S) is a lysosomal cysteine protease,
     NOT a receptor. Itgax (CD11c) and Clec4d are receptors.
   - Ctss inconsistency resolved: was listed as "immune receptor" in volcano
     section but "ECM/remodeling" in heatmap section. Ctss functions in both
     antigen processing (immune) and ECM degradation. Now labeled consistently:
     "immune effector" in volcano context, "ECM/remodeling" in heatmap context
     (where it clusters with other proteases and matrix genes).
4. **Misleading gWAT heatmap subset sentence**: Rewrote to clarify that apparent
   trends in gWAT mitochondrial genes were driven by sex-dependent variation rather
   than a consistent AKO genotype effect.
5. **Last heatmap claim**: Removed "thermogenic" — the heatmap genes (Atp5pb,
   Ndufab1, Ndufb10, Sdhd, mt-Nd2, mt-Nd5, Mrps10, Mrpl15) are ETC and
   mitoribosomal genes, not thermogenic. Ucp1 appears only in the volcano section.
   Softened "largely consistent between depots" to "broadly consistent
   directionality" for the genes that DO show agreement.
6. **Results vs. conclusion**: Removed interpretive/mechanistic claims from
   results ("depot-selective" implication). Kept descriptive language only.
   Moved functional interpretation to conclusion.
7. **Conclusion claims verified against literature**:
   - Visceral fat (gWAT) having greater inflammatory tone / macrophage infiltration:
     supported (Weisberg et al. 2003, Xu et al. 2003, Harman-Boehm et al. 2007).
   - Subcutaneous fat (iWAT) having greater thermogenic / lipid-oxidative capacity:
     supported (Wu et al. 2012, Frontini & Cinti 2010).
   - Added specific literature-aligned language about these depot distinctions.
   - Softened "indicates KAT8 is required" to "is consistent with a role for KAT8
     in maintaining" to avoid overclaiming causality from transcriptomic data.
8. **Dot plot order**: iWAT upregulated described first, then iWAT downregulated,
   then gWAT upregulated, then gWAT downregulated — matching figure panel order.

---

## Corrected Results Section (v2)

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
loss (Fig. A). In iWAT, upregulated genes were enriched in ECM-receptor
interaction (p.adj = 8.3 x 10^-8), neuroactive ligand-receptor interaction,
cytokine-cytokine receptor interaction, and PI3K-Akt signaling. Downregulated
genes in iWAT were strongly enriched in metabolic pathways including PPAR
signaling (p.adj = 8.3 x 10^-11), carbon metabolism, branched-chain amino acid
degradation, AMPK signaling, pyruvate metabolism, fatty acid metabolism,
regulation of lipolysis in adipocytes, and insulin signaling. In gWAT,
upregulated genes were enriched in cytokine-cytokine receptor interaction
(p.adj = 2.0 x 10^-7), IL-17 signaling, neutrophil extracellular trap formation,
and chemokine signaling. No KEGG pathways reached significance among gWAT
downregulated genes under our analysis thresholds, which may in part reflect the
smaller number of downregulated DEGs identified in gWAT relative to iWAT.

Volcano plot analysis confirmed that iWAT exhibited a broader range of fold
changes among DEGs (Fig. B). Downregulated genes in iWAT included those involved
in lipid metabolism (Cidea, Acacb, Ppara, Pck1, Acaca, Acox1, Ucp1),
glutathione metabolism (Gstk1, Mgst3, Gsta4), and peroxisomal function (Pxmp2,
Eci3). Upregulated genes included immune signaling and stress-response genes
(Il1rn, Gdf15, Eda2r), ECM components (Cd44, Sdc1), and calcium signaling genes
(Slc8a1, Ryr2). In contrast, gWAT displayed a more prominent upregulation of
immune and inflammatory genes, including chemokines (Ccl8, Ccl2, Cxcl2, Cxcl10,
Ccl3), cytokines (Il1b, Il12b), and immune effectors and receptors (Ctss, Itgax,
Clec4d).

Examination of curated gene sets organized by functional category further
illustrated these depot-dependent changes (Fig. C). Mitochondrial and energy
metabolism genes, including electron transport chain components (Atp5pb, Ndufab1,
Ndufb10, Sdhd), mitochondrial-encoded genes (mt-Nd2, mt-Nd5), and mitochondrial
ribosomal proteins (Mrps10, Mrpl15), were consistently downregulated in KAT8 AKO
iWAT. In gWAT, these mitochondrial genes did not show a consistent
genotype-driven pattern; observed variation was primarily attributable to sex
rather than a uniform AKO effect. Immune-related genes, including Nckap1l, Myo1f,
Cd84, Lat2, Fcer1g, Syk, Spi1, and Pycard, were upregulated in KAT8 AKO adipose
tissue in both depots. Lipid metabolism genes (Lipa, Scarb2, Acsl3, Hilpda, Pltp)
and ECM/remodeling genes (Ctss, Adam8, Gpnmb, Ctsk, Timp1, Lgals3) were also
differentially expressed between genotypes, with broadly consistent directionality
across depots. Overall, iWAT showed more pronounced suppression of mitochondrial
and metabolic gene programs, while gWAT was characterized by a stronger
inflammatory gene signature.

### Conclusion

Collectively, these data demonstrate that adipocyte-specific loss of KAT8 results
in widespread transcriptional dysregulation in both subcutaneous and visceral
white adipose depots. The predominant downregulation of genes governing
mitochondrial function, energy metabolism, and lipid catabolism in iWAT, coupled
with the upregulation of immune and inflammatory gene programs in both depots, is
consistent with a role for KAT8 in maintaining normal adipocyte transcriptional
programs. The depot-dependent nature of these changes -- with iWAT exhibiting
more pronounced metabolic and mitochondrial gene suppression and gWAT displaying
a stronger inflammatory and chemokine response -- aligns with the known
functional distinctions between subcutaneous and visceral adipose tissue,
including the greater thermogenic and lipid-oxidative capacity of subcutaneous fat
and the established susceptibility of visceral fat to inflammatory macrophage
infiltration. While the overall KAT8 AKO response was concordant between sexes,
a small subset of genes (2-4%) exhibited significant Sex x Genotype interactions
that may warrant further investigation in follow-up studies.

---

## Methods: RNA-seq Analysis

Sequencing reads were aligned and quantified to generate a gene-level count
matrix. Following quality assessment, two outlier samples were removed, yielding
38 samples for downstream analysis. Differential expression analysis was
performed using DESeq2 with a unified model across both depots, and
depot-specific KAT8 AKO effects were extracted via pairwise contrasts. Genes
with |log2FC| > 1 and FDR < 0.05 were considered differentially expressed. KEGG
pathway over-representation analysis was performed on DEG lists, and gene set
enrichment analysis (fGSEA) was performed on ranked gene lists. Both male and
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
