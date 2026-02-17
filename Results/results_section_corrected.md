# KAT8 Results Section - Corrected and Improved

## Data Audit Summary

The following issues were identified by cross-referencing the results text
against actual CSV data, generated figures, and analysis code.

### Errors Found

1. **Gene name error: "Lcp2" should be "Lat2"**
   - Text stated: "Nckap1l, Myo1f, Cd84, **Lcp2**, Fcer1g, Syk, Spi1, and Pycard"
   - Both the iWAT and gWAT heatmaps (Heatmap_avg4col_*.png) clearly show **Lat2**
     (Linker for Activation of T Cells Family Member 2), not Lcp2 (SLP-76)

2. **iWAT upregulated KEGG pathways were selectively reported**
   - Text listed: "ECM-receptor interaction, cornified envelope formation, and PI3K-Akt signaling"
   - PI3K-Akt signaling is ranked #15 (last significant, p.adj = 0.019)
   - Between #2 and #15, there are 12 other significant pathways including:
     - Neuroactive ligand-receptor interaction (#3, p.adj = 3.8e-06)
     - Dilated cardiomyopathy (#4, p.adj = 2.4e-05)
     - Hypertrophic cardiomyopathy (#5, p.adj = 2.4e-05)
     - Cytokine-cytokine receptor interaction (#12, p.adj = 0.013)
   - The muscle contraction/cardiomyopathy pathways reflect shared integrin and
     calcium channel genes rather than a true cardiac phenotype, but omitting all
     mention of them is misleading

3. **Heatmap gene lists are incomplete**
   - Immune/Inflammatory category in the actual heatmap also includes: Lat2 (not
     Lcp2), Hmox1, Cd300a, Itgb2, Tyrobp
   - Lipid Metabolism category also includes: Nus1, P2rx7, Abcg1
   - ECM/Remodeling category also includes: Mdm2, Tcirg1, Lgals3, Tnfaip8l2

4. **No gWAT downregulated KEGG ORA results exist**
   - No CSV file for ORA_kegg_gWAT_KD_vs_CTL_Down was generated, meaning no KEGG
     pathways reached significance among gWAT downregulated genes
   - The fGSEA GO:BP for gWAT does show downregulation of cilium assembly
     (NES = -1.96) and tRNA metabolic process (NES = -2.12)
   - This depot asymmetry is not addressed in the original text

5. **PCA variance values (27% and 14.7%) are unverifiable**
   - These are computed dynamically in the code (part1_main_analysis.R:569)
   - No saved PCA output exists in the repository to confirm these exact values
   - Recommend verifying from the actual PCA plot or saved run output

6. **Missing GO:BP ORA and fGSEA context**
   - The GO:BP barplots show important pathways not mentioned in the text:
     - iWAT Down: fatty acid metabolic process, brown fat cell differentiation,
       adaptive thermogenesis, temperature homeostasis
     - iWAT Up: regulation of membrane potential, muscle contraction
     - gWAT Up: granulocyte migration, myeloid leukocyte activation,
       leukocyte mediated immunity
   - fGSEA GO:BP reveals the strongest signal is mitochondrial dysfunction in iWAT:
     mitochondrial ATP synthesis (NES = -3.20), mitochondrial translation
     (NES = -3.17), respiratory chain complex I assembly (NES = -3.09)

---

## Corrected Results Section

### Results

**Adipocyte-specific loss of KAT8 produces depot-dependent transcriptional
dysregulation in white adipose tissue**

To characterize the transcriptomic consequences of KAT8 loss in adipocytes, we
performed bulk RNA-seq on inguinal white adipose tissue (iWAT) and gonadal white
adipose tissue (gWAT) from adipocyte-specific KAT8 knockdown (KD) and control
(CTL) mice (n = 5 per group per depot, both sexes). Principal component analysis
of variance-stabilized counts demonstrated clear separation of samples by adipose
depot along PC1 and by genotype along PC2, confirming that loss of KAT8 produces
a consistent transcriptional shift in both depots (Supplemental Fig. A).
Hierarchical clustering of the top differentially expressed genes (DEGs; |log2FC|
> 1, FDR < 0.05) revealed robust separation between CTL and KAT8 KD samples in
both iWAT and gWAT (Supplemental Fig. B).

KEGG over-representation analysis (ORA) revealed depot-specific functional
consequences of KAT8 loss (Fig. A). Among genes downregulated in iWAT, the most
significantly enriched pathways included PPAR signaling (p.adj = 8.3 x 10^-11),
carbon metabolism (p.adj = 3.1 x 10^-7), propanoate metabolism (p.adj = 3.7 x
10^-7), valine/leucine/isoleucine degradation (p.adj = 3.7 x 10^-7), AMPK
signaling (p.adj = 1.3 x 10^-6), fatty acid metabolism (p.adj = 4.8 x 10^-5),
regulation of lipolysis in adipocytes (p.adj = 8.7 x 10^-5), and insulin
signaling (p.adj = 1.6 x 10^-4). Metabolism of xenobiotics by cytochrome P450
(p.adj = 1.0 x 10^-6) and pyruvate metabolism (p.adj = 7.8 x 10^-6) were also
significantly enriched among downregulated genes. Upregulated genes in iWAT were
most significantly enriched in ECM-receptor interaction (p.adj = 8.3 x 10^-8) and
cornified envelope formation (p.adj = 8.5 x 10^-8), followed by neuroactive
ligand-receptor interaction (p.adj = 3.8 x 10^-6) and several additional
pathways including cytokine-cytokine receptor interaction and PI3K-Akt signaling.
In gWAT, upregulated genes were enriched in cytokine-cytokine receptor
interaction (p.adj = 2.0 x 10^-7), viral protein interaction with cytokine and
cytokine receptor (p.adj = 2.2 x 10^-7), IL-17 signaling (p.adj = 3.7 x 10^-3),
chemokine signaling (p.adj = 0.013), and neutrophil extracellular trap formation
(p.adj = 0.011), reflecting a pronounced inflammatory transcriptional signature.
Notably, no KEGG pathways reached significance among gWAT downregulated genes by
ORA, indicating that the metabolic gene suppression observed in iWAT was more
depot-selective.

Gene Ontology Biological Process (GO:BP) ORA further supported these findings
(Supplemental Fig. C). Downregulated genes in iWAT were enriched for fatty acid
metabolic process, carboxylic acid catabolic process, lipid catabolic process,
brown fat cell differentiation, adaptive thermogenesis, and temperature
homeostasis. Upregulated genes in gWAT were enriched for granulocyte migration,
myeloid leukocyte activation, leukocyte mediated immunity, and positive
regulation of cytokine production. Gene set enrichment analysis (fGSEA) on ranked
gene lists corroborated these results, identifying mitochondrial ATP synthesis
(NES = -3.20), mitochondrial translation (NES = -3.17), and respiratory chain
complex I assembly (NES = -3.09) as the most strongly downregulated gene sets in
iWAT, while leukocyte degranulation and myeloid cell activation were among the
most strongly upregulated gene sets in gWAT.

Volcano plot analysis confirmed that iWAT exhibited a broader dynamic range of
fold changes among DEGs, with numerous significantly downregulated genes involved
in lipid metabolism such as Cidea, Ucp1, Ppara, and Acsm5, while upregulated
genes included ECM-associated and integrin family members (Fig. B). In contrast,
gWAT displayed a more prominent upregulation of immune and inflammatory genes,
including Eda2r, Ctss, Il1rn, and Itgb2.

Examination of curated gene sets organized by functional category further
illustrated these depot-dependent changes (Fig. C). Mitochondrial and energy
metabolism genes, including electron transport chain components (Atp5pb, Ndufab1,
Ndufb10, Sdhd, Ndufb8), mitochondrial-encoded genes (mt-Nd2, mt-Nd5), and
mitochondrial ribosomal proteins (Mrps10, Mrpl15, Mrpl17), were consistently
downregulated in KAT8 KD samples in iWAT, and this pattern was also evident in
gWAT. Immune and inflammatory genes, including Nckap1l, Myo1f, Cd84, Lat2,
Fcer1g, Hmox1, Cd300a, Itgb2, Syk, Spi1, Tyrobp, and Pycard, were upregulated
in KAT8 KD adipose tissue in both depots. Lipid metabolism genes (Soat1, Lipa,
Scarb2, Nus1, Acsl3, Hilpda, Pltp, Abcg1) and ECM/remodeling genes (Ctss,
Adam8, Tcirg1, Gpnmb, Ctsk, Timp1, Lgals3, Cd44, Tnfaip8l2) were also
differentially expressed between genotypes. While the directionality of these
expression changes was largely consistent between iWAT and gWAT, the relative
magnitude differed by depot, with iWAT displaying more pronounced metabolic and
thermogenic gene suppression and gWAT exhibiting a stronger inflammatory gene
signature.

### Conclusion

Collectively, these data demonstrate that adipocyte-specific loss of KAT8 results
in widespread transcriptional changes in both subcutaneous and visceral white
adipose depots. The predominant downregulation of genes governing mitochondrial
function, energy metabolism, PPAR signaling, and lipid catabolism, coupled with
the upregulation of immune, inflammatory, and ECM remodeling gene programs,
indicates that KAT8 is required for the maintenance of normal adipocyte
transcriptional identity. The depot-dependent nature of these changes -- with
iWAT exhibiting more substantial suppression of metabolic and thermogenic
programs including brown fat cell differentiation and adaptive thermogenesis, and
gWAT displaying a more robust inflammatory response characterized by chemokine,
cytokine, and myeloid activation pathways -- is consistent with the known
functional and developmental distinctions between subcutaneous and visceral
adipose tissue. These findings identify KAT8 as a previously unrecognized
regulator of adipocyte gene expression and suggest that its loss promotes a
transcriptional environment favoring adipose tissue dysfunction.
