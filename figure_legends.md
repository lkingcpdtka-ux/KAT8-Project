# Figure Legends

## ORA KEGG Dot Plots

### Figure: KEGG Pathway Enrichment -- Upregulated Genes in gWAT (Dotplot_KEGG_gWAT_KD_vs_CTL_Up)

Over-representation analysis (ORA) of KEGG pathways enriched among significantly upregulated genes (FDR < 0.05, log2FC > 1.0) in gonadal white adipose tissue (gWAT) of KAT8 knockdown (KD) mice relative to controls (CTL). Each dot represents a single KEGG pathway. The x-axis shows the Gene Ratio, defined as the number of differentially expressed genes annotated to a given pathway divided by the total number of upregulated genes with KEGG annotations. Dot size corresponds to the number of genes (Count) mapping to each pathway, and dot color reflects the statistical significance as --log10(FDR), with warmer colors (red) indicating greater significance. Pathways are ordered by adjusted p-value (top = most significant). The most significantly enriched pathway was Cytokine-cytokine receptor interaction (Gene Ratio ~0.10, 29 genes). Other notable enriched pathways include Chemokine signaling pathway, IL-17 signaling pathway, and Neutrophil extracellular trap formation, consistent with an inflammatory signature in gWAT upon KAT8 loss. ORA was performed using a hypergeometric test with all DESeq2-tested genes as the background universe; pathway gene set sizes were restricted to 5--500 genes, and only pathways with FDR < 0.05 are shown (up to 15 pathways).

### Figure: KEGG Pathway Enrichment -- Downregulated Genes in iWAT (Dotplot_KEGG_iWAT_KD_vs_CTL_Down)

Over-representation analysis (ORA) of KEGG pathways enriched among significantly downregulated genes (FDR < 0.05, log2FC < --1.0) in inguinal white adipose tissue (iWAT) of KAT8 KD mice relative to CTL. Each dot represents a KEGG pathway, with the x-axis showing the Gene Ratio (number of DEGs in the pathway / total downregulated genes with KEGG annotations), dot size indicating gene count, and dot color representing --log10(FDR). The most significantly enriched pathway was PPAR signaling pathway (23 genes, including *Pparg*, *Ppara*, *Adipoq*, *Ucp1*, *Fabp4*, and *Plin1*), followed by Carbon metabolism, Propanoate metabolism, and Valine/leucine/isoleucine degradation. Other enriched pathways include AMPK signaling, Fatty acid metabolism, Insulin signaling pathway, Pyruvate metabolism, Regulation of lipolysis in adipocytes, and Glucagon signaling pathway, indicating broad suppression of lipid and energy metabolism programs in iWAT upon KAT8 loss. Statistical approach as described above.

### Figure: KEGG Pathway Enrichment -- Upregulated Genes in iWAT (Dotplot_KEGG_iWAT_KD_vs_CTL_Up)

Over-representation analysis (ORA) of KEGG pathways enriched among significantly upregulated genes (FDR < 0.05, log2FC > 1.0) in iWAT of KAT8 KD mice relative to CTL. Dot plot conventions as described above. The most significantly enriched pathways were ECM-receptor interaction (26 genes) and Cornified envelope formation (34 genes), followed by Neuroactive ligand-receptor interaction. Several structural and signaling pathways were enriched, including Dilated cardiomyopathy, Hypertrophic cardiomyopathy, PI3K-Akt signaling pathway (the largest pathway by gene count, ~40 genes), Integrin signaling, Cell adhesion molecule (CAM) interaction, and Cytokine-cytokine receptor interaction. These results suggest activation of extracellular matrix remodeling, cell adhesion, and epithelial differentiation programs in iWAT following KAT8 knockdown. Statistical approach as described above.

---

## GO:BP ORA Bar Plots

### Figure: GO Biological Process Enrichment -- Upregulated Genes in gWAT (ORA_GOBP_gWAT_KD_vs_CTL_Up)

Over-representation analysis (ORA) of Gene Ontology Biological Process (GO:BP) terms enriched among significantly upregulated genes (FDR < 0.05, log2FC > 1.0) in gWAT of KAT8 KD mice versus CTL. Bar length represents the Gene Ratio (number of DEGs annotated to a GO term / total upregulated genes with GO annotations), and bar color indicates statistical significance as --log10(FDR), with red denoting higher significance. Redundant GO terms were reduced using semantic similarity-based simplification (cutoff = 0.7). The top enriched terms are dominated by innate immune and inflammatory processes, including granulocyte migration, myeloid leukocyte activation, leukocyte mediated immunity, granulocyte chemotaxis, and leukocyte activation involved in immune response. Additional enriched terms include adaptive immune response, chemotaxis, defense response to bacterium, and chemokine-mediated signaling pathway. Up to 15 terms with FDR < 0.1 are displayed; gene set sizes restricted to 5--500.

### Figure: GO Biological Process Enrichment -- Downregulated Genes in iWAT (ORA_GOBP_iWAT_KD_vs_CTL_Down)

ORA of GO:BP terms enriched among significantly downregulated genes in iWAT of KAT8 KD mice versus CTL. Plot conventions as described above. The most significantly enriched term was fatty acid metabolic process (Gene Ratio ~0.10, --log10(FDR) > 25), followed by carboxylic acid catabolic process, small molecule catabolic process, lipid catabolic process, and cellular catabolic process. Other enriched processes include brown fat cell differentiation, adaptive thermogenesis, temperature homeostasis, triglyceride metabolic process, and regulation of lipid metabolic process, indicating a broad reduction of adipocyte-specific metabolic and thermogenic gene programs in iWAT upon KAT8 loss.

### Figure: GO Biological Process Enrichment -- Upregulated Genes in iWAT (ORA_GOBP_iWAT_KD_vs_CTL_Up)

ORA of GO:BP terms enriched among significantly upregulated genes in iWAT of KAT8 KD mice versus CTL. Plot conventions as described above. The most significantly enriched term was regulation of membrane potential, followed by muscle contraction and regulation of muscle contraction. Other enriched terms include signal release, modulation of chemical synaptic transmission, ion transport processes (sodium and potassium), skin development, sensory perception, and regulated exocytosis. These findings suggest ectopic activation of neuronal, muscular, and epithelial gene programs in iWAT upon KAT8 knockdown.

---

## Gene of Interest (GOI) Heatmaps (Heatmap_avg4col)

### Figure: Gene of Interest Heatmap -- gWAT KD vs CTL (Heatmap_avg4col_gWAT_KD_vs_CTL)

Heatmap displaying expression of curated genes of interest in gWAT, averaged across biological replicates within each genotype-by-sex group (4 columns: CTL Female, CTL Male, KAT8KD Female, KAT8KD Male). Expression values are row-wise z-score normalized VST (variance-stabilizing transformation) counts, clipped at +/- 2.0. Color scale ranges from blue (low relative expression, z-score = --2) through white (z-score = 0) to red/orange (high relative expression, z-score = 2). Column annotations indicate Genotype (CTL vs. KAT8KD) and Sex (F vs. M). Genes are organized into four functional categories derived from GO:BP fGSEA leading-edge analysis: (1) **Mitochondrial / Energy** -- includes mitochondrial complex subunits (*Atp5pb*, *Ndufab1*, *Ndufb10*, *Atp5f1d*, *Sdhd*, *Ndufb8*, *Atp5me*), mitochondrial-encoded genes (*mt-Nd2*, *mt-Nd5*), mitochondrial ribosomal proteins (*Mrps10*, *Mrpl15*, *Mrpl17*), and mitochondrial translation/assembly factors (*Nsun4*, *Uqcc2*, *Ptcd1*, *Mtg2*, *Dap3*); (2) **Immune / Inflammatory** -- includes myeloid and leukocyte activation genes (*Nckap1l*, *Myo1f*, *Cd84*, *Lat2*, *Fcer1g*, *Hmox1*, *Cd300a*, *Itgb2*, *Syk*, *Spi1*, *Tyrobp*, *Pycard*) from the leukocyte degranulation leading edge; (3) **Lipid Metabolism** -- includes lipid transport and processing genes (*Soat1*, *Lipa*, *Scarb2*, *Nus1*, *Acsl3*, *P2rx7*, *Hilpda*, *Pltp*, *Abcg1*); (4) **ECM / Remodeling** -- includes tissue remodeling genes (*Mdm2*, *Ctss*, *Adam8*, *Tcirg1*, *Gpnmb*, *Ctsk*, *Timp1*, *Lgals3*, *Cd44*, *Tnfaip8l2*). Rows are not clustered; gene order within each category is preserved from the fGSEA leading-edge ranking.

### Figure: Gene of Interest Heatmap -- iWAT KD vs CTL (Heatmap_avg4col_iWAT_KD_vs_CTL)

Heatmap displaying expression of the same curated gene set in iWAT, averaged across biological replicates within each genotype-by-sex group (4 columns). All plot conventions, gene groupings, normalization, and color scales are identical to the gWAT heatmap described above. This plot enables direct tissue comparison of KAT8 knockdown effects on genes of interest across the same functional categories (Mitochondrial/Energy, Immune/Inflammatory, Lipid Metabolism, ECM/Remodeling).

---

## GOI Volcano Plots

### Figure: Volcano Plot -- gWAT KD vs CTL (Volcano_gWAT_KD_vs_CTL)

Volcano plot of differential gene expression in gWAT comparing KAT8 KD to CTL. Each point represents a gene. The x-axis shows the log2 fold change (KD/CTL) and the y-axis shows --log10(FDR-adjusted p-value) from DESeq2. Genes are colored by significance category: orange indicates significantly upregulated genes (FDR < 0.05 and log2FC >= 1.0), cyan indicates significantly downregulated genes (FDR < 0.05 and log2FC <= --1.0), and grey indicates non-significant genes. Dashed horizontal and vertical lines mark the FDR and fold-change thresholds, respectively. The top 10 most significant upregulated genes, top 10 most significant downregulated genes, and top 10 genes by absolute fold change are labeled. *Acsm3* is highlighted with a yellow-filled point and bold label as a gene of particular interest. Gene labels are positioned using repulsive text placement to minimize overlap.

### Figure: Volcano Plot -- iWAT KD vs CTL (Volcano_iWAT_KD_vs_CTL)

Volcano plot of differential gene expression in iWAT comparing KAT8 KD to CTL. All plot conventions, thresholds, and labeling strategies are identical to the gWAT volcano plot described above. The x-axis shows log2 fold change (KD/CTL), the y-axis shows --log10(FDR), and genes are colored orange (upregulated), cyan (downregulated), or grey (non-significant) based on the same significance criteria (FDR < 0.05, |log2FC| >= 1.0). *Acsm3* is highlighted as a gene of interest.
