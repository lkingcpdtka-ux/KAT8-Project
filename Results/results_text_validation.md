# Validation notes for drafted RNA-seq Results text

This note checks whether the drafted narrative is supported by the exported result tables and presentation assets currently in `Results/`.

## Overall assessment
- The core interpretation is directionally consistent with the data exports:
  - iWAT downregulated genes are strongly enriched for metabolic pathways (PPAR, carbon, AMPK, pyruvate, fatty acid metabolism).
  - iWAT upregulated genes are enriched for ECM and signaling categories.
  - gWAT upregulated genes are enriched for cytokine/immune-inflammatory pathways.
- Several numeric pathway p.adjust values in the draft exactly match the KEGG ORA outputs.

## Specific points that are well-supported
- iWAT-down KEGG pathways and p.adjust values (PPAR, carbon metabolism, propanoate, BCAA degradation, xenobiotic metabolism, AMPK, pyruvate, fatty acid metabolism, lipolysis regulation, insulin signaling) are present in `ORA_kegg_iWAT_KD_vs_CTL_Down...csv`.
- iWAT-up KEGG pathways (ECM-receptor interaction, cornified envelope formation, neuroactive ligand-receptor interaction, cytokine-cytokine receptor interaction, PI3K-Akt signaling) are present in `ORA_kegg_iWAT_KD_vs_CTL_Up...csv`.
- gWAT-up KEGG pathways and p.adjust values (cytokine-cytokine receptor interaction, viral protein interaction with cytokine/cytokine receptor, IL-17, NET formation, chemokine signaling) are present in `ORA_kegg_gWAT_KD_vs_CTL_Up...csv`.
- "No significant KEGG pathways among gWAT downregulated genes" is consistent with the absence of a `ORA_kegg_gWAT_KD_vs_CTL_Down...csv` output in the Results directory.

## Points that should be softened/clarified
- PCA variance percentages (PC1 27%, PC2 14.7%) and specific clustering statements cannot be directly verified from the currently exported files listed in `Results/`; the PPT contains mostly figure panels and minimal embedded text.
- Claims that mitochondrial/energy programs are "also evident" as downregulated in gWAT are weaker than in iWAT based on available fgsea exports; gWAT fgsea is predominantly immune-up and includes mixed non-immune signals.
- Detailed volcano gene lists (e.g., exact genes called out as significantly changed) should be cross-checked against the exact DEG table used to build the plotted volcano labels before finalizing manuscript text.

## Recommendation
Use the current draft with minor edits to reduce over-specificity where direct support is not visible in the exported tables (especially PCA percentages and strength of gWAT mitochondrial suppression), and cite supplemental figures/tables explicitly for those claims.
