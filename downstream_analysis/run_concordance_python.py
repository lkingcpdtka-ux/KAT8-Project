#!/usr/bin/env python3
"""Cell <-> tissue concordance (part4b), Python port. Direct/indirect & cell-autonomy.
Threshold-free on the cell side. No external data needed."""
import os, numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "downstream_analysis/data"
OUT  = "downstream_analysis/results"
os.makedirs(OUT, exist_ok=True)
np.random.seed(12345)

FILES = {
    "cells": "DE_cells_KAT8KD_vs_CTL.csv",
    "iWAT":  "DE_tissue_iWAT_KD_vs_CTL.csv",
    "gWAT":  "DE_tissue_gWAT_KD_vs_CTL.csv",
}
MARK_METAB = ["Ucp1","Cidea","Ppara","Pparg","Adipoq","Acox1","Pck1","Acacb","Cd36","Plin1"]
MARK_IMM   = ["Ccl2","Ccl8","Cxcl2","Cxcl10","Il1b","Il12b","Itgax","Ctss","Cd44","Gdf15"]
ANCHOR     = ["Kat8"]
TISSUE_FDR = 0.05

def load(path):
    df = pd.read_csv(path)
    gene = df["gene_name"].astype(str) if "gene_name" in df.columns else df.iloc[:,0].astype(str)
    padj = df["padj"] if "padj" in df.columns else df.get("adj.P.Val")
    out = pd.DataFrame({
        "gene": gene.values,
        "stat": pd.to_numeric(df["stat"], errors="coerce").values,
        "lfc":  pd.to_numeric(df["log2FoldChange"], errors="coerce").values,
        "padj": pd.to_numeric(padj, errors="coerce").values,
    })
    out = out.dropna(subset=["stat"]).sort_values("padj")
    out = out.drop_duplicates("gene", keep="first").reset_index(drop=True)
    return out

de = {k: load(os.path.join(DATA, v)) for k, v in FILES.items()}
for k, d in de.items():
    print(f"[INFO] {k:6s}: {len(d)} genes, {int((d.padj<0.05).sum())} at padj<0.05")

def gsea_es(rank_genes, gene_set, stats_by_gene, n_perm=2000):
    """Weighted KS enrichment score (GSEA) + permutation p (gene-label shuffle)."""
    genes = list(rank_genes)
    in_set = np.array([g in gene_set for g in genes])
    if in_set.sum() < 5:
        return np.nan, np.nan, np.nan
    w = np.abs(np.array([stats_by_gene[g] for g in genes]))
    def es(mask):
        hit = mask * w
        nr = hit.sum()
        if nr == 0: return 0.0
        miss = (~mask).astype(float)
        run = np.cumsum(hit/nr - miss/max((~mask).sum(),1))
        return run[np.argmax(np.abs(run))]
    obs = es(in_set)
    k = int(in_set.sum())
    null = np.empty(n_perm)
    idx = np.arange(len(genes))
    for i in range(n_perm):
        m = np.zeros(len(genes), bool); m[np.random.choice(idx, k, replace=False)] = True
        null[i] = es(m)
    if obs >= 0:
        p = (1 + (null >= obs).sum()) / (1 + n_perm)
        nes = obs / (np.mean(null[null>0]) if (null>0).any() else 1)
    else:
        p = (1 + (null <= obs).sum()) / (1 + n_perm)
        nes = obs / (abs(np.mean(null[null<0])) if (null<0).any() else 1)
    return obs, nes, p

summary, dir_rows = [], []
for tis in ["iWAT","gWAT"]:
    c, t = de["cells"], de[tis]
    m = c.merge(t, on="gene", suffixes=("_cell","_tis"))
    print(f"\n===== cells vs {tis}: {len(m)} shared genes =====")

    rho_all = stats.spearmanr(m.stat_cell, m.stat_tis).statistic
    tsig = m[(m.padj_tis < TISSUE_FDR)]
    rho_sig = stats.spearmanr(tsig.stat_cell, tsig.stat_tis).statistic if len(tsig) > 5 else np.nan
    print(f"[RESULT] Spearman(Wald) all={rho_all:.3f} | within {len(tsig)} tissue DEGs={rho_sig:.3f}")

    # directional concordance, borrowing power from the well-powered tissue side
    for lab, sub in [("DOWN(direct?)", tsig[tsig.lfc_tis<0]), ("UP(indirect?)", tsig[tsig.lfc_tis>0])]:
        if len(sub) < 3: continue
        same = int((np.sign(sub.stat_cell)==np.sign(sub.stat_tis)).sum())
        p = stats.binomtest(same, len(sub), 0.5, alternative="greater").pvalue
        frac = same/len(sub)
        print(f"[RESULT] tissue {lab:14s}: {same}/{len(sub)} same sign in cells ({100*frac:.0f}%), binom p={p:.2e}")
        dir_rows.append(dict(comparison=f"cells_vs_{tis}", set=lab, n=len(sub),
                             same_sign=same, frac=frac, binom_p=p))

    # signature GSEA: is the tissue program coordinately shifted in the FULL cell ranking?
    ranked = m.sort_values("stat_cell", ascending=False)
    sbg = dict(zip(ranked.gene, ranked.stat_cell))
    down_set = set(tsig[tsig.lfc_tis<0].gene); up_set = set(tsig[tsig.lfc_tis>0].gene)
    es_d, nes_d, p_d = gsea_es(ranked.gene, down_set, sbg)
    es_u, nes_u, p_u = gsea_es(ranked.gene, up_set, sbg)
    print(f"[RESULT] signature-GSEA in cells: tissue_DOWN NES={nes_d:.2f} p={p_d:.2e} (n={len(down_set)}) | "
          f"tissue_UP NES={nes_u:.2f} p={p_u:.2e} (n={len(up_set)})")

    summary.append(dict(comparison=f"cells_vs_{tis}", n_shared=len(m),
                        rho_all=rho_all, rho_tissueDEGs=rho_sig, n_tissueDEGs=len(tsig),
                        NES_tissueDOWN=nes_d, p_tissueDOWN=p_d,
                        NES_tissueUP=nes_u, p_tissueUP=p_u))
    m.to_csv(f"{OUT}/concordance_genes_{tis}.csv", index=False)

    # ---- quadrant scatter ----
    labset = set(MARK_METAB+MARK_IMM+ANCHOR)
    def quad(r):
        if r.lfc_cell<0 and r.lfc_tis<0: return "concordant down (cell-autonomous)"
        if r.lfc_cell>0 and r.lfc_tis>0: return "concordant up"
        if abs(r.lfc_cell)<0.2 and abs(r.lfc_tis)>0.5: return "tissue-only (non-cell-autonomous)"
        return "other / discordant"
    m["quad"] = m.apply(quad, axis=1)
    colors = {"concordant down (cell-autonomous)":"#2166AC","concordant up":"#B2182B",
              "tissue-only (non-cell-autonomous)":"#F1A340","other / discordant":"#CCCCCC"}
    fig, ax = plt.subplots(figsize=(7.5,7))
    for q,col in colors.items():
        sub=m[m["quad"]==q]; ax.scatter(sub.lfc_tis, sub.lfc_cell, s=6, c=col, alpha=0.45, label=q, linewidths=0)
    lim = np.nanpercentile(np.abs(np.r_[m.lfc_tis, m.lfc_cell]), 99)
    ax.axhline(0,color="grey",lw=.5); ax.axvline(0,color="grey",lw=.5)
    ax.plot([-lim,lim],[-lim,lim],"--",color="grey",lw=.7)
    lab = m[m.gene.isin(labset)]
    ax.scatter(lab.lfc_tis, lab.lfc_cell, s=22, c="black", zorder=5)
    for _,r in lab.iterrows():
        ax.annotate(r.gene, (r.lfc_tis, r.lfc_cell), fontsize=7,
                    xytext=(3,3), textcoords="offset points")
    ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    ax.set_xlabel(f"{tis} log2FC (KO vs CTL)"); ax.set_ylabel("3T3-L1 cells log2FC (KAT8KD vs CTL)")
    ax.set_title(f"KAT8-KD concordance: cells vs {tis}  (Spearman={rho_all:.2f})\n"
                 "down-left = direct/cell-autonomous candidates", fontsize=10)
    ax.legend(fontsize=6, loc="upper left", framealpha=.9)
    fig.tight_layout(); fig.savefig(f"{OUT}/concordance_scatter_{tis}.png", dpi=200); plt.close(fig)

pd.DataFrame(summary).to_csv(f"{OUT}/concordance_SUMMARY.csv", index=False)
pd.DataFrame(dir_rows).to_csv(f"{OUT}/concordance_directional_tests.csv", index=False)
print("\n===== SUMMARY =====")
print(pd.DataFrame(summary).to_string(index=False))
print("\n[DONE] outputs in", OUT)
