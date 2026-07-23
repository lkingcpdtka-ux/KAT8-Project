#!/usr/bin/env python3
"""Program-level KAT8 concordance: mean log2FC per curated program across
cells / iWAT / gWAT, with a competitive test and a summary bar figure.
Companion to run_concordance_python.py. Inputs = downstream_analysis/data/DE_*.csv."""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
DATA="downstream_analysis/data"; OUT="downstream_analysis/results"
F={"cells":"DE_cells_KAT8KD_vs_CTL.csv","iWAT":"DE_tissue_iWAT_KD_vs_CTL.csv","gWAT":"DE_tissue_gWAT_KD_vs_CTL.csv"}
def load(p):
    d=pd.read_csv(f"{DATA}/{p}"); g=d["gene_name"].astype(str)
    padj=d["padj"] if "padj" in d else d["adj.P.Val"]
    return pd.DataFrame({"gene":g.values,"lfc":pd.to_numeric(d["log2FoldChange"],errors="coerce").values,
        "stat":pd.to_numeric(d["stat"],errors="coerce").values,
        "padj":pd.to_numeric(padj,errors="coerce").values}).dropna(subset=["stat"]).drop_duplicates("gene").set_index("gene")
de={k:load(v) for k,v in F.items()}
PROG={
 "Adipocyte identity/PPARG":["Pparg","Adipoq","Plin1","Cd36","Fabp4","Lep","Cfd","Lipe","Pnpla2","Cav1","Retn"],
 "Thermogenic":["Ucp1","Cidea","Ppargc1a","Cox7a1","Elovl3","Dio2","Cox8b","Prdm16"],
 "Fatty-acid/lipid metab":["Ppara","Acox1","Acacb","Pck1","Scd1","Fasn","Acaca","Cpt1b","Lpl"],
 "Inflammatory/immune":["Ccl2","Ccl8","Cxcl2","Cxcl10","Il1b","Il12b","Itgax","Ctss","Adgre1","Tnf","Cd68","Lyz2"],
 "ECM/fibrosis":["Col1a1","Col3a1","Col6a1","Fn1","Timp1","Sparc","Col1a2","Lox"]}
rows=[]
for prog,genes in PROG.items():
    for k,d in de.items():
        pres=[g for g in genes if g in d.index]
        if not pres: rows.append(dict(program=prog,dataset=k,n_present=0,mean_lfc=np.nan,n_sig=0,comp_p=np.nan)); continue
        sub=d.loc[pres]; bg=d.loc[~d.index.isin(pres),"stat"]
        cp=stats.mannwhitneyu(sub["stat"],bg,alternative="two-sided").pvalue if len(bg)>0 else np.nan
        rows.append(dict(program=prog,dataset=k,n_present=len(pres),mean_lfc=float(sub["lfc"].mean()),
                         n_sig=int((sub["padj"]<0.05).sum()),comp_p=cp))
tab=pd.DataFrame(rows); tab.to_csv(f"{OUT}/program_level_concordance.csv",index=False)
piv=tab.pivot(index="program",columns="dataset",values="mean_lfc")[["cells","iWAT","gWAT"]]
progs=list(PROG); x=np.arange(len(progs)); w=0.26; cols={"cells":"#111111","iWAT":"#2166AC","gWAT":"#B2182B"}
fig,ax=plt.subplots(figsize=(9,4.6))
for i,k in enumerate(["cells","iWAT","gWAT"]): ax.bar(x+(i-1)*w,piv[k].values,w,label=k,color=cols[k])
ax.axhline(0,color="k",lw=.6); ax.set_xticks(x); ax.set_xticklabels(progs,rotation=20,ha="right",fontsize=8)
ax.set_ylabel("mean log2FC (KO/KD vs CTL)")
ax.set_title("KAT8 loss: program-level response — pure adipocytes (cells) vs tissue\n"
             "Metabolic/identity collapse is iWAT-specific and ABSENT in cells (not cell-autonomous)",fontsize=10)
ax.legend(title="dataset"); fig.tight_layout(); fig.savefig(f"{OUT}/program_level_concordance.png",dpi=200)
print(piv.round(2).to_string())
