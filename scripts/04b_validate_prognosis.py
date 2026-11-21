"""Prognostic cross-cohort validation: TCGA candidates re-tested in GSE14520."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / "results" / "validation"
OUT = ROOT / "results" / "survival"
sys.path.insert(0, str(ROOT))

from stade.io import read_matrix  # noqa: E402
from stade.plotting import plot_forest, setup_style  # noqa: E402
from stade.survival import cox_univariate  # noqa: E402


def main() -> int:
    # TCGA prognostic genes with symbols
    cox = pd.read_csv(OUT / "cox_univariate.tsv", sep="\t")
    gm = read_matrix(ROOT / "data" / "processed" / "gene_meta.tsv", index_col=None)
    sym = gm.dropna(subset=["gene_name"]).drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
    cox["symbol"] = cox["symbol"].fillna(pd.Series(sym.reindex(cox["gene"]).to_numpy(), index=cox.index))
    sig = cox[(cox["padj"] < 0.05) & cox["symbol"].notna()].drop_duplicates("symbol")
    print(f"TCGA prognostic genes (FDR<0.05 with symbol): {len(sig)}")

    # GSE14520 cohort
    expr = read_matrix(VAL / "gse14520_gene_matrix.tsv.gz")
    surv = read_matrix(VAL / "gse14520_survival.tsv")
    meta = read_matrix(VAL / "gse14520_sample_meta.tsv")
    tumor = meta.index[meta["group"] == "Tumor"]
    st = surv.reindex(tumor).dropna(subset=["OS", "OS.time"])
    st = st[st["OS.time"] > 0]
    expr_t = expr.loc[:, [s for s in st.index if s in expr.columns]]
    st = st.loc[expr_t.columns]
    print(f"GSE14520 survival cohort: {expr_t.shape[1]} tumour samples, events={int(st['OS'].sum())}")

    panel = [g for g in sig["symbol"] if g in expr_t.index]
    print(f"candidates present in GSE14520: {len(panel)}")
    sub_expr = expr_t.loc[[g for g in dict.fromkeys(panel)]]
    geo = cox_univariate(st, sub_expr)
    geo["symbol"] = geo["gene"]
    tcga_map = sig.set_index("symbol")["HR"]
    geo["tcga_HR"] = tcga_map.reindex(geo["symbol"]).to_numpy()
    geo["same_direction"] = np.sign(geo["HR"] - 1) == np.sign(geo["tcga_HR"] - 1)
    geo.to_csv(VAL / "geo_cox_candidates.tsv", sep="\t", index=False)
    rep = geo[geo["padj"] < 0.05]
    print(f"tested {len(geo)}; replicated (FDR<0.05): {len(rep)}; direction concordance: {geo['same_direction'].mean():.3f}")
    if len(rep):
        setup_style()
        plot_forest(rep, VAL / "figs" / "geo_forest_replicated.png", top_n=15, title="Replicated prognostic genes (GSE14520)")
    summary = {
        "n_tcga_prognostic": int(len(sig)),
        "n_tested_in_geo": int(len(geo)),
        "n_replicated": int(len(rep)),
        "direction_concordance": float(geo["same_direction"].mean()),
    }
    (VAL / "prognostic_validation.json").write_text(pd.Series(summary).to_json(indent=2), encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
