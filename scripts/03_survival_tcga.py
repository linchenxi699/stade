"""TCGA-LIHC survival analysis on candidate genes.

* univariate Cox across all expressed genes (z-scored expression per SD)
* Kaplan-Meier + log-rank for the top candidates
* multivariate Cox for a compact signature
* joint DE x prognosis candidate ranking
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "survival"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT))

from stade.io import Provenance, read_matrix  # noqa: E402
from stade.plotting import (  # noqa: E402
    plot_forest,
    plot_km,
    save,
    setup_style,
)
from stade.survival import candidate_signature, cox_multivariate, cox_univariate, km_analysis  # noqa: E402


def main() -> int:
    prov = Provenance(command="survival_tcga_lihc", parameters={"cut": "median"})
    expr = read_matrix(ROOT / "work" / "logcpm.tsv.gz")
    counts = read_matrix(ROOT / "work" / "filtered_counts.tsv.gz")
    meta = read_matrix(ROOT / "data" / "processed" / "sample_meta.tsv", index_col=None)
    clin = read_matrix(ROOT / "data" / "processed" / "clinical.tsv")
    prov.add_input("logcpm", ROOT / "work" / "logcpm.tsv.gz")
    prov.add_input("clinical", ROOT / "data" / "processed" / "clinical.tsv")

    # sample -> case mapping (one RNA-seq sample per tumour case here)
    smeta = meta.set_index("file_id")
    tumor = smeta[smeta["group"] == "Tumor"]
    sample_case = tumor["case_id"]
    clin = clin.reindex(sample_case.unique())
    expr_t = expr.loc[:, [s for s in sample_case.index if s in expr.columns]]
    expr_t.columns = [sample_case[s] for s in expr_t.columns]
    expr_t = expr_t.loc[:, ~expr_t.columns.duplicated()]
    print(f"tumour samples with expression: {expr_t.shape[1]}; clinical cases: {len(clin)}")
    valid = clin.dropna(subset=["OS", "OS.time"])
    valid = valid[valid["OS.time"] > 0].reindex(expr_t.columns).dropna(subset=["OS.time"])
    expr_t = expr_t.loc[:, valid.index]
    print(f"survival cohort: {expr_t.shape[1]} samples, events={int(valid['OS'].sum())}, median follow-up={valid['OS.time'].median()/30.44:.1f} months")

    # ---- univariate Cox on all genes -----------------------------------
    cox = cox_univariate(valid, expr_t)
    # map Ensembl -> symbol
    gm = read_matrix(ROOT / "data" / "processed" / "gene_meta.tsv", index_col=None)
    sym = gm.dropna(subset=["gene_name"]).drop_duplicates("gene_id").set_index("gene_id")["gene_name"]
    cox["symbol"] = sym.reindex(cox["gene"]).to_numpy()
    cox = cox.drop_duplicates("symbol")
    cox.to_csv(OUT / "cox_univariate.tsv", sep="\t", index=False)
    print(f"Cox tested {len(cox)} genes; padj<0.05: {int((cox['padj']<0.05).sum())}")
    print(cox.head(12)[["symbol", "HR", "pvalue", "padj"]].to_string(index=False))

    # ---- KM for the top candidates -------------------------------------
    setup_style()
    top = cox.dropna(subset=["symbol"]).sort_values("pvalue").head(6)
    km_rows = []
    km_curves = {}
    for _, row in top.iterrows():
        gene = row["gene"]
        score = expr_t.loc[gene]
        res = km_analysis(valid["OS.time"], valid["OS"], score, cut="median")
        km_rows.append(
            {
                "gene": row["symbol"],
                "ensembl": gene,
                "HR": row["HR"],
                "cox_p": row["pvalue"],
                "logrank_p": res.logrank_p,
                "median_high_months": res.median_high,
                "median_low_months": res.median_low,
                "n_high": res.n_high,
                "n_low": res.n_low,
            }
        )
        km_curves[row["symbol"]] = (
            valid["OS.time"],
            valid["OS"],
            score,
        )
    pd.DataFrame(km_rows).to_csv(OUT / "km_summary.tsv", sep="\t", index=False)
    # individual KM panels
    from lifelines import KaplanMeierFitter
    import matplotlib.pyplot as plt

    for name, (t, e, score) in km_curves.items():
        cut = score.median()
        hi = score > cut
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        f = KaplanMeierFitter()
        f.fit(t[hi], e[hi], label=f"{name} high")
        f.plot_survival_function(ax=ax, ci_show=True, lw=1.1, color="#C0392B", show_censors=True, censor_styles={"ms": 2})
        f2 = KaplanMeierFitter()
        f2.fit(t[~hi], e[~hi], label=f"{name} low")
        f2.plot_survival_function(ax=ax, ci_show=True, lw=1.1, color="#2E86C1", show_censors=True, censor_styles={"ms": 2})
        p = km_analysis(t, e, score, cut="median").logrank_p
        ax.set_title(f"{name}  log-rank p={p:.2e}", fontsize=8)
        ax.set_xlabel("Time (months)")
        ax.set_ylabel("Survival probability")
        save(fig, OUT / "figs" / f"km_{name}.png")

    plot_forest(cox.rename(columns={"symbol": "gene"}), OUT / "figs" / "cox_forest.png", top_n=18)

    # ---- multivariate signature ----------------------------------------
    panel = [g for g in top["gene"].tolist() if g in expr_t.index]
    msum, cph = cox_multivariate(valid, expr_t, panel)
    msum_out = msum.copy()
    msum_out.index = [sym.get(g, g) for g in msum_out.index]
    msum_out.to_csv(OUT / "cox_multivariate.tsv", sep="\t")
    print("multivariate params:")
    print(msum_out[["coef", "exp(coef)", "p"]].to_string())

    # ---- joint DE x survival ranking -----------------------------------
    de = pd.read_csv(ROOT / "results" / "de_stability.tsv", sep="\t", index_col=0)
    de["symbol"] = sym.reindex(de.index).to_numpy()
    de_s = de.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    cox_s = cox.dropna(subset=["symbol"]).set_index("symbol")
    cox_s["gene"] = cox_s.index
    cand = candidate_signature(de_s, cox_s.reset_index(), top_n=25, require_stability=0.7)
    cand.to_csv(OUT / "candidates_joint.tsv", sep="\t")
    print("\njoint DE + prognosis candidates:")
    print(cand[["logFC", "adj.P.Val", "stability", "HR", "cox_p"]].head(15).to_string())

    summary = {
        "n_samples": int(expr_t.shape[1]),
        "events": int(valid["OS"].sum()),
        "median_followup_months": float(valid["OS.time"].median() / 30.44),
        "n_genes_tested": int(len(cox)),
        "n_prognostic": int((cox["padj"] < 0.05).sum()),
        "top_prognostic": cox.head(10)[["symbol", "HR", "pvalue", "padj"]].to_dict(orient="records"),
        "multivariate": msum_out[["exp(coef)", "p"]].to_dict(orient="index"),
        "n_joint_candidates": int(len(cand)),
    }
    (OUT / "survival_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    prov.finish(**{k: v for k, v in summary.items() if isinstance(v, (int, float))})
    prov.write(OUT / "provenance_survival.json")
    prov.add_input("expr", ROOT / "work" / "logcpm.tsv.gz")
    print("\nsummary written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
