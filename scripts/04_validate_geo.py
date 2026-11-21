"""TCGA-LIHC <-> GSE14520 cross-cohort validation.

1. parse GSE14520 (GPL3921) series matrix + GPL3921 annotation (probe -> symbol)
2. differential expression (tumour vs non-tumour) with the same EB model
3. concordance of log2FC with TCGA-LIHC (all genes, robust genes, top genes)
4. survival validation of TCGA-derived prognostic candidates in GSE14520
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# repository root: derived from this file so the scripts work from any checkout
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "geo"
PROC = ROOT / "data" / "processed"
OUT = ROOT / "results" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT))

from stade.geo import load_gse14520  # noqa: E402
from stade.stats import fit_ebayes, significance, bh_adjust  # noqa: E402
from stade.plotting import setup_style, plot_validation, plot_forest, save  # noqa: E402
from stade.survival import cox_univariate  # noqa: E402


def load_annotation(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        lines = [line.rstrip("\n") for line in fh if not line.startswith(("#", "!"))]
    # the file starts with a "^Annotation" marker line; drop non-tabular lines
    lines = [line for line in lines if "\t" in line]
    if not lines:
        raise ValueError("annotation file has no tabular content")
    head = [h.strip().strip('"') for h in lines[0].split("\t")]
    rows = [[v.strip().strip('"') for v in line.split("\t")] for line in lines[1:] if line.strip()]
    rows = [r for r in rows if len(r) == len(head)]
    ann = pd.DataFrame(rows, columns=head)
    return ann


def main() -> int:
    expr, sm, surv = load_gse14520(
        RAW / "GSE14520-GPL3921_series_matrix.txt.gz",
        RAW / "GSE14520_Extra_Supplement.txt.gz",
    )
    print("GSE14520 matrix:", expr.shape, sm["group"].value_counts().to_dict())
    ann = load_annotation(RAW / "GPL3921.annot.gz")
    print("annotation:", ann.shape, ann.columns[:6].tolist())
    col = "Gene symbol"
    ann = ann[[ann.columns[0], col]].rename(columns={ann.columns[0]: "probe", col: "symbol"})
    ann["symbol"] = ann["symbol"].replace({"": np.nan, "---": np.nan})
    ann = ann.dropna().drop_duplicates("probe")

    # collapse probes -> genes (mean of probes, max-mean heuristic is common; use mean)
    m = expr.join(ann.set_index("probe"), how="inner")
    m = m.drop(columns=["symbol"]).groupby(m["symbol"]).mean(numeric_only=True)
    print("gene-level matrix:", m.shape)
    m.to_csv(OUT / "gse14520_gene_matrix.tsv.gz", sep="\t", compression="gzip")
    sm.to_csv(OUT / "gse14520_sample_meta.tsv", sep="\t")
    surv.to_csv(OUT / "gse14520_survival.tsv", sep="\t")

    groups = sm["group"].reindex(m.columns)
    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=m.columns)
    res = fit_ebayes(m, design, "group", trend=True, robust=True)
    de = significance(res.table, alpha=0.05, min_abs_logfc=1.0)
    de["gene"] = de.index
    de.to_csv(OUT / "gse14520_de.tsv", sep="\t")
    print(
        f"GSE14520 DE: d0={res.d0:.2f}  significant={int(de['significant'].sum())} "
        f"(up={int(((de.significant)&(de.logFC>0)).sum())}, down={int(((de.significant)&(de.logFC<0)).sum())})"
    )

    # ---- concordance with TCGA ------------------------------------------
    tcga = pd.read_csv(ROOT / "results" / "de_stability.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(PROC / "gene_meta.tsv", sep="\t")
    if "gene_id" not in meta.columns:
        meta = meta.reset_index().rename(columns={meta.index.name or "index": "gene_id"})
    sym = meta.dropna(subset=["gene_symbol"]).drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]
    tcga["symbol"] = sym.reindex(tcga.index).to_numpy()
    tcga_s = tcga.dropna(subset=["symbol"]).drop_duplicates("symbol").set_index("symbol")
    common = tcga_s.index.intersection(de.index)
    print(f"genes in common: {len(common)}")
    cmp = pd.DataFrame(
        {
            "tcga_logFC": tcga_s.loc[common, "logFC"],
            "geo_logFC": de.loc[common, "logFC"],
        }
    )
    cmp["tcga_sig"] = tcga_s.loc[common, "significant"]
    cmp["tcga_robust"] = tcga_s.loc[common, "robust"]
    cmp["geo_sig"] = (de.loc[common, "adj.P.Val"] < 0.05) & (de.loc[common, "logFC"].abs() >= 1)
    cmp.to_csv(OUT / "concordance_pairs.tsv", sep="\t")

    rho_all = stats.spearmanr(cmp["tcga_logFC"], cmp["geo_logFC"]).statistic
    rho_sig = stats.spearmanr(cmp.loc[cmp["tcga_sig"], "tcga_logFC"], cmp.loc[cmp["tcga_sig"], "geo_logFC"]).statistic
    same_robust = float((np.sign(cmp.loc[cmp["tcga_robust"], "tcga_logFC"]) == np.sign(cmp.loc[cmp["tcga_robust"], "geo_logFC"])).mean())
    geo_sig_in_robust = float(cmp.loc[cmp["tcga_robust"], "geo_sig"].mean())
    # confusion-style replication table
    tab = pd.crosstab(cmp["tcga_robust"], cmp["geo_sig"], margins=True)
    tab.to_csv(OUT / "replication_table.tsv", sep="\t")
    print(f"Spearman rho all={rho_all:.3f} significant={rho_sig:.3f}; robust sign agreement={same_robust:.3f}; GEO sig among robust={geo_sig_in_robust:.3f}")

    setup_style()
    tcga_plot = pd.DataFrame({"gene": cmp.index, "logFC": cmp["tcga_logFC"]})
    geo_plot = pd.DataFrame({"gene": cmp.index, "logFC": cmp["geo_logFC"]})
    plot_validation(tcga_plot, geo_plot, "gene", OUT / "figs" / "validation_logfc.png")
    plot_validation(tcga_plot.loc[cmp["tcga_robust"].to_numpy()], geo_plot.loc[cmp["tcga_robust"].to_numpy()], "gene", OUT / "figs" / "validation_logfc_robust.png")

    # ---- survival validation in GSE14520 --------------------------------
    st = surv.reindex(m.columns).dropna(subset=["OS", "OS.time"])
    st = st[st["OS.time"] > 0]
    st = st.join(groups)
    st = st[st["group"] == "Tumor"]
    print(f"GSE14520 survival cohort: {len(st)} tumour samples, events={int(st['OS'].sum())}")
    expr_t = m.loc[:, st.index]
    cox = cox_univariate(st, expr_t)
    cox.to_csv(OUT / "geo_cox_univariate.tsv", sep="\t", index=False)
    if not cox.empty:
        print(f"GSE14520 Cox: tested {len(cox)}, padj<0.05: {int((cox['padj']<0.05).sum())}")

        # restrict to TCGA-derived candidate genes
        cand_path = ROOT / "results" / "survival" / "cox_univariate.tsv"
        if cand_path.exists():
            tcga_cox = pd.read_csv(cand_path, sep="\t")
            cand = set(tcga_cox.loc[tcga_cox["padj"] < 0.05, "gene"])
            sub = cox[cox["gene"].isin(cand)]
            sub.to_csv(OUT / "geo_cox_candidates.tsv", sep="\t", index=False)
            both = sub[(sub["padj"] < 0.05)]
            print(f"TCGA prognostic candidates tested in GSE14520: {len(sub)}, replicated (padj<0.05): {len(both)}")
            if len(sub):
                print(f"  sign agreement: {(np.sign(sub['HR']-1)==np.sign(tcga_cox.set_index('gene').loc[sub.gene,'HR'].to_numpy()-1)).mean():.3f}")
            if len(both):
                plot_forest(both, OUT / "figs" / "geo_forest_replicated.png", top_n=15, title="Replicated prognostic genes (GSE14520)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
