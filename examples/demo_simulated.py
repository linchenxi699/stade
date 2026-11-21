"""End-to-end demo on simulated data (runs in < 1 minute, no downloads).

Usage::

    python examples/demo_simulated.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from stade.enrichment import GeneSet, GeneSetCollection, gsea, ora_weighted  # noqa: E402
from stade.normalization import filter_low_expression, log_cpm  # noqa: E402
from stade.plotting import plot_stability_scatter, plot_volcano, save, setup_style  # noqa: E402
from stade.stability import bootstrap_de, consensus_table, gene_stability  # noqa: E402
from stade.stats import fit_ebayes, significance, welch_test, wilcoxon_test  # noqa: E402

OUT = ROOT / "examples" / "output"
OUT.mkdir(parents=True, exist_ok=True)


def simulate(n_genes=3000, n_case=15, n_ctrl=12, n_de=150, seed=42):
    rng = np.random.default_rng(seed)
    lib = rng.uniform(5e6, 9e6, size=n_case + n_ctrl)
    base = rng.gamma(0.8, 80, size=n_genes)
    counts = np.outer(base, lib / 1e6) * rng.gamma(6, 0.18, size=(n_genes, n_case + n_ctrl))
    effect = np.zeros(n_genes)
    effect[:n_de] = rng.choice([-1, 1], n_de) * rng.uniform(1.2, 2.4, n_de)
    is_case = np.r_[np.ones(n_case, bool), np.zeros(n_ctrl, bool)]
    counts[:, is_case] *= 2 ** effect[:, None]
    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:02d}" for i in range(n_case + n_ctrl)]
    df = pd.DataFrame(np.round(counts).astype(int), index=genes, columns=samples)
    groups = pd.Series(np.where(is_case, "Tumor", "Normal"), index=samples, name="group")
    return df, groups, set(genes[:n_de])


def main() -> int:
    t0 = time.time()
    counts, groups, truth = simulate()
    print(f"simulated matrix: {counts.shape}, planted DE genes: {len(truth)}")
    kept = filter_low_expression(counts, min_count=10, min_prop=0.5, groups=groups)
    expr = log_cpm(kept, use_tmm=True)
    print(f"after filtering: {expr.shape}")

    design = pd.DataFrame({"Intercept": 1.0, "group": (groups == "Tumor").astype(float)}, index=expr.columns)
    res = fit_ebayes(expr, design, "group")
    de = significance(res.table, alpha=0.05, min_abs_logfc=1.0)
    print(f"EB: d0={res.d0:.2f}, significant={int(de['significant'].sum())}")

    boot = bootstrap_de(expr, groups, "Tumor", n_boot=60, seed=7)
    stab = gene_stability(boot, de["logFC"])
    de = de.join(stab)

    cons = consensus_table(res.table, welch_test(expr, groups, "Tumor"), wilcoxon_test(expr, groups, "Tumor"))
    cons = cons.join(stab[["stability", "sel_freq"]])
    cons.to_csv(OUT / "demo_de_consensus.tsv", sep="\t")

    tp = len(truth & set(de.index[de["significant"]]))
    fp = int(de["significant"].sum()) - tp
    print(f"recovery: TP={tp}/{len(truth)}, false positives={fp}")
    print(f"median stability: planted={de.loc[list(truth & set(de.index)), 'stability'].median():.3f}, "
          f"others={de.loc[~de.index.isin(truth), 'stability'].median():.3f}")

    # toy gene sets: one true set enriched with the planted genes
    gs = GeneSetCollection(
        {
            "planted": GeneSet("planted", frozenset(truth)),
            "random": GeneSet("random", frozenset(sorted(set(de.index) - truth)[:150])),
        }
    )
    w = de["stability"]
    enr = ora_weighted(w, gs, universe=set(de.index), min_size=20, max_size=500)
    print(enr[["term", "set_size", "weight_ratio", "pvalue", "padj"]].to_string(index=False))
    enr.to_csv(OUT / "demo_enrichment.tsv", sep="\t", index=False)

    setup_style()
    plot_volcano(de, OUT / "demo_volcano.png", label_genes=list(truth)[:10])
    plot_stability_scatter(de, stab, OUT / "demo_stability.png")
    de.to_csv(OUT / "demo_de.tsv", sep="\t")
    print(f"done in {time.time()-t0:.1f}s -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
