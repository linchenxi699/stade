"""Publication-oriented figures (matplotlib, no seaborn dependency).

All functions take an explicit ``out`` path and return it, so the analysis
driver can record every figure in the provenance record.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, linkage  # noqa: E402

__all__ = [
    "setup_style",
    "save",
    "plot_pca",
    "plot_library_density",
    "plot_ma",
    "plot_volcano",
    "plot_stability_scatter",
    "plot_heatmap",
    "plot_gene_boxplots",
    "plot_enrichment_dot",
    "plot_km",
    "plot_forest",
    "plot_method_concordance",
    "plot_bootstrap_fdr",
    "plot_validation",
    "plot_enrichment_running_sum",
]

STADE_PALETTE = {
    "up": "#C0392B",
    "down": "#2E86C1",
    "neutral": "#B3B6B7",
    "tumor": "#D35400",
    "normal": "#1F618D",
    "accent": "#117A65",
    "warn": "#B7950B",
    "grid": "#E5E7E9",
}


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Microsoft YaHei"],
            "axes.unicode_minus": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "axes.grid": False,
            "figure.constrained_layout.use": True,
        }
    )


def save(fig, out: str | Path, close: bool = True) -> Path:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    if close:
        plt.close(fig)
    return out


def _pca_basis(expr: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    x = expr.T.to_numpy(dtype=float)
    x = x - np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0)
    sd[sd == 0] = 1.0
    x = x / sd
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    var = s**2 / np.sum(s**2)
    return u[:, :2] * s[:2], var[:2]


def plot_pca(expr: pd.DataFrame, groups: pd.Series, out: str | Path, title: str = "PCA of log-CPM") -> Path:
    scores, var = _pca_basis(expr)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    for i, g in enumerate(pd.unique(groups.reindex(expr.columns))):
        m = (groups.reindex(expr.columns) == g).to_numpy()
        ax.scatter(scores[m, 0], scores[m, 1], s=14, alpha=0.75, label=str(g), edgecolor="none")
    ax.set_xlabel(f"PC1 ({var[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.1f}% var)")
    ax.set_title(title)
    ax.legend(loc="best")
    return save(fig, out)


def plot_library_density(log_expr: pd.DataFrame, out: str | Path) -> Path:
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    for i in range(0, log_expr.shape[1], max(1, log_expr.shape[1] // 40)):
        vals = log_expr.iloc[:, i].dropna()
        ax.plot(np.sort(vals), np.linspace(0, 1, vals.size), lw=0.6, alpha=0.35, color="#566573")
    ax.set_xlabel("log2 CPM")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Expression distributions after normalisation")
    return save(fig, out)


def plot_ma(de: pd.DataFrame, out: str | Path, alpha: float = 0.05, min_abs_logfc: float = 1.0) -> Path:
    sig = (de["adj.P.Val"] < alpha) & (de["logFC"].abs() >= min_abs_logfc)
    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.scatter(de.loc[~sig, "AveExpr"], de.loc[~sig, "logFC"], s=2, c=STADE_PALETTE["neutral"], alpha=0.4, rasterized=True)
    up = sig & (de["logFC"] > 0)
    dn = sig & (de["logFC"] < 0)
    ax.scatter(de.loc[up, "AveExpr"], de.loc[up, "logFC"], s=3, c=STADE_PALETTE["up"], alpha=0.7, rasterized=True, label=f"Up ({int(up.sum())})")
    ax.scatter(de.loc[dn, "AveExpr"], de.loc[dn, "logFC"], s=3, c=STADE_PALETTE["down"], alpha=0.7, rasterized=True, label=f"Down ({int(dn.sum())})")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("Average log2 CPM")
    ax.set_ylabel("log2 fold change")
    ax.set_title("Mean-difference plot")
    ax.legend(loc="upper right", markerscale=4)
    return save(fig, out)


def plot_volcano(
    de: pd.DataFrame,
    out: str | Path,
    alpha: float = 0.05,
    min_abs_logfc: float = 1.0,
    label_genes: list[str] | None = None,
) -> Path:
    y = -np.log10(np.clip(de["adj.P.Val"], 1e-300, 1))
    fig, ax = plt.subplots(figsize=(4.8, 3.9))
    sig = (de["adj.P.Val"] < alpha) & (de["logFC"].abs() >= min_abs_logfc)
    ax.scatter(de.loc[~sig, "logFC"], y[~sig], s=2, c=STADE_PALETTE["neutral"], alpha=0.45, rasterized=True)
    ax.scatter(de.loc[sig & (de["logFC"] > 0), "logFC"], y[sig & (de["logFC"] > 0)], s=3, c=STADE_PALETTE["up"], alpha=0.7, rasterized=True)
    ax.scatter(de.loc[sig & (de["logFC"] < 0), "logFC"], y[sig & (de["logFC"] < 0)], s=3, c=STADE_PALETTE["down"], alpha=0.7, rasterized=True)
    if label_genes:
        for g in label_genes:
            if g in de.index:
                ax.annotate(g, (de.at[g, "logFC"], -np.log10(max(de.at[g, "adj.P.Val"], 1e-300))), fontsize=6, xytext=(2, 2), textcoords="offset points")
    ax.axhline(-np.log10(alpha), ls="--", lw=0.7, c="k")
    ax.axvline(min_abs_logfc, ls="--", lw=0.7, c="k")
    ax.axvline(-min_abs_logfc, ls="--", lw=0.7, c="k")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("Differential expression (tumour vs normal)")
    return save(fig, out)


def plot_stability_scatter(de: pd.DataFrame, stab: pd.DataFrame, out: str | Path) -> Path:
    extra = [c for c in ("stability", "sel_freq") if c not in de.columns]
    d = de.join(stab[extra], how="inner") if extra else de.copy()
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
    ax = axes[0]
    ax.scatter(d["t"].abs(), d["stability"], s=2, c=STADE_PALETTE["accent"], alpha=0.35, rasterized=True)
    ax.set_xlabel("|moderated t|")
    ax.set_ylabel("stability score")
    ax.set_title("Stability vs. effect size")
    ax = axes[1]
    hist = np.histogram(d["stability"].dropna(), bins=40)
    ax.bar(hist[1][:-1], hist[0], width=np.diff(hist[1])[0] * 0.95, color="#5D6D7E")
    ax.set_xlabel("stability score")
    ax.set_ylabel("genes")
    ax.set_title("Distribution of stability")
    return save(fig, out)


def plot_heatmap(
    expr: pd.DataFrame,
    genes: list[str],
    groups: pd.Series,
    out: str | Path,
    zscore: bool = True,
    figsize: tuple[float, float] = (6.2, 5.0),
) -> Path:
    genes = [g for g in genes if g in expr.index]
    sub = expr.loc[genes]
    if zscore:
        mu = sub.mean(axis=1)
        sd = sub.std(axis=1).replace(0, 1)
        sub = sub.sub(mu, axis=0).div(sd, axis=0)
    order = dendrogram(linkage(sub.T.to_numpy(), method="average", metric="correlation"), no_plot=True)["leaves"]
    cols = sub.columns[order]
    order_r = dendrogram(linkage(sub.to_numpy(), method="average", metric="correlation"), no_plot=True)["leaves"]
    rows = sub.index[order_r]
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(sub.loc[rows, cols].to_numpy(), aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-2.5, vmax=2.5), interpolation="nearest")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=5)
    gs = groups.reindex(cols).astype(str)
    ax.set_xticks([])
    for i, g in enumerate(pd.unique(gs)):
        idx = np.flatnonzero(gs.to_numpy() == g)
        ax.add_patch(plt.Rectangle((idx.min() - 0.5, len(rows) - 0.2), np.ptp(idx) + 1, 1.0, clip_on=False, color=STADE_PALETTE["tumor"] if i else STADE_PALETTE["normal"], alpha=0.8))
    ax.set_title("Top differential genes (row z-score)")
    fig.colorbar(im, ax=ax, shrink=0.6, label="z-score")
    return save(fig, out)


def plot_gene_boxplots(
    expr: pd.DataFrame,
    groups: pd.Series,
    genes: list[str],
    out: str | Path,
    ncols: int = 4,
) -> Path:
    genes = [g for g in genes if g in expr.index]
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 1.9, nrows * 2.0), squeeze=False)
    gl = groups.reindex(expr.columns)
    order = list(pd.unique(gl.dropna()))
    for k, g in enumerate(genes):
        ax = axes[k // ncols][k % ncols]
        data = [expr.loc[g, gl == gg].dropna() for gg in order]
        bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False)
        for i, box in enumerate(bp["boxes"]):
            box.set_facecolor(STADE_PALETTE["tumor"] if i == 0 else STADE_PALETTE["normal"])
            box.set_alpha(0.75)
        ax.set_xticks(range(1, len(order) + 1), order, fontsize=7)
        ax.set_title(g, fontsize=8)
        ax.tick_params(axis="y", labelsize=6)
    for k in range(len(genes), nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")
    return save(fig, out)


def plot_enrichment_dot(
    df: pd.DataFrame,
    out: str | Path,
    top_n: int = 15,
    score: str = "padj",
    size: str = "overlap",
    x: str | None = None,
    color_by: str | None = None,
) -> Path:
    """Dot plot of enrichment results.

    ``x`` chooses the horizontal axis (default: -log10 of ``score``); ``color_by``
    maps the colour scale to another column (default: same as the x axis).
    """
    d = df.head(top_n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.4, max(2.6, 0.28 * len(d) + 0.8)))
    if x is not None and x in d:
        xv = d[x].to_numpy(dtype=float)
        xlabel = {
            "fold_enrichment": "Fold enrichment (observed / expected)",
            "evidence": "Integrated evidence score",
        }.get(x, x)
    elif score in d:
        xv = -np.log10(np.clip(d[score].to_numpy(dtype=float), 1e-300, 1))
        xlabel = f"-log10({score})"
    else:
        xv = d.get("evidence", pd.Series(np.ones(len(d)))).to_numpy()
        xlabel = "evidence"
    csrc = color_by if (color_by and color_by in d) else score
    if csrc in d:
        cv = -np.log10(np.clip(d[csrc].to_numpy(dtype=float), 1e-300, 1))
        clabel = f"-log10(padj)"
    else:
        cv = xv
        clabel = xlabel
    sizes = d[size].to_numpy() if size in d else np.full(len(d), 10)
    span = np.ptp(sizes) if sizes.size else 0
    sizes = 12 + 130 * (sizes - sizes.min()) / (span if span > 0 else 1)
    sc = ax.scatter(xv, np.arange(len(d)), s=sizes, c=cv, cmap="viridis", edgecolor="k", linewidth=0.3)
    ax.set_yticks(np.arange(len(d)), [str(t)[:42] for t in d["term"]], fontsize=6.5)
    ax.set_xlabel(xlabel)
    ax.set_title("Functional enrichment (top terms)")
    fig.colorbar(sc, ax=ax, shrink=0.5, label=clabel)
    return save(fig, out)


def plot_km(km, out: str | Path, title: str = "Overall survival by expression") -> Path:
    """Kaplan-Meier curves from a fitted ``lifelines`` KaplanMeierFitter dict.

    ``km`` maps group label -> ``(time, event)`` arrays.
    """
    from lifelines import KaplanMeierFitter

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    colors = [STADE_PALETTE["tumor"], STADE_PALETTE["normal"], STADE_PALETTE["accent"]]
    for i, (label, (t, e)) in enumerate(km.items()):
        f = KaplanMeierFitter()
        f.fit(np.asarray(t, float), np.asarray(e, float), label=str(label))
        f.plot_survival_function(ax=ax, ci_show=True, color=colors[i % len(colors)], lw=1.1, show_censors=True, censor_styles={"ms": 2})
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.set_title(title)
    return save(fig, out)


def plot_forest(coef: pd.DataFrame, out: str | Path, top_n: int = 15, title: str = "Univariate Cox (HR per SD)") -> Path:
    d = coef.sort_values("pvalue").head(top_n).iloc[::-1]
    label_col = "gene" if "gene" in d.columns else ("symbol" if "symbol" in d.columns else d.columns[0])
    fig, ax = plt.subplots(figsize=(5.2, max(2.4, 0.3 * len(d) + 0.8)))
    y = np.arange(len(d))
    ax.errorbar(d["HR"], y, xerr=[d["HR"] - d["HR_lower"], d["HR_upper"] - d["HR"]], fmt="o", ms=3, color="#1F618D", ecolor="#7FB3D5", lw=0.9)
    ax.axvline(1, ls="--", lw=0.7, c="k")
    ax.set_yticks(y, [str(v) for v in d[label_col]], fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("Hazard ratio (95% CI, log scale)")
    ax.set_title(title)
    return save(fig, out)


def plot_method_concordance(cons: pd.DataFrame, out: str | Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))
    ax = axes[0]
    methods = ["adj.P.Val", "adj.welch", "adj.wilcoxon"]
    counts = [(cons[m] < 0.05).sum() for m in methods]
    ax.bar(range(3), counts, color=["#117A65", "#B7950B", "#7D3C98"], width=0.6)
    ax.set_xticks(range(3), ["limma-EB", "Welch", "Wilcoxon"], fontsize=8)
    ax.set_ylabel("significant genes")
    ax.set_title("Significance calls by method")
    ax = axes[1]
    ax.hist(cons["n_methods_sig"], bins=[-0.5, 0.5, 1.5, 2.5, 3.5], color="#566573", rwidth=0.8)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xlabel("methods agreeing")
    ax.set_ylabel("genes")
    ax.set_title("Method agreement")
    return save(fig, out)


def plot_bootstrap_fdr(stab: pd.DataFrame, de: pd.DataFrame, out: str | Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.3))
    ax = axes[0]
    sig = de["adj.P.Val"] < 0.05
    ax.hist(stab.loc[sig, "sel_freq"], bins=25, color=STADE_PALETTE["up"], alpha=0.8, label="significant")
    ax.hist(stab.loc[~sig, "sel_freq"], bins=25, color=STADE_PALETTE["neutral"], alpha=0.7, label="not significant")
    ax.set_xlabel("bootstrap selection frequency")
    ax.set_ylabel("genes")
    ax.legend()
    ax.set_title("Reproducibility of DE calls")
    ax = axes[1]
    bins = np.linspace(0, 1, 21)
    idx = np.digitize(stab["stability"].to_numpy(), bins) - 1
    frac = np.array([(de["adj.P.Val"].to_numpy()[idx == b] < 0.05).mean() if np.any(idx == b) else np.nan for b in range(20)])
    ax.plot(0.5 * (bins[:-1] + bins[1:]), frac, "o-", ms=3, color=STADE_PALETTE["accent"])
    ax.set_xlabel("stability score")
    ax.set_ylabel("fraction significant (FDR<0.05)")
    ax.set_title("Stability predicts significance")
    return save(fig, out)


def plot_validation(
    tcga: pd.DataFrame,
    geo: pd.DataFrame,
    gene_col: str,
    out: str | Path,
    title: str = "TCGA-LIHC vs GSE14520 log2FC",
) -> Path:
    m = tcga.set_index("gene")["logFC"].to_frame("tcga").join(geo.set_index(gene_col)["logFC"].to_frame("geo"), how="inner")
    if m.empty:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.text(0.5, 0.5, "no overlapping genes", ha="center")
        return save(fig, out)
    r = float(m["tcga"].corr(m["geo"], method="spearman"))
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    ax.scatter(m["tcga"], m["geo"], s=3, alpha=0.35, c="#1F618D", rasterized=True)
    lim = np.nanpercentile(np.abs(np.r_[m["tcga"], m["geo"]]), 99.5)
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=0.8, c="k")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("TCGA-LIHC log2FC")
    ax.set_ylabel("GSE14520 log2FC")
    ax.set_title(f"{title}\nSpearman rho = {r:.3f} (n={len(m)})")
    return save(fig, out)


def plot_enrichment_running_sum(
    terms: list[dict],
    out: str | Path,
    title: str = "GSEA running enrichment score",
) -> Path:
    fig, axes = plt.subplots(len(terms), 1, figsize=(6.0, 1.6 * len(terms)), squeeze=False)
    for ax, t in zip(axes[:, 0], terms):
        run = np.asarray(t["run"])
        ax.plot(run, lw=1.0, color=t.get("color", "#117A65"))
        ax.axhline(0, lw=0.5, c="k")
        ax.set_title(t["term"], fontsize=8)
        ax.set_ylabel("ES", fontsize=7)
        if t.get("es") is not None:
            ax.annotate(f"ES={t['es']:.2f}\np={t['p']:.2e}", (0.02, 0.05), xycoords="axes fraction", fontsize=6)
    axes[-1, 0].set_xlabel("ranked gene list position")
    fig.suptitle(title, fontsize=9)
    return save(fig, out)
