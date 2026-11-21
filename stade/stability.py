"""Stability, uncertainty and multiple-testing calibration utilities.

Three complementary resampling procedures are provided:

``bootstrap_de``
    Non-parametric bootstrap over samples (stratified by group). For every
    replicate the DE model is refitted, yielding the empirical distribution of
    log fold-changes and adjusted p-values under sample perturbation.

``gene_stability``
    Turns those distributions into three per-gene stability components -
    selection frequency, sign consistency and effect-size consistency - and a
    composite stability score in [0, 1].

``permutation_de`` / ``empirical_null``
    Group-label permutation for null calibration and an Efron-style empirical
    null giving local FDRs as an orthogonal, model-free FDR estimate.

``consensus_table``
    Merges the empirical-Bayes result with Welch and Wilcoxon tests and adds
    per-gene method concordance.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

from .stats import bh_adjust, fit_ebayes

__all__ = [
    "BootstrapResult",
    "bootstrap_de",
    "gene_stability",
    "permutation_de",
    "empirical_null",
    "consensus_table",
]


@dataclass
class BootstrapResult:
    logfc: pd.DataFrame  # genes x replicates
    padj: pd.DataFrame
    n_ok: int
    seeds: list[int]
    tstat: pd.DataFrame | None = None  # genes x replicates, useful for ranking


def bootstrap_de(
    log_expr: pd.DataFrame,
    groups: pd.Series,
    group_a: str,
    n_boot: int = 200,
    seed: int = 20261215,
    stratified: bool = True,
    method: str = "ebayes",
) -> BootstrapResult:
    """Bootstrap the DE analysis; returns logFC, adjusted-p and t matrices.

    ``stratified`` keeps the original group sizes in every replicate, which is
    the appropriate resampling scheme when comparing two conditions.
    """
    groups = groups.reindex(log_expr.columns)
    rng = np.random.default_rng(seed)
    idx_a = np.flatnonzero((groups == group_a).to_numpy())
    idx_b = np.flatnonzero((groups != group_a).to_numpy())
    y = log_expr.to_numpy(dtype=float)
    n_genes = y.shape[0]
    logfc = np.full((n_genes, n_boot), np.nan)
    padj = np.full((n_genes, n_boot), np.nan)
    tstat = np.full((n_genes, n_boot), np.nan)
    seeds: list[int] = []
    for b in range(n_boot):
        if stratified:
            take_a = rng.choice(idx_a, size=idx_a.size, replace=True)
            take_b = rng.choice(idx_b, size=idx_b.size, replace=True)
        else:
            take = rng.choice(y.shape[1], size=y.shape[1], replace=True)
            g_take = groups.iloc[take].to_numpy()
            take_a = take[g_take == group_a]
            take_b = take[g_take != group_a]
        cols = np.concatenate([take_a, take_b])
        n_a = take_a.size
        sub = pd.DataFrame(y[:, cols], index=log_expr.index, columns=np.arange(cols.size))
        if method == "ebayes":
            design = pd.DataFrame(
                {"Intercept": 1.0, "group": np.r_[np.ones(n_a), np.zeros(cols.size - n_a)]},
                index=sub.columns,
            )
            res = fit_ebayes(sub, design, "group", trend=True, robust=True)
            logfc[:, b] = res.table["logFC"].to_numpy()
            padj[:, b] = res.table["adj.P.Val"].to_numpy()
            tstat[:, b] = res.table["t"].to_numpy()
        else:  # fast Welch fallback
            a = sub.to_numpy()[:, :n_a]
            bb = sub.to_numpy()[:, n_a:]
            t, p = stats.ttest_ind(a, bb, axis=1, equal_var=False)
            logfc[:, b] = a.mean(axis=1) - bb.mean(axis=1)
            padj[:, b] = bh_adjust(p)
            tstat[:, b] = t
        seeds.append(int(rng.integers(0, 2**31 - 1)))
    return BootstrapResult(
        logfc=pd.DataFrame(logfc, index=log_expr.index),
        padj=pd.DataFrame(padj, index=log_expr.index),
        n_ok=n_boot,
        seeds=seeds,
        tstat=pd.DataFrame(tstat, index=log_expr.index),
    )


def gene_stability(
    boot: BootstrapResult,
    base_logfc: pd.Series,
    alpha: float = 0.05,
    min_abs_logfc: float = 1.0,
    weights: tuple[float, float] = (0.7, 0.3),
) -> pd.DataFrame:
    """Per-gene stability components from bootstrap replicates.

    Columns
    -------
    sel_freq          share of replicates that stay significant in the same direction
    sign_consistency  share of replicates whose logFC keeps the observed sign
    effect_cv         coefficient of variation of the bootstrap logFC
    stability         ``sel_freq * (w1*sign_consistency + w2*cv_score)`` with
                      default weights (0.7, 0.3); it is 0 when a gene is never
                      re-selected, so the score separates stable from unstable
                      calls instead of compressing everything into [0.3, 0.6].
    """
    lf = boot.logfc.to_numpy()
    pa = boot.padj.to_numpy()
    sign = np.sign(base_logfc.reindex(boot.logfc.index).fillna(0.0).to_numpy())
    sig = (pa < alpha) & (np.abs(lf) >= min_abs_logfc)
    agree = np.sign(lf) == sign[:, None]
    sel_freq = (sig & agree).mean(axis=1)
    majority = np.where(sign[:, None] >= 0, lf >= 0, lf < 0).mean(axis=1)
    mu = np.nanmean(lf, axis=1)
    sd = np.nanstd(lf, axis=1, ddof=1)
    cv = np.where(np.abs(mu) > 0, sd / np.abs(mu), np.inf)
    cv_score = np.clip(1.0 - cv, 0.0, 1.0)
    w1, w2 = weights
    confidence = w1 * majority + w2 * cv_score
    stability = sel_freq * confidence
    return pd.DataFrame(
        {
            "sel_freq": sel_freq,
            "sign_consistency": majority,
            "effect_cv": cv,
            "effect_sd": sd,
            "logFC_boot_mean": mu,
            "logFC_boot_q05": np.nanquantile(lf, 0.05, axis=1),
            "logFC_boot_q95": np.nanquantile(lf, 0.95, axis=1),
            "cv_score": cv_score,
            "confidence": confidence,
            "stability": stability,
            "n_boot": boot.n_ok,
        },
        index=boot.logfc.index,
    )


def permutation_de(
    log_expr: pd.DataFrame,
    groups: pd.Series,
    group_a: str,
    n_perm: int = 200,
    seed: int = 7,
) -> pd.DataFrame:
    """Group-label permutation null for the (Welch) t-statistic per gene."""
    rng = np.random.default_rng(seed)
    y = log_expr.to_numpy(dtype=float)
    is_a = (groups.reindex(log_expr.columns) == group_a).to_numpy()
    n_a = int(is_a.sum())
    null_t = np.empty((y.shape[0], n_perm))
    for k in range(n_perm):
        perm = np.zeros(y.shape[1], dtype=bool)
        perm[rng.choice(y.shape[1], size=n_a, replace=False)] = True
        t, _ = stats.ttest_ind(y[:, perm], y[:, ~perm], axis=1, equal_var=False)
        null_t[:, k] = t
    return pd.DataFrame(null_t, index=log_expr.index)


def _smooth_kde(eval_at: np.ndarray, samples: np.ndarray) -> np.ndarray:
    """Gaussian KDE with Silverman bandwidth evaluated at ``eval_at``."""
    n = samples.size
    if n < 5:
        return np.full(eval_at.shape, np.nan)
    sd = float(np.std(samples))
    iqr = float(np.subtract(*np.quantile(samples, [0.75, 0.25])))
    scale = min(sd, iqr / 1.34) if iqr > 0 else sd
    bw = 0.9 * scale * n ** (-0.2) if scale > 0 else 0.1
    try:
        return stats.gaussian_kde(samples, bw_method=bw / max(np.std(samples), 1e-9))(eval_at)
    except Exception:  # noqa: BLE001 - singular covariance on degenerate input
        return np.full(eval_at.shape, np.nan)


def empirical_null(t: np.ndarray, n_bins: int = 60) -> pd.DataFrame:
    """Efron-style empirical null and local FDR for a t-statistic vector.

    The null centre is the median of the central 90% of the statistics; the null
    scale is obtained by scaling the robust IQR-based standard deviation of the
    central bulk by the inverse normal quantile of the assumed central mass
    (0.95 by default under a two-sided alternative). Local FDRs are
    ``lfdr_i = pi0 * f0(z_i) / f(z_i)`` with a Gaussian KDE for ``f``.
    """
    t = np.asarray(t, dtype=float)
    ok = np.isfinite(t)
    z = t[ok]
    out = pd.DataFrame({"z": t, "lfdr": np.nan, "pi0": np.nan, "null_mean": np.nan, "null_sd": np.nan})
    if z.size < 100:
        return out
    q90 = float(np.quantile(np.abs(z), 0.90))
    central = z[np.abs(z) <= q90]
    mu0 = float(np.median(central))
    # robust IQR of the central bulk, rescaled to a standard deviation
    q25, q75 = np.quantile(central, [0.25, 0.75])
    sigma_rob = (q75 - q25) / 1.349
    if not np.isfinite(sigma_rob) or sigma_rob <= 1e-6:
        sigma_rob = float(np.std(central))
    # the central mass is 0.90 of the data, i.e. |z| < q90 under the null
    mass = 0.90
    null_sd = float((q90 - mu0) / stats.norm.ppf(0.5 + mass / 2)) if q90 > abs(mu0) else sigma_rob
    null_sd = null_sd if np.isfinite(null_sd) and null_sd > 1e-6 else sigma_rob
    f_emp = _smooth_kde(z, z)
    f0 = stats.norm.pdf(z, mu0, null_sd)
    ratio = np.clip(f_emp / np.maximum(f0, 1e-12), 0, 5)
    keep = np.abs(z - mu0) <= np.quantile(np.abs(z - mu0), 0.5)
    pi0 = float(np.clip(np.median(ratio[keep]), 0.05, 1.0))
    lfdr = np.clip(pi0 * f0 / np.maximum(f_emp, 1e-12), 0.0, 1.0)
    out.loc[ok, "lfdr"] = lfdr
    out.loc[ok, "pi0"] = pi0
    out.loc[ok, "null_mean"] = mu0
    out.loc[ok, "null_sd"] = null_sd
    return out


def consensus_table(
    eb: pd.DataFrame,
    welch: pd.DataFrame,
    wilcoxon: pd.DataFrame,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Merge three DE results and quantify per-gene method concordance."""
    idx = eb.index.intersection(welch.index).intersection(wilcoxon.index)
    out = pd.DataFrame(index=idx)
    for col in ("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"):
        out[col] = eb.loc[idx, col]
    out["P.welch"] = welch.loc[idx, "P.Value"]
    out["P.wilcoxon"] = wilcoxon.loc[idx, "P.Value"]
    out["adj.welch"] = bh_adjust(out["P.welch"].to_numpy())
    out["adj.wilcoxon"] = bh_adjust(out["P.wilcoxon"].to_numpy())
    out["logFC.welch"] = welch.loc[idx, "logFC"]
    out["logFC.wilcoxon"] = wilcoxon.loc[idx, "logFC"]
    out["n_methods_sig"] = (
        (out["adj.P.Val"] < alpha).astype(int)
        + (out["adj.welch"] < alpha).astype(int)
        + (out["adj.wilcoxon"] < alpha).astype(int)
    )
    ranks = pd.DataFrame(
        {
            "eb": out["t"].rank(pct=True),
            "welch": out["P.welch"].rank(pct=True, ascending=False),
            "wilcoxon": out["P.wilcoxon"].rank(pct=True, ascending=False),
        }
    )
    out["rank_concordance"] = ranks.mean(axis=1)
    return out
