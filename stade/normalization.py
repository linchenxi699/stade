"""Count filtering, normalisation and expression transformations.

All functions are pure ``numpy``/``pandas`` and operate on gene x sample
matrices (genes in rows), matching the limma/edgeR convention.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "filter_low_expression",
    "cpm",
    "tmm_factors",
    "log_cpm",
    "quantile_normalize",
    "detect_outliers",
]


def filter_low_expression(
    counts: pd.DataFrame,
    min_count: float = 10.0,
    min_prop: float = 0.25,
    groups: pd.Series | None = None,
) -> pd.DataFrame:
    """Drop genes that are not expressed in a meaningful fraction of samples.

    A gene is kept when, within at least one group (or overall when ``groups``
    is ``None``), it reaches ``min_count`` raw counts in at least ``min_prop``
    of the samples. This mirrors the intent of ``edgeR::filterByExpr`` without
    depending on R.
    """
    if groups is None:
        keep = (counts >= min_count).sum(axis=1) >= max(1, int(np.ceil(min_prop * counts.shape[1])))
    else:
        groups = groups.reindex(counts.columns)
        keep = pd.Series(False, index=counts.index)
        for _, cols in groups.groupby(groups):
            sub = counts.loc[:, cols.index]
            thr = max(1, int(np.ceil(min_prop * sub.shape[1])))
            keep |= (sub >= min_count).sum(axis=1) >= thr
    return counts.loc[keep]


def cpm(counts: pd.DataFrame, log: bool = False, prior_count: float = 0.5) -> pd.DataFrame:
    """Counts-per-million with an optional ``log2`` transform."""
    lib = counts.sum(axis=0).to_numpy(dtype=float)
    lib[lib <= 0] = np.nan
    out = counts.div(lib, axis=1) * 1e6
    if log:
        out = np.log2(out + prior_count)
    return out


def tmm_factors(
    counts: pd.DataFrame,
    ref_col: str | None = None,
    logratio_trim: float = 0.3,
    sum_trim: float = 0.05,
    prior: float = 1.0,
) -> pd.Series:
    """Trimmed mean of M-values (TMM) normalisation factors (edgeR-style).

    Factors are rescaled so their geometric mean is one, hence only relative
    library sizes are changed and expression levels stay interpretable.
    """
    arr = counts.to_numpy(dtype=float)
    lib = arr.sum(axis=0)
    lib[lib <= 0] = np.nan
    lib_ratio = lib / np.exp(np.mean(np.log(lib)))
    if ref_col is None:
        ref_j = int(np.argmin(np.abs(lib_ratio - 1.0)))
    else:
        ref_j = counts.columns.get_loc(ref_col)
    ref = arr[:, ref_j]
    ref_eff = lib[ref_j]
    factors = np.ones(arr.shape[1])
    for j in range(arr.shape[1]):
        if j == ref_j:
            continue
        obs = arr[:, j]
        eff = lib[j]
        ok = (obs > 0) & (ref > 0) & np.isfinite(obs) & np.isfinite(ref)
        if ok.sum() < 100:
            factors[j] = np.nan
            continue
        m = np.log2((obs[ok] / eff) / (ref[ok] / ref_eff))
        a = 0.5 * np.log2((obs[ok] / eff) * (ref[ok] / ref_eff))
        w = 1.0 / (1.0 / (obs[ok] / eff + prior) + 1.0 / (ref[ok] / ref_eff + prior))
        lo_m, hi_m = np.quantile(m, [logratio_trim / 2, 1 - logratio_trim / 2])
        lo_a, hi_a = np.quantile(a, [sum_trim / 2, 1 - sum_trim / 2])
        sel = (m >= lo_m) & (m <= hi_m) & (a >= lo_a) & (a <= hi_a)
        if sel.sum() < 50:
            factors[j] = 1.0
            continue
        factors[j] = float(2 ** (np.sum(w[sel] * m[sel]) / np.sum(w[sel])))
    f = pd.Series(factors, index=counts.columns)
    return f / np.exp(np.mean(np.log(f)))


def log_cpm(counts: pd.DataFrame, prior_count: float = 0.5, use_tmm: bool = True) -> pd.DataFrame:
    """``log2`` CPM after optional TMM normalisation (the limma-voom scale)."""
    norm = counts.div(tmm_factors(counts), axis=1) if use_tmm else counts
    return cpm(norm, log=True, prior_count=prior_count)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalise columns of a DataFrame (ties keep their sample order)."""
    arr = df.to_numpy(dtype=float)
    order = np.argsort(arr, axis=0, kind="stable")
    ranked = np.empty_like(arr)
    for j in range(arr.shape[1]):
        ranked[order[:, j], j] = np.arange(1, arr.shape[0] + 1)
    long = pd.DataFrame({"r": ranked.ravel(), "v": arr.ravel()})
    mean_by_rank = long.groupby("r")["v"].mean().to_numpy()
    out = np.empty_like(arr)
    for j in range(arr.shape[1]):
        out[:, j] = mean_by_rank[ranked[:, j].astype(int) - 1]
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def detect_outliers(expr: pd.DataFrame, z_threshold: float = 3.5) -> pd.Series:
    """Flag samples whose library size is a robust (MAD) outlier."""
    lib = expr.sum(axis=0).replace(0, np.nan)
    loglib = np.log2(lib)
    med = loglib.median()
    mad = np.median(np.abs(loglib - med)) * 1.4826
    mad = mad if mad > 0 else 1e-9
    return ((loglib - med) / mad).abs() > z_threshold
